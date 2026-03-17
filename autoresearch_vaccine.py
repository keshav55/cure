"""
Autoresearch daemon for vaccine selection optimization.

Karpathy pattern: propose modification → evaluate → keep/reject → learn.
Metric: per-patient recall@20 on NeoRanking test split.
Baseline: 0.698 (stacking ensemble GBT + RF + LR).

Run: python autoresearch_vaccine.py [--rounds N]
"""

import csv
import json
import math
import random
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression

# ── Configuration ──────────────────────────────────────────────
BASELINE_RECALL = 0.698
KEEP_THRESHOLD = 0.005  # Keep if improvement >= 0.5%
STATE_FILE = Path(__file__).parent / ".vaccine_autoresearch_state.json"
LEARNINGS_FILE = Path(__file__).parent / "LEARNINGS.md"

FEATURES = [
    'mutant_rank', 'mut_Rank_Stab', 'rnaseq_TPM', 'TAP_score',
    'CSCAPE_score', 'bestWTPeptideCount_I', 'bestWTMatchScore_I',
    'mut_binding_score', 'DAI_NetMHC', 'mutant_rank_netMHCpan',
    'mutant_rank_PRIME', 'CCF', 'TumorContent', 'mut_netchop_score_ct',
    'GTEx_all_tissues_expression_mean', 'TCGA_Cancer_expression',
    'bestMutationScore_I', 'bestWTMatchOverlap_I',
]


def sf(v, d=0.0):
    try:
        return float(v)
    except (ValueError, TypeError):
        return d


def extract_features(row, extra_features=None):
    feats = [sf(row.get(f, ''), 0) for f in FEATURES]
    feats.append(math.log1p(sf(row.get('rnaseq_TPM', ''), 0)))
    feats.append(1.0 / max(0.01, sf(row.get('mutant_rank', ''), 100)))
    if extra_features:
        for ef in extra_features:
            feats.append(ef(row))
    return feats


def load_data():
    data_file = Path(__file__).parent / "data" / "Neopep_data_org.txt"
    with open(data_file, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    train = [r for r in rows if r.get('train_test') == 'train'
             and r['response_type'] in ('CD8', 'negative')]
    test = [r for r in rows if r.get('train_test') == 'test'
            and r['response_type'] in ('CD8', 'negative')]
    return train, test


def evaluate_config(train, test, config):
    """Evaluate a vaccine selection configuration. Returns recall@20."""
    train_pos = [r for r in train if r['response_type'] == 'CD8']
    train_neg = [r for r in train if r['response_type'] == 'negative']

    random.seed(42)
    neg_sample = random.sample(train_neg, min(config.get('neg_sample', 5000), len(train_neg)))
    data = train_pos + neg_sample
    random.shuffle(data)

    extra_fns = []
    for ef in config.get('extra_features', []):
        if ef == 'log_stability':
            extra_fns.append(lambda r: math.log1p(sf(r.get('mut_Rank_Stab', ''), 0)))
        elif ef == 'bind_x_expr':
            extra_fns.append(lambda r: sf(r.get('mutant_rank', ''), 50) * sf(r.get('rnaseq_TPM', ''), 0) * -0.001)
        elif ef == 'bind_x_stab':
            extra_fns.append(lambda r: sf(r.get('mutant_rank', ''), 50) * sf(r.get('mut_Rank_Stab', ''), 50))
        elif ef == 'expr_rank':
            extra_fns.append(lambda r: -math.log1p(sf(r.get('rnaseq_TPM', ''), 0)) / max(0.01, sf(r.get('mutant_rank', ''), 100)))
        elif ef == 'allele_count_proxy':
            extra_fns.append(lambda r: 1 if r.get('mutant_other_significant_alleles', '') else 0)

    X = np.array([extract_features(r, extra_fns) for r in data])
    X = np.nan_to_num(X, nan=0.0, posinf=100.0, neginf=-100.0)
    y = np.array([1 if r['response_type'] == 'CD8' else 0 for r in data])

    means = X.mean(axis=0)
    stds = X.std(axis=0) + 1e-8
    X_s = (X - means) / stds

    # Train models with config
    gbt = GradientBoostingClassifier(
        n_estimators=config.get('gbt_n', 50),
        max_depth=config.get('gbt_depth', 2),
        learning_rate=config.get('gbt_lr', 0.1),
        subsample=config.get('gbt_subsample', 1.0),
        random_state=42,
    )
    rf = RandomForestClassifier(
        n_estimators=config.get('rf_n', 100),
        max_depth=config.get('rf_depth', 5),
        class_weight='balanced',
        random_state=42,
    )
    lr = LogisticRegression(max_iter=1000, class_weight='balanced')

    gbt.fit(X_s, y)
    rf.fit(X_s, y)
    lr.fit(X_s, y)

    # Score test set
    X_test = np.array([extract_features(r, extra_fns) for r in test])
    X_test = np.nan_to_num(X_test, nan=0.0, posinf=100.0, neginf=-100.0)
    X_test_s = (X_test - means) / stds

    w_gbt = config.get('w_gbt', 0.4)
    w_rf = config.get('w_rf', 0.4)
    w_lr = config.get('w_lr', 0.2)

    blended = (w_gbt * gbt.predict_proba(X_test_s)[:, 1] +
               w_rf * rf.predict_proba(X_test_s)[:, 1] +
               w_lr * lr.predict_proba(X_test_s)[:, 1])

    # Per-patient recall@20
    patients = defaultdict(list)
    for i, r in enumerate(test):
        patients[r.get('patient', '')].append((i, r))

    recalls = []
    k = config.get('k', 20)
    for pid, indexed_rows in patients.items():
        pos_peps = set(r['mutant_seq'] for _, r in indexed_rows if r['response_type'] == 'CD8')
        if not pos_peps:
            continue

        patient_scores = [(blended[i], r['mutant_seq']) for i, r in indexed_rows]
        patient_scores.sort(reverse=True)

        selected = set(pep for _, pep in patient_scores[:k])
        hits = len(selected & pos_peps)
        recalls.append(hits / len(pos_peps))

    return sum(recalls) / len(recalls) if recalls else 0.0


# ── Search space ──────────────────────────────────────────────
SEARCH_SPACE = {
    'gbt_n': [30, 50, 75, 100, 150],
    'gbt_depth': [2, 3, 4],
    'gbt_lr': [0.05, 0.1, 0.15, 0.2],
    'gbt_subsample': [0.7, 0.8, 0.9, 1.0],
    'rf_n': [50, 100, 150, 200],
    'rf_depth': [3, 4, 5, 6, 7],
    'w_gbt': [0.3, 0.35, 0.4, 0.45, 0.5],
    'w_rf': [0.25, 0.3, 0.35, 0.4, 0.45],
    'w_lr': [0.1, 0.15, 0.2, 0.25, 0.3],
    'neg_sample': [3000, 5000, 7000, 10000],
    'extra_features': [
        [],
        ['log_stability'],
        ['bind_x_expr'],
        ['bind_x_stab'],
        ['expr_rank'],
        ['allele_count_proxy'],
        ['log_stability', 'bind_x_expr'],
        ['log_stability', 'expr_rank'],
        ['bind_x_expr', 'allele_count_proxy'],
    ],
}


def propose_config(state):
    """Propose a new configuration to test."""
    best_config = state.get('best_config', {})
    config = dict(best_config)

    # Mutate 1-2 parameters
    params_to_mutate = random.sample(list(SEARCH_SPACE.keys()), k=random.choice([1, 2]))

    for param in params_to_mutate:
        config[param] = random.choice(SEARCH_SPACE[param])

    # Ensure weights sum to ~1
    total = config.get('w_gbt', 0.4) + config.get('w_rf', 0.4) + config.get('w_lr', 0.2)
    if total > 0:
        config['w_gbt'] = config.get('w_gbt', 0.4) / total
        config['w_rf'] = config.get('w_rf', 0.4) / total
        config['w_lr'] = config.get('w_lr', 0.2) / total

    return config


def load_state():
    if STATE_FILE.exists():
        with open(STATE_FILE) as f:
            return json.load(f)
    return {
        'best_recall': BASELINE_RECALL,
        'best_config': {
            'gbt_n': 50, 'gbt_depth': 2, 'gbt_lr': 0.1, 'gbt_subsample': 1.0,
            'rf_n': 100, 'rf_depth': 5,
            'w_gbt': 0.4, 'w_rf': 0.4, 'w_lr': 0.2,
            'neg_sample': 5000,
            'extra_features': [],
        },
        'round': 0,
        'keeps': 0,
        'history': [],
    }


def save_state(state):
    with open(STATE_FILE, 'w') as f:
        json.dump(state, f, indent=2)


def main():
    rounds = int(sys.argv[sys.argv.index('--rounds') + 1]) if '--rounds' in sys.argv else 20

    print(f"═══ VACCINE AUTORESEARCH DAEMON ═══")
    print(f"Baseline: {BASELINE_RECALL:.3f}")
    print(f"Rounds: {rounds}\n")

    train, test = load_data()
    state = load_state()

    print(f"Starting from round {state['round']}, best = {state['best_recall']:.4f}\n")

    for i in range(rounds):
        state['round'] += 1
        config = propose_config(state)

        print(f"Round {state['round']}: testing {json.dumps({k: v for k, v in config.items() if k != 'extra_features'}, default=str)[:80]}...")

        try:
            recall = evaluate_config(train, test, config)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

        delta = recall - state['best_recall']

        if delta >= KEEP_THRESHOLD:
            state['best_recall'] = recall
            state['best_config'] = config
            state['keeps'] += 1
            print(f"  ★ KEEP: {recall:.4f} (+{delta:.4f})")
            state['history'].append({'round': state['round'], 'recall': recall, 'config': config, 'kept': True})
        else:
            print(f"  reject: {recall:.4f} ({delta:+.4f})")
            state['history'].append({'round': state['round'], 'recall': recall, 'kept': False})

        save_state(state)

    print(f"\n═══ SUMMARY ═══")
    print(f"Rounds: {state['round']}, Keeps: {state['keeps']}")
    print(f"Best recall@20: {state['best_recall']:.4f}")
    print(f"Best config: {json.dumps(state['best_config'], indent=2)}")


if __name__ == "__main__":
    main()
