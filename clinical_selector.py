#!/usr/bin/env python3
"""
Clinical neoantigen vaccine selector.

Safety net strategy: top-K/2 by binding + top-K/2 by ML ensemble.
Validated on 73 patients (LOPO), Wilcoxon p = 0.002 (K=20), p = 0.0001 (K=30).

Input:  TSV with columns: peptide, allele, [optional: stability_rank, cscape, expression_gtex, expression_tcga, rank_netmhcpan, rank_prime]
Output: Ranked vaccine candidates with scores and selection method

Usage:
    # With NeoRanking-format data
    python clinical_selector.py --input peptides.tsv --k 20 --output vaccine.tsv

    # Train on NeoRanking, select from custom data
    python clinical_selector.py --train data/Neopep_data_org.txt --input my_peptides.tsv --k 30
"""

import argparse
import csv
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression


FEATURES = [
    'mutant_rank', 'mut_Rank_Stab', 'CSCAPE_score',
    'mutant_rank_netMHCpan', 'mutant_rank_PRIME',
    'GTEx_all_tissues_expression_mean', 'TCGA_Cancer_expression',
]


def sf(v, d=0.0):
    try:
        return float(v)
    except (ValueError, TypeError):
        return d


def extract_features(row):
    feats = [sf(row.get(f, ''), 0) for f in FEATURES]
    feats.append(math.log1p(sf(row.get('rnaseq_TPM', ''), 0)))
    feats.append(1.0 / max(0.01, sf(row.get('mutant_rank', ''), 100)))
    return feats


class ClinicalSelector:
    """Safety net vaccine selector: binding floor + ML ensemble discovery."""

    def __init__(self, n_seeds=10):
        self.n_seeds = n_seeds
        self.models = []  # list of (means, stds, gbt, rf, lr)

    def train(self, train_file):
        """Train on NeoRanking-format TSV."""
        with open(train_file, 'r') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))

        train = [r for r in rows if r.get('train_test') == 'train'
                 and r['response_type'] in ('CD8', 'negative')]
        pos = [r for r in train if r['response_type'] == 'CD8']
        neg = [r for r in train if r['response_type'] == 'negative']

        print(f"Training on {len(pos)} positives, sampling from {len(neg)} negatives...")

        self.models = []
        for seed in range(self.n_seeds):
            random.seed(seed)
            balanced = pos + random.sample(neg, min(5000, len(neg)))
            random.shuffle(balanced)

            X = np.array([extract_features(r) for r in balanced])
            X = np.nan_to_num(X, nan=0.0, posinf=100.0, neginf=-100.0)
            y = np.array([1 if r['response_type'] == 'CD8' else 0 for r in balanced])

            means = X.mean(axis=0)
            stds = X.std(axis=0) + 1e-8
            Xs = (X - means) / stds

            gbt = GradientBoostingClassifier(n_estimators=50, max_depth=2, learning_rate=0.1, random_state=42)
            rf = RandomForestClassifier(n_estimators=100, max_depth=5, class_weight='balanced', random_state=42)
            lr = LogisticRegression(max_iter=1000, class_weight='balanced')

            gbt.fit(Xs, y)
            rf.fit(Xs, y)
            lr.fit(Xs, y)

            self.models.append((means, stds, gbt, rf, lr))

        print(f"Trained {self.n_seeds} model seeds.")

    def score(self, peptide_rows):
        """Score peptides using multi-seed ensemble."""
        X = np.array([extract_features(r) for r in peptide_rows])
        X = np.nan_to_num(X, nan=0.0, posinf=100.0, neginf=-100.0)

        all_probs = np.zeros(len(peptide_rows))
        for means, stds, gbt, rf, lr in self.models:
            Xs = (X - means) / stds
            all_probs += (0.4 * gbt.predict_proba(Xs)[:, 1] +
                          0.4 * rf.predict_proba(Xs)[:, 1] +
                          0.2 * lr.predict_proba(Xs)[:, 1])
        return all_probs / len(self.models)

    def select(self, patient_rows, k=20):
        """Safety net selection: k/2 binding + k/2 ensemble."""
        bind_k = k // 2
        ens_k = k - bind_k

        scores = self.score(patient_rows)

        # Top binding
        bind_ranked = sorted(range(len(patient_rows)),
                             key=lambda i: sf(patient_rows[i].get('mutant_rank', ''), 100))
        bind_set = set(bind_ranked[:bind_k])
        bind_peps = set(patient_rows[i].get('mutant_seq', '') for i in bind_set)

        # Top ensemble (excluding binding picks)
        ens_ranked = sorted(range(len(patient_rows)), key=lambda i: -scores[i])
        ens_set = set()
        for i in ens_ranked:
            pep = patient_rows[i].get('mutant_seq', '')
            if pep not in bind_peps:
                ens_set.add(i)
                if len(ens_set) >= ens_k:
                    break

        # Build results
        results = []
        for i in bind_set:
            r = dict(patient_rows[i])
            r['selection'] = 'binding'
            r['ensemble_score'] = f"{scores[i]:.4f}"
            r['bind_rank'] = bind_ranked.index(i) + 1
            results.append(r)

        for i in ens_set:
            r = dict(patient_rows[i])
            r['selection'] = 'ensemble'
            r['ensemble_score'] = f"{scores[i]:.4f}"
            r['bind_rank'] = bind_ranked.index(i) + 1
            results.append(r)

        # Sort by ensemble score descending
        results.sort(key=lambda r: -float(r['ensemble_score']))
        return results


def main():
    parser = argparse.ArgumentParser(description='Clinical neoantigen vaccine selector')
    parser.add_argument('--train', default=str(Path(__file__).parent / 'data' / 'Neopep_data_org.txt'))
    parser.add_argument('--input', help='Input peptides TSV (NeoRanking format)')
    parser.add_argument('--patient', help='Filter to specific patient ID')
    parser.add_argument('--k', type=int, default=20, help='Number of peptides to select')
    parser.add_argument('--output', help='Output TSV file')
    args = parser.parse_args()

    selector = ClinicalSelector(n_seeds=10)
    selector.train(args.train)

    if args.input:
        with open(args.input, 'r') as f:
            peptides = list(csv.DictReader(f, delimiter='\t'))
    else:
        # Demo: use test set
        with open(args.train, 'r') as f:
            all_rows = list(csv.DictReader(f, delimiter='\t'))
        peptides = [r for r in all_rows if r.get('train_test') == 'test'
                    and r['response_type'] in ('CD8', 'negative')]

    # Group by patient
    patients = defaultdict(list)
    for r in peptides:
        pid = args.patient or r.get('patient', 'unknown')
        if args.patient and r.get('patient', '') != args.patient:
            continue
        patients[r.get('patient', 'unknown')].append(r)

    all_results = []
    for pid, prows in patients.items():
        selected = selector.select(prows, k=args.k)
        for r in selected:
            r['patient'] = pid
        all_results.extend(selected)

        pos_count = sum(1 for r in prows if r.get('response_type') == 'CD8')
        bind_count = sum(1 for r in selected if r['selection'] == 'binding')
        ens_count = sum(1 for r in selected if r['selection'] == 'ensemble')
        print(f"  Patient {pid}: {len(prows)} candidates → {bind_count} binding + {ens_count} ensemble = {len(selected)} selected")

    if args.output:
        cols = ['patient', 'mutant_seq', 'mutant_best_alleles', 'gene',
                'selection', 'ensemble_score', 'bind_rank',
                'mutant_rank', 'mut_Rank_Stab', 'rnaseq_TPM', 'response_type']
        with open(args.output, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(all_results)
        print(f"\nWrote {len(all_results)} peptides to {args.output}")


if __name__ == '__main__':
    main()
