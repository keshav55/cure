"""
NeoRanking Validation — 4,307 mutations from 131 patients.

The definitive benchmark. Real wildtype sequences, real expression data,
proper train/test split, 106 CD8+ immunogenic mutations vs 4,201 negatives.

Source: Muller et al. 2023 Immunity, NCI + TESLA + HiTIDE harmonized dataset.
"""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import warnings
warnings.filterwarnings("ignore")


def run():
    data_file = Path(__file__).parent / "data" / "Mutation_data_org.txt"
    if not data_file.exists():
        print("ERROR: Mutation_data_org.txt not found.")
        sys.exit(1)

    from scorer import score_peptide

    with open(data_file, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    # Test set only — no training data contamination
    test = [r for r in rows if r.get('train_test') == 'test'
            and r['response_type'] in ('CD8', 'negative')]

    print(f"═══ NEORANKING VALIDATION — {len(test)} mutations ═══")
    pos_total = sum(1 for r in test if r['response_type'] == 'CD8')
    neg_total = len(test) - pos_total
    print(f"Positive (CD8): {pos_total}  Negative: {neg_total}\n")

    scores_pos = []
    scores_neg = []

    for r in test:
        mutant_seq = r.get('mutant_seq', '').strip()
        wt_seq = r.get('wt_seq', '').strip()

        if not mutant_seq or len(mutant_seq) < 8:
            continue

        # Extract 9-mer centered on mutation from the 25-mer context
        # The mutation is at the center of the 25-mer
        center = len(mutant_seq) // 2
        for length in [9, 10, 8]:
            start = max(0, center - length // 2)
            end = start + length
            if end <= len(mutant_seq) and start >= 0:
                peptide = mutant_seq[start:end]
                wildtype = wt_seq[start:end] if len(wt_seq) >= end else peptide
                if peptide != wildtype:
                    break
        else:
            peptide = mutant_seq[center-4:center+5] if len(mutant_seq) >= center+5 else mutant_seq[:9]
            wildtype = wt_seq[center-4:center+5] if len(wt_seq) >= center+5 else wt_seq[:9]

        if len(peptide) < 8:
            continue

        actual = 1 if r['response_type'] == 'CD8' else 0
        our_score = score_peptide(peptide, wildtype)

        if actual == 1:
            scores_pos.append(our_score)
        else:
            scores_neg.append(our_score)

    # AUC
    if scores_pos and scores_neg:
        auc = sum(1 for p in scores_pos for n in scores_neg if p > n) / (len(scores_pos) * len(scores_neg))
    else:
        auc = 0.5

    # CI
    n_pos, n_neg = len(scores_pos), len(scores_neg)
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    se = math.sqrt((auc * (1-auc) + (n_pos-1)*(q1-auc**2) + (n_neg-1)*(q2-auc**2)) / (n_pos * n_neg))

    avg_pos = sum(scores_pos) / len(scores_pos) if scores_pos else 0
    avg_neg = sum(scores_neg) / len(scores_neg) if scores_neg else 0

    print(f"Scored: {len(scores_pos)} positive, {len(scores_neg)} negative")
    print(f"ROC-AUC: {auc:.4f}  95% CI [{auc-1.96*se:.4f}, {auc+1.96*se:.4f}]")
    print(f"Avg immunogenic: {avg_pos:.4f}  Avg negative: {avg_neg:.4f}  Separation: {avg_pos-avg_neg:.4f}")

    if auc > 0.7:
        print(f"\n✓ SIGNAL on {n_pos+n_neg} NeoRanking mutations (AUC={auc:.3f})")
    elif auc > 0.6:
        print(f"\n~ WEAK SIGNAL on {n_pos+n_neg} NeoRanking mutations")
    else:
        print(f"\n✗ NO SIGNAL on {n_pos+n_neg} NeoRanking mutations")

    return {"auc": auc, "ci_low": auc-1.96*se, "ci_high": auc+1.96*se, "n": n_pos+n_neg}


if __name__ == "__main__":
    run()
