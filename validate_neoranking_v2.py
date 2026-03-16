"""
NeoRanking Validation v2 — with patient-matched HLA alleles.

Uses the allele column from NeoRanking neopeptide data to test
each peptide against the CORRECT allele. This is what matters.
"""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import warnings
warnings.filterwarnings("ignore")


def run():
    data_file = Path(__file__).parent / "data" / "Neopep_data_org.txt"
    if not data_file.exists():
        print("ERROR: Neopep_data_org.txt not found.")
        sys.exit(1)

    from scorer import score_peptide

    with open(data_file, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    test = [r for r in rows if r.get('train_test') == 'test'
            and r['response_type'] in ('CD8', 'negative')]

    pos_total = sum(1 for r in test if r['response_type'] == 'CD8')
    neg_total = len(test) - pos_total
    print(f"═══ NEORANKING v2 — PATIENT-MATCHED ALLELES ═══")
    print(f"Test set: {len(test)} peptides ({pos_total} CD8+, {neg_total} negative)\n")

    scores_pos = []
    scores_neg = []
    processed = 0

    for r in test:
        peptide = r.get('mutant_seq', '').strip()
        wildtype = r.get('wt_seq', '').strip()

        if not peptide or len(peptide) < 8 or len(peptide) > 11:
            continue

        if not wildtype or len(wildtype) != len(peptide):
            wildtype = peptide

        # Patient-matched allele from NeoRanking
        allele_raw = r.get('mutant_best_alleles', '').strip()
        if allele_raw:
            # Normalize: B0801 → HLA-B*08:01
            if len(allele_raw) >= 4 and '*' not in allele_raw:
                gene = allele_raw[0]
                digits = allele_raw[1:]
                if len(digits) == 4:
                    allele = f"HLA-{gene}*{digits[:2]}:{digits[2:]}"
                else:
                    allele = f"HLA-{allele_raw}"
            else:
                allele = allele_raw if allele_raw.startswith('HLA-') else f"HLA-{allele_raw}"
            alleles = [allele]
        else:
            alleles = None  # fall back to defaults

        actual = 1 if r['response_type'] == 'CD8' else 0

        score = score_peptide(peptide, wildtype, alleles=alleles)

        if actual == 1:
            scores_pos.append(score)
        else:
            scores_neg.append(score)

        processed += 1
        if processed % 5000 == 0:
            print(f"  processed {processed}/{len(test)}...", flush=True)

    # AUC
    if scores_pos and scores_neg:
        auc = sum(1 for p in scores_pos for n in scores_neg if p > n) / (len(scores_pos) * len(scores_neg))
    else:
        auc = 0.5

    n_pos, n_neg = len(scores_pos), len(scores_neg)
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    se = math.sqrt((auc*(1-auc) + (n_pos-1)*(q1-auc**2) + (n_neg-1)*(q2-auc**2)) / (n_pos * n_neg))

    avg_pos = sum(scores_pos) / len(scores_pos) if scores_pos else 0
    avg_neg = sum(scores_neg) / len(scores_neg) if scores_neg else 0

    # Count gated (score = 0.05)
    gated_pos = sum(1 for s in scores_pos if s <= 0.05)
    gated_neg = sum(1 for s in scores_neg if s <= 0.05)

    print(f"\nScored: {n_pos} positive, {n_neg} negative")
    print(f"Gated (score ≤ 0.05): {gated_pos}/{n_pos} positive ({gated_pos/n_pos*100:.0f}%), {gated_neg}/{n_neg} negative ({gated_neg/n_neg*100:.0f}%)")
    print(f"ROC-AUC: {auc:.4f}  95% CI [{auc-1.96*se:.4f}, {auc+1.96*se:.4f}]")
    print(f"Avg immunogenic: {avg_pos:.4f}  Avg negative: {avg_neg:.4f}  Separation: {avg_pos-avg_neg:.4f}")

    if auc > 0.7:
        print(f"\n✓ SIGNAL on {n_pos+n_neg} peptides (AUC={auc:.3f})")
    elif auc > 0.6:
        print(f"\n~ WEAK SIGNAL on {n_pos+n_neg} peptides")
    else:
        print(f"\n✗ NO SIGNAL on {n_pos+n_neg} peptides")

    return {"auc": auc, "n": n_pos + n_neg, "gated_pos_pct": gated_pos/n_pos*100}


if __name__ == "__main__":
    run()
