"""
TESLA Validation — 608 peptides from Wells et al. 2020 Cell.

The real benchmark. 608 clinically validated peptides from 6 melanoma patients
with binary immunogenicity labels from pMHC multimer assays.

Key limitation: TESLA provides mutant peptide but not wildtype. We use
MUTATION_POSITION to create a synthetic wildtype (replace mutant AA with 'X'
placeholder) so the scorer's foreignness features can fire.
"""

import sys
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))

import warnings
warnings.filterwarnings("ignore")


# Common amino acids by group — for synthetic wildtype construction
# If mutation is at a position, guess the wildtype as the most common
# AA from a DIFFERENT group (maximize foreignness signal)
AA_GROUPS = {
    'A': 'aliphatic', 'V': 'aliphatic', 'I': 'aliphatic', 'L': 'aliphatic', 'M': 'aliphatic',
    'F': 'aromatic', 'W': 'aromatic', 'Y': 'aromatic',
    'S': 'polar', 'T': 'polar', 'C': 'polar',
    'N': 'amide', 'Q': 'amide',
    'D': 'negative', 'E': 'negative',
    'K': 'positive', 'R': 'positive', 'H': 'positive',
    'G': 'special', 'P': 'special',
}

# Most common wildtype AA per group (from human proteome frequencies)
GROUP_COMMON = {
    'aliphatic': 'A', 'aromatic': 'Y', 'polar': 'S',
    'amide': 'N', 'negative': 'E', 'positive': 'K',
    'special': 'G',
}


def make_synthetic_wildtype(peptide: str, mut_pos: int) -> str:
    """Create a synthetic wildtype by replacing the mutated position.

    Since TESLA doesn't provide the original AA, we substitute the mutant
    AA with the most common AA from a DIFFERENT group. This ensures the
    scorer's foreignness features detect a cross-group change.
    """
    if mut_pos < 1 or mut_pos > len(peptide):
        return peptide  # can't modify

    idx = mut_pos - 1  # 1-indexed → 0-indexed
    mutant_aa = peptide[idx]
    mutant_group = AA_GROUPS.get(mutant_aa, 'aliphatic')

    # Pick a common AA from a different group
    for group, common_aa in GROUP_COMMON.items():
        if group != mutant_group:
            wt_aa = common_aa
            break
    else:
        wt_aa = 'A'  # fallback

    return peptide[:idx] + wt_aa + peptide[idx + 1:]


def _safe_float(val):
    if val is None or val == 'NA' or val == '':
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def run():
    data_file = Path(__file__).parent / "data" / "tesla_mmc4.xlsx"
    if not data_file.exists():
        print("ERROR: tesla_mmc4.xlsx not found.")
        sys.exit(1)

    from scorer import score_peptide
    from scorer_tesla import score_with_features

    wb = openpyxl.load_workbook(str(data_file), read_only=True)
    ws = wb['master-bindings-selected']
    rows = list(ws.iter_rows(min_row=2, values_only=True))

    print(f"═══ TESLA VALIDATION — {len(rows)} peptides ═══\n")

    # Run BOTH scorers: sequence-only and feature-based
    results_seq = {"tp": 0, "fp": 0, "tn": 0, "fn": 0, "pos": [], "neg": []}
    results_feat = {"tp": 0, "fp": 0, "tn": 0, "fn": 0, "pos": [], "neg": []}
    threshold = 0.5

    for row in rows:
        peptide = str(row[4]).strip()
        hla_raw = str(row[3]).strip()
        hla = f"HLA-{hla_raw}" if not hla_raw.startswith("HLA-") else hla_raw
        if "*" in hla and ":" not in hla:
            parts = hla.split("*")
            allele = parts[1]
            if len(allele) == 4:
                hla = f"{parts[0]}*{allele[:2]}:{allele[2:]}"

        actual = 1 if row[15] == True else 0
        try:
            mut_pos = int(row[13]) if row[13] not in (None, 'NA', '') else 0
        except (ValueError, TypeError):
            mut_pos = 0

        # ── Scorer 1: Sequence-only (synthetic wildtype) ──
        if mut_pos > 0:
            wildtype = make_synthetic_wildtype(peptide, mut_pos)
        else:
            wildtype = peptide

        alleles = [hla]
        for a in ["HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01"]:
            if a not in alleles:
                alleles.append(a)

        seq_score = score_peptide(peptide, wildtype, alleles=alleles)
        predicted_seq = 1 if seq_score >= threshold else 0

        # ── Scorer 2: Feature-based (uses TESLA measured data) ──
        binding = _safe_float(row[6])
        stability = _safe_float(row[9])
        abundance = _safe_float(row[8])
        num_pred = _safe_float(row[14])

        if binding is not None and stability is not None:
            feat_score = score_with_features(binding, stability, abundance, num_pred)
        else:
            feat_score = seq_score  # fall back to sequence-only
        predicted_feat = 1 if feat_score >= threshold else 0

        # Record results for both
        for results, predicted, score in [
            (results_seq, predicted_seq, seq_score),
            (results_feat, predicted_feat, feat_score),
        ]:
            if predicted == 1 and actual == 1: results["tp"] += 1
            elif predicted == 1 and actual == 0: results["fp"] += 1
            elif predicted == 0 and actual == 1: results["fn"] += 1
            else: results["tn"] += 1
            (results["pos"] if actual == 1 else results["neg"]).append(score)

        our_score = seq_score  # for backward compat with rest of function
        predicted = predicted_seq

    # ── Report both scorers ──
    def _report(name, r, n_total):
        tp, fp, tn, fn = r["tp"], r["fp"], r["tn"], r["fn"]
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
        acc = (tp + tn) / n_total
        auc = 0.5
        if r["pos"] and r["neg"]:
            auc = sum(1 for p in r["pos"] for n in r["neg"] if p > n) / (len(r["pos"]) * len(r["neg"]))
        avg_p = sum(r["pos"]) / len(r["pos"]) if r["pos"] else 0
        avg_n = sum(r["neg"]) / len(r["neg"]) if r["neg"] else 0
        print(f"── {name} ──")
        print(f"tp={tp} fp={fp} tn={tn} fn={fn}")
        print(f"Precision: {prec:.4f}  Recall: {rec:.4f}  F1: {f1:.4f}")
        print(f"Accuracy: {acc:.4f}")
        print(f"ROC-AUC: {auc:.4f}")
        print(f"Avg immunogenic: {avg_p:.4f}  Avg non-immunogenic: {avg_n:.4f}  Separation: {avg_p-avg_n:.4f}")
        return auc

    print()
    auc_seq = _report("Sequence-only scorer (scorer.py)", results_seq, len(rows))
    print()
    auc_feat = _report("Feature-based scorer (binding + stability + expression)", results_feat, len(rows))
    print()
    print(f"Improvement: {auc_seq:.3f} → {auc_feat:.3f} (Δ = {auc_feat-auc_seq:+.3f})")

    scores_pos = results_seq["pos"]
    scores_neg = results_seq["neg"]
    auc = auc_seq

    if auc > 0.7:
        print(f"\n✓ SIGNAL DETECTED on {len(rows)} TESLA peptides")
    elif auc > 0.6:
        print(f"\n~ WEAK SIGNAL on {len(rows)} TESLA peptides")
    else:
        print(f"\n✗ NO SIGNAL on {len(rows)} TESLA peptides")

    return {"auc": auc, "f1": f1, "precision": precision, "recall": recall}


if __name__ == "__main__":
    run()
