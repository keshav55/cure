"""
Cancer Vaccine Pipeline — from tumor mutations to vaccine candidates.

Pi Day 2026. This is the tool that connects:
  1. Tumor mutation data (VCF/MAF or raw mutations)
  2. Peptide candidate generation (sliding window over mutated proteins)
  3. MHC binding prediction (MHCflurry neural network)
  4. Neoantigen immunogenicity scoring
  5. mRNA sequence design (codon-optimized)

Input: list of mutations (protein, position, wildtype_aa, mutant_aa)
Output: ranked vaccine candidates with mRNA sequences

Usage:
    python pipeline.py                    # run demo with example mutations
    python pipeline.py --mutations file   # run with custom mutation file
"""

import warnings
warnings.filterwarnings("ignore")

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Add project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "backend"))


@dataclass
class Mutation:
    """A somatic mutation found in the tumor."""
    gene: str
    protein_position: int
    wildtype_aa: str
    mutant_aa: str
    protein_context: str = ""  # surrounding amino acid context (21+ residues)


@dataclass
class PeptideCandidate:
    """A candidate neoantigen peptide."""
    peptide: str
    wildtype: str
    gene: str
    mutation_position: int
    length: int
    mhc_affinity_nm: float = 0.0
    mhc_presentation: float = 0.0
    immunogenicity_score: float = 0.0
    mrna_sequence: str = ""
    rank: int = 0


# ── Standard genetic code ──
CODON_TABLE = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    '*': ['TAA', 'TAG', 'TGA'],  # stop codons
}

# Human codon usage frequency (per thousand) — for codon optimization
# Source: Kazusa codon usage database, Homo sapiens
HUMAN_CODON_FREQ = {
    'TTT': 17.6, 'TTC': 20.3, 'TTA': 7.7, 'TTG': 12.9,
    'CTT': 13.2, 'CTC': 19.6, 'CTA': 7.2, 'CTG': 39.6,
    'ATT': 16.0, 'ATC': 20.8, 'ATA': 7.5, 'ATG': 22.0,
    'GTT': 11.0, 'GTC': 14.5, 'GTA': 7.1, 'GTG': 28.1,
    'TCT': 15.2, 'TCC': 17.7, 'TCA': 12.2, 'TCG': 4.4,
    'AGT': 12.1, 'AGC': 19.5, 'CCT': 17.5, 'CCC': 19.8,
    'CCA': 16.9, 'CCG': 6.9, 'ACT': 13.1, 'ACC': 18.9,
    'ACA': 15.1, 'ACG': 6.1, 'GCT': 18.4, 'GCC': 27.7,
    'GCA': 15.8, 'GCG': 7.4, 'TAT': 12.2, 'TAC': 15.3,
    'CAT': 10.9, 'CAC': 15.1, 'CAA': 12.3, 'CAG': 34.2,
    'AAT': 17.0, 'AAC': 19.1, 'AAA': 24.4, 'AAG': 31.9,
    'GAT': 21.8, 'GAC': 25.1, 'GAA': 29.0, 'GAG': 39.6,
    'TGT': 10.6, 'TGC': 12.6, 'TGG': 13.2,
    'CGT': 4.5, 'CGC': 10.4, 'CGA': 6.2, 'CGG': 11.4,
    'AGA': 12.2, 'AGG': 12.0, 'GGT': 10.8, 'GGC': 22.2,
    'GGA': 16.5, 'GGG': 16.5,
}


def generate_peptide_candidates(mutation: Mutation, lengths=(8, 9, 10, 11)) -> list[PeptideCandidate]:
    """Generate all possible peptides containing the mutation.

    For each peptide length, slide a window over the protein context
    such that the mutation position falls within the peptide.
    """
    context = mutation.protein_context
    if not context:
        return []

    # Find the mutation position within the context
    # Context was built by resolve_protein_context: window=12 around the mutation
    # But if mutation is near protein start/end, it may not be centered
    # The mutation position in the full protein tells us where it is in the context
    from protein_db import get_protein_sequence
    seq = get_protein_sequence(mutation.gene)
    if seq:
        idx = mutation.protein_position - 1  # 0-indexed in full protein
        start = max(0, idx - 12)
        mut_idx = idx - start  # position within context
    else:
        mut_idx = len(context) // 2  # fallback: assume centered

    candidates = []
    for length in lengths:
        # Slide window: mutation must be within the peptide
        for start in range(max(0, mut_idx - length + 1), min(mut_idx + 1, len(context) - length + 1)):
            end = start + length
            peptide = context[start:end]
            # Build wildtype by substituting back
            wt_context = context[:mut_idx] + mutation.wildtype_aa + context[mut_idx + 1:]
            wildtype = wt_context[start:end]

            if len(peptide) == length and peptide != wildtype:
                candidates.append(PeptideCandidate(
                    peptide=peptide,
                    wildtype=wildtype,
                    gene=mutation.gene,
                    mutation_position=mutation.protein_position,
                    length=length,
                ))

    # Deduplicate
    seen = set()
    unique = []
    for c in candidates:
        key = (c.peptide, c.wildtype)
        if key not in seen:
            seen.add(key)
            unique.append(c)

    return unique


def predict_mhc_binding(candidates: list[PeptideCandidate], allele: str = "HLA-A*02:01") -> list[PeptideCandidate]:
    """Score MHC-I binding using MHCflurry."""
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    for c in candidates:
        try:
            df = predictor.predict(peptides=[c.peptide], alleles=[allele], verbose=0)
            c.mhc_affinity_nm = float(df["affinity"].values[0])
            c.mhc_presentation = float(df["presentation_score"].values[0])
        except Exception:
            c.mhc_affinity_nm = 50000.0
            c.mhc_presentation = 0.0

    return candidates


def score_immunogenicity(candidates: list[PeptideCandidate]) -> list[PeptideCandidate]:
    """Score immunogenicity using the neoantigen scorer (daemon-optimized)."""
    # Import the scorer that the daemon is optimizing
    scorer_path = Path(__file__).parent.parent / "experiments" / "neoantigen-scorer"
    if str(scorer_path) not in sys.path:
        sys.path.insert(0, str(scorer_path))
    from target import score_peptide

    for c in candidates:
        c.immunogenicity_score = score_peptide(c.peptide, c.wildtype)

    return candidates


def optimize_codons(peptide: str) -> str:
    """Design codon-optimized mRNA sequence for a peptide.

    Uses human codon usage frequencies to select the most efficiently
    translated codon for each amino acid. Avoids rare codons that
    slow translation and destabilize mRNA.
    """
    mrna = []
    for aa in peptide:
        codons = CODON_TABLE.get(aa, [])
        if not codons:
            continue
        # Pick the codon with highest human usage frequency
        best_codon = max(codons, key=lambda c: HUMAN_CODON_FREQ.get(c, 0))
        mrna.append(best_codon)
    return ''.join(mrna)


def design_mrna_construct(peptide: str, include_utr: bool = True) -> str:
    """Design a full mRNA construct encoding the neoantigen peptide.

    Structure: 5'UTR — Kozak — START — peptide codons — STOP — 3'UTR

    For a real vaccine, multiple epitopes would be concatenated with
    linker sequences (GPGPG or AAY), flanked by signal peptides for
    MHC-I presentation.
    """
    # Codon-optimized coding sequence
    cds = optimize_codons(peptide)

    if include_utr:
        # Kozak consensus for efficient translation initiation
        kozak = "GCCACC"
        start = "ATG"  # methionine start codon
        stop = "TGA"   # most common stop in human genes

        # Minimal 5'UTR (based on BNT162b2 design principles)
        utr5 = "AGAATAAACTAGTATTCTTCTGGTCCCCACAGACTCAGAGAG"

        # 3'UTR with stabilizing elements (simplified from BNT162b2)
        utr3 = "GCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCC"

        # Poly-A tail (critical for mRNA stability)
        poly_a = "A" * 100

        return f"{utr5}{kozak}{start}{cds}{stop}{utr3}{poly_a}"
    else:
        return f"ATG{cds}TGA"


def rank_candidates(candidates: list[PeptideCandidate], top_n: int = 20) -> list[PeptideCandidate]:
    """Rank candidates by combined score and return top N."""
    # Filter: must bind MHC (affinity < 2000 nM)
    binders = [c for c in candidates if c.mhc_affinity_nm < 2000]

    # Sort by immunogenicity score (primary), then presentation score (tiebreak)
    binders.sort(key=lambda c: (c.immunogenicity_score, c.mhc_presentation), reverse=True)

    # Add mRNA sequences for top candidates
    for i, c in enumerate(binders[:top_n]):
        c.rank = i + 1
        c.mrna_sequence = design_mrna_construct(c.peptide)

    return binders[:top_n]


def run_pipeline(mutations: list[Mutation], allele: str = "HLA-A*02:01", top_n: int = 20) -> list[PeptideCandidate]:
    """Run the full cancer vaccine pipeline.

    mutation data → peptide candidates → MHC binding → immunogenicity → mRNA design
    """
    print(f"═══ CANCER VACCINE PIPELINE ═══")
    print(f"Mutations: {len(mutations)}")
    print(f"HLA allele: {allele}")
    print()

    # Step 1: Generate peptide candidates
    all_candidates = []
    for mut in mutations:
        candidates = generate_peptide_candidates(mut)
        all_candidates.extend(candidates)
        print(f"  {mut.gene} {mut.wildtype_aa}{mut.protein_position}{mut.mutant_aa}: {len(candidates)} peptides")

    print(f"\nTotal candidates: {len(all_candidates)}")

    if not all_candidates:
        print("No candidates generated. Check protein context sequences.")
        return []

    # Step 2: MHC binding prediction
    print("\nPredicting MHC-I binding (MHCflurry)...")
    all_candidates = predict_mhc_binding(all_candidates, allele)
    binders = [c for c in all_candidates if c.mhc_affinity_nm < 2000]
    print(f"  Binders (affinity < 2000 nM): {len(binders)} / {len(all_candidates)}")

    # Step 3: Immunogenicity scoring
    print("\nScoring immunogenicity...")
    all_candidates = score_immunogenicity(all_candidates)

    # Step 4: Rank and design mRNA
    print("\nRanking and designing mRNA constructs...")
    top = rank_candidates(all_candidates, top_n)

    # Print results
    print(f"\n{'═' * 80}")
    print(f"  TOP {len(top)} VACCINE CANDIDATES")
    print(f"{'═' * 80}")
    for c in top:
        print(f"\n  #{c.rank} — {c.gene} pos {c.mutation_position}")
        print(f"  Peptide:  {c.peptide} (wt: {c.wildtype})")
        print(f"  Length:   {c.length}-mer")
        print(f"  Binding:  {c.mhc_affinity_nm:.0f} nM (presentation: {c.mhc_presentation:.4f})")
        print(f"  Immuno:   {c.immunogenicity_score:.4f}")
        print(f"  mRNA CDS: {optimize_codons(c.peptide)}")

    # Multi-epitope construct
    if len(top) >= 3:
        print(f"\n{'═' * 80}")
        print(f"  MULTI-EPITOPE mRNA VACCINE CONSTRUCT")
        print(f"{'═' * 80}")
        # Concatenate top 5 epitopes with AAY linkers (cleavage sites)
        epitopes = [c.peptide for c in top[:5]]
        multi = "AAY".join(epitopes)
        multi_mrna = design_mrna_construct(multi)
        print(f"  Epitopes: {' — AAY — '.join(epitopes)}")
        print(f"  Total length: {len(multi)} amino acids")
        print(f"  mRNA length: {len(multi_mrna)} nucleotides")
        print(f"  mRNA (first 200 nt): {multi_mrna[:200]}...")

    return top


# ── Protein context resolution ──

def resolve_protein_context(mutations: list[Mutation]) -> list[Mutation]:
    """Fill in protein_context for mutations using the protein database."""
    from protein_db import get_protein_sequence

    resolved = []
    for mut in mutations:
        if mut.protein_context:
            resolved.append(mut)
            continue

        seq = get_protein_sequence(mut.gene)
        if not seq:
            print(f"  WARNING: No protein sequence for {mut.gene}")
            resolved.append(mut)
            continue

        idx = mut.protein_position - 1  # 1-indexed → 0-indexed
        if idx < 0 or idx >= len(seq):
            print(f"  WARNING: {mut.gene} pos {mut.protein_position} out of range (len={len(seq)})")
            resolved.append(mut)
            continue

        actual_aa = seq[idx]
        if actual_aa != mut.wildtype_aa:
            print(f"  WARNING: {mut.gene} pos {mut.protein_position}: "
                  f"expected {mut.wildtype_aa}, found {actual_aa} — trying anyway")

        # Extract context window (12 residues on each side)
        window = 12
        start = max(0, idx - window)
        end = min(len(seq), idx + window + 1)

        # Build mutated context: replace the wildtype AA with the mutant
        mutated_seq = seq[:idx] + mut.mutant_aa + seq[idx + 1:]
        mut.protein_context = mutated_seq[start:end]
        resolved.append(mut)

    return resolved


# ── VCF integration ──

def mutations_from_vcf(vcf_path: str) -> list[Mutation]:
    """Parse a VCF file and return Mutation objects with protein context."""
    from vcf_parser import parse_vcf, filter_missense, parse_protein_change

    variants = parse_vcf(vcf_path)
    missense = filter_missense(variants)
    print(f"VCF: {len(variants)} variants → {len(missense)} missense mutations")

    mutations = []
    for v in missense:
        if v.protein_change:
            parsed = parse_protein_change(v.protein_change)
            if parsed:
                wt, pos, mut = parsed
                mutations.append(Mutation(
                    gene=v.gene,
                    protein_position=pos,
                    wildtype_aa=wt,
                    mutant_aa=mut,
                ))

    return resolve_protein_context(mutations)


# ── Example mutations (common cancer hotspots) ──

def example_mutations() -> list[Mutation]:
    """Generate example mutations for the most common cancer driver genes."""
    mutations = [
        Mutation("KRAS", 12, "G", "V"),    # NSCLC, pancreatic
        Mutation("KRAS", 12, "G", "D"),    # pancreatic, colorectal
        Mutation("KRAS", 12, "G", "C"),    # NSCLC (sotorasib target)
        Mutation("TP53", 175, "R", "H"),   # most common TP53 hotspot
        Mutation("TP53", 248, "R", "W"),   # second most common TP53
        Mutation("BRAF", 600, "V", "E"),   # melanoma (vemurafenib target)
        Mutation("PIK3CA", 545, "E", "K"), # breast, colorectal
        Mutation("EGFR", 858, "L", "R"),   # NSCLC (erlotinib target)
    ]
    return resolve_protein_context(mutations)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Cancer Vaccine Pipeline")
    parser.add_argument("--vcf", help="Path to tumor VCF file")
    parser.add_argument("--mutations", help="Path to mutations JSON file")
    parser.add_argument("--allele", default="HLA-A*02:01", help="HLA allele")
    parser.add_argument("--top", type=int, default=15, help="Number of top candidates")
    args = parser.parse_args()

    if args.vcf:
        mutations = mutations_from_vcf(args.vcf)
    elif args.mutations:
        with open(args.mutations) as f:
            data = json.load(f)
        mutations = resolve_protein_context([Mutation(**m) for m in data])
    else:
        mutations = example_mutations()

    results = run_pipeline(mutations, top_n=args.top)
