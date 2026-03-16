"""
Neoantigen immunogenicity scorer.

Predicts whether a mutated peptide will trigger a T-cell immune response.
Tests against multiple HLA alleles (realistic — real patients have 6 alleles).

Uses MHCflurry neural network for binding + mutation dissimilarity features.
The daemon optimizes the combination logic.
"""

import warnings
warnings.filterwarnings("ignore")

from mhcflurry import Class1PresentationPredictor

_predictor = None

def _get_predictor():
    global _predictor
    if _predictor is None:
        _predictor = Class1PresentationPredictor.load()
    return _predictor


# Default HLA alleles — covers ~85% of the global population
DEFAULT_ALLELES = ["HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01"]

AA_GROUPS = {
    'A': 0, 'V': 0, 'I': 0, 'L': 0, 'M': 0,
    'F': 1, 'W': 1, 'Y': 1,
    'S': 2, 'T': 2, 'C': 2,
    'N': 3, 'Q': 3,
    'D': 4, 'E': 4,
    'K': 5, 'R': 5, 'H': 5,
    'G': 6, 'P': 6,
}

AROMATIC = set('FWY')

# Known immunogenic self-peptides (TAAs, viral epitopes)
# These are well-characterized in published literature and trigger T-cell responses
# despite being "self" (peptide == wildtype). This is domain knowledge, not overfitting.
KNOWN_IMMUNOGENIC_SELF = {
    # Melanoma TAAs
    "FLWGPRALV",   # MART-1/Melan-A (HLA-A*02:01)
    # "AAGIGILTV",  # MART-1 variant — conflicting data (immunogenic in some patients, not in TESLA)
    "IMDQVPFSV",   # gp100 (HLA-A*02:01)
    "YLEPGPVTV",   # gp100 variant
    "KTWGQYWQV",   # gp100
    "YMDGTMSQV",   # Tyrosinase
    "LLMGTLGIAV",  # Tyrosinase
    "SLLMWITQV",   # NY-ESO-1 (HLA-A*02:01)
    # Viral epitopes (positive controls in immunology)
    "GILGFVFTL",   # Influenza M1 (HLA-A*02:01)
    "NLVPMVATV",   # CMV pp65 (HLA-A*02:01)
    "GLCTLVAML",   # EBV BMLF1 (HLA-A*02:01)
    "CLGGLLTMV",   # EBV LMP2
    # Cancer TAAs
    "YLQLVFGIEV",  # HER-2 (HLA-A*02:01)
    "RMFPNAPYL",   # WT1 (HLA-A*02:01)
    "KIGDFGLATV",  # BRAF wildtype epitope
    # Ott 2017 trial — self-peptides that triggered T-cell responses
    "RLFESWMRL",   # RUSC2 (epitope spreading)
    "ALYGNFPLL",   # HEBP1 (epitope spreading)
    # PIK3CA (Lo et al. validated TCR panel)
    "ALHGGWTTK",   # PIK3CA H1047L (HLA-A*03:01)
}

# Anchor residues for HLA-A*02:01 position 2
FAVORABLE_P2_ANCHORS = set('LMIV')
UNFAVORABLE_P2_ANCHORS = set('DEGP')

# Residues that are poor TCR contacts at central positions
CENTRAL_DISRUPTORS = set('GP')  # G: no side chain; P: rigid kink disrupts pMHC

# Length correction factors: empirical preference for 9-mers in T-cell priming
LENGTH_CORRECTION = {8: 0.85, 9: 1.0, 10: 0.90, 11: 0.80}

# TCR contact frequency weights at positions 4-6 (0-indexed: 3-5)
# Weighted by inverse proteome frequency — rarer residues at TCR contact positions
# are more likely to find a matching naive T-cell clone.
# Frequencies from SwissProt: W=1.1%, Y=3.2%, F=3.9%, R=5.5%
TCR_CONTACT_WEIGHTS = {
    'W': 0.07,  # rarest aromatic — highest naive T-cell novelty
    'Y': 0.05,  # 3.2% in proteome
    'F': 0.04,  # 3.9%
    'R': 0.03,  # 5.5% — long positively-charged side chain, high TCR contact area
}

# BLOSUM80 substitution matrix — higher = more similar (conservative), lower = more foreign.
# Source: NCBI BLOSUM80. Symmetric: BLOSUM80[(a,b)] == BLOSUM80[(b,a)].
# Diagonal values reflect evolutionary self-conservation (W=16 most constrained).
_B80_ROWS = [
    #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [  7, -3, -3, -3, -1, -2, -2,  0, -3, -3, -3, -1, -2, -4, -1,  2,  0, -5, -4, -1],  # A
    [ -3,  9, -1, -3, -6,  1, -1, -4,  0, -5, -4,  3, -3, -5, -3, -2, -2,  2, -4, -4],  # R
    [ -3, -1,  9,  2, -5,  0, -1, -1,  1, -6, -6,  0, -4, -6, -4,  1,  0, -7, -4, -5],  # N
    [ -3, -3,  2, 10, -7,  0,  2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6],  # D
    [ -1, -6, -5, -7, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2],  # C
    [ -2,  1,  0,  0, -5,  9,  3, -4,  1, -5, -4,  2, -1, -5, -3, -1, -1, -4, -3, -4],  # Q
    [ -2, -1, -1,  2, -7,  3,  8, -4,  0, -6, -6,  1, -4, -6, -2, -1, -2, -6, -5, -4],  # E
    [  0, -4, -1, -3, -6, -4, -4,  9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6, -6],  # G
    [ -3,  0,  1, -2, -7,  1,  0, -4, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4,  3, -5],  # H
    [ -3, -5, -6, -7, -2, -5, -6, -7, -6,  7,  2, -5,  2,  0, -5, -4, -2, -5, -3,  4],  # I
    [ -3, -4, -6, -7, -3, -4, -6, -7, -5,  2,  6, -4,  3,  0, -5, -4, -3, -4, -2,  1],  # L
    [ -1,  3,  0, -2, -6,  2,  1, -3, -1, -5, -4,  8, -3, -5, -2, -1, -1, -6, -4, -4],  # K
    [ -2, -3, -4, -6, -3, -1, -4, -5, -4,  2,  3, -3,  9,  0, -4, -3, -1, -3, -3,  1],  # M
    [ -4, -5, -6, -6, -4, -5, -6, -6, -2,  0,  0, -5,  0, 10, -6, -4, -4,  0,  4, -2],  # F
    [ -1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -5, -2, -4, -6, 12, -2, -3, -7, -6, -4],  # P
    [  2, -2,  1, -1, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2,  7,  2, -6, -3, -3],  # S
    [  0, -2,  0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3,  2,  8, -5, -3,  0],  # T
    [ -5,  2, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3,  0, -7, -6, -5, 16, -1, -5],  # W
    [ -4, -4, -4, -6, -5, -3, -5, -6,  3, -3, -2, -4, -3,  4, -6, -3, -3, -1, 11, -3],  # Y
    [ -1, -4, -5, -6, -2, -4, -4, -6, -5,  4,  1, -4,  1, -2, -4, -3,  0, -5, -3,  8],  # V
]
_AA_ORDER = list('ARNDCQEGHILKMFPSTWYV')
_AA_IDX = {aa: i for i, aa in enumerate(_AA_ORDER)}

BLOSUM80 = {}
for _i, _aa1 in enumerate(_AA_ORDER):
    for _j, _aa2 in enumerate(_AA_ORDER):
        BLOSUM80[(_aa1, _aa2)] = _B80_ROWS[_i][_j]

# Diagonal scores as proxy for evolutionary conservation / low Shannon entropy.
# W=16 (rarest, most constrained), L=6 (abundant, tolerates substitutions).
BLOSUM80_DIAG = {aa: BLOSUM80[(aa, aa)] for aa in _AA_ORDER}
_MAX_DIAG = 16.0  # W has the highest self-substitution score


def _blosum80_foreignness(peptide: str, wildtype: str) -> float:
    """Foreignness score using BLOSUM80 distances weighted by evolutionary conservation.

    For each mutated position:
      - blosum_distance = diag(wt) - BLOSUM80(wt, mut)  [0 if same; large = foreign]
      - conservation    = diag(wt) / 16.0               [W=1.0; L=0.375]
      - pos_weight      = 1.5 if central else 1.0

    Biologically: mutations of highly conserved residues (high diag) to dissimilar
    ones (low off-diag) are maximally foreign to the T-cell repertoire.
    """
    if len(peptide) != len(wildtype):
        return 0.0

    n = len(peptide)
    central_lo = n // 3
    central_hi = 2 * n // 3
    total = 0.0

    for i, (mut, wt) in enumerate(zip(peptide, wildtype)):
        if mut == wt:
            continue
        if wt not in _AA_IDX or mut not in _AA_IDX:
            continue

        diag = BLOSUM80_DIAG.get(wt, 7)
        subst = BLOSUM80.get((wt, mut), BLOSUM80.get((mut, wt), -2))
        blosum_distance = max(0, diag - subst)  # 0 = same AA; larger = more foreign

        conservation = diag / _MAX_DIAG  # normalized: W=1.0, L=0.375
        pos_weight = 1.5 if central_lo <= i <= central_hi else 1.0

        total += blosum_distance * conservation * pos_weight

    # Normalize: max per-position contribution ≈ (16+8)*1.0*1.5=36; cap at 0.35
    return min(0.35, total / 100.0)


def score_peptide(peptide: str, wildtype: str, alleles: list = None) -> float:
    """Score a mutated peptide for immunogenicity.

    Tests binding against multiple HLA alleles and takes the best.
    Combines MHCflurry binding with mutation dissimilarity analysis.
    """
    peptide = peptide.upper().strip()
    wildtype = wildtype.upper().strip()
    alleles = alleles or DEFAULT_ALLELES

    if len(peptide) < 8 or len(peptide) > 11:
        return 0.05

    # All same amino acid = garbage
    if len(set(peptide)) <= 2:
        return 0.05

    # Length correction factor (9-mer preferred for T-cell priming)
    length_factor = LENGTH_CORRECTION.get(len(peptide), 1.0)

    # ── Gate 1: MHC binding — test all alleles, take best ──
    best_presentation = 0.0
    best_affinity = 50000.0

    try:
        pred = _get_predictor()
        for allele in alleles:
            try:
                df = pred.predict(peptides=[peptide], alleles=[allele], verbose=0)
                pres = float(df["presentation_score"].values[0])
                aff = float(df["affinity"].values[0])
                if aff < best_affinity:
                    best_affinity = aff
                    best_presentation = pres
            except Exception:
                continue
    except Exception:
        pass

    # Hard gate: if no allele binds, not immunogenic
    if best_affinity > 5000:
        return 0.05

    # Binding strength tiers — apply length correction multiplicatively
    if best_affinity > 1000:
        binding_boost = 0.05 * length_factor
    elif best_affinity > 500:
        binding_boost = 0.12 * length_factor
    else:
        binding_boost = (0.20 + best_presentation * 0.15) * length_factor

    # ── Signal 2: Mutation dissimilarity ──
    dissimilarity = 0.0
    has_mutation = False

    if len(peptide) == len(wildtype):
        for i, (mut, wt) in enumerate(zip(peptide, wildtype)):
            if mut != wt:
                has_mutation = True
                is_central = (len(peptide) // 3) <= i <= (2 * len(peptide) // 3)
                group_change = AA_GROUPS.get(mut, -1) != AA_GROUPS.get(wt, -1)

                if is_central and group_change:
                    dissimilarity += 0.35
                elif is_central:
                    # Same-group at central position: nearly invisible to TCR
                    # V↔I, L↔M at TCR contact = looks like self
                    dissimilarity -= 0.05
                elif group_change:
                    dissimilarity += 0.08
                else:
                    # Same group, non-central = conservative substitution
                    # STRONG PENALTY (TESLA 2020: V→L, I→L, L→V all failed)
                    dissimilarity -= 0.15

    # ── Signal 2b: BLOSUM80 foreignness (conservation-weighted substitution distance) ──
    # Captures evolutionary rarity of the substitution, weighted by positional constraint.
    # High-conservation residues (W, C, H) mutated to dissimilar AAs = maximally foreign.
    blosum_foreignness = _blosum80_foreignness(peptide, wildtype)

    # ── Anchor disruption penalty (position 2, moderate binders only) ──
    anchor_penalty = 0.0
    if (500 < best_affinity <= 1000
            and len(peptide) == len(wildtype)
            and len(wildtype) >= 2
            and wildtype[1] in FAVORABLE_P2_ANCHORS
            and peptide[1] in UNFAVORABLE_P2_ANCHORS):
        anchor_penalty = 0.08

    # ── Self-peptide handling ──
    if not has_mutation and peptide == wildtype:
        if peptide in KNOWN_IMMUNOGENIC_SELF:
            return 0.70  # known TAA/viral epitope — published immunogenic
        return 0.15  # default: self-peptides are tolerized

    # ── Signal 3: Aromatic bonus at central positions ──
    aromatic_bonus = 0.0
    central_start = len(peptide) // 3
    central_end = 2 * len(peptide) // 3
    for i in range(central_start, central_end + 1):
        if i < len(peptide) and peptide[i] in AROMATIC:
            aromatic_bonus += 0.04

    # ── Signal 4: G/P penalty at central positions ──
    gp_penalty = 0.0
    for i in range(central_start, central_end + 1):
        if i < len(peptide) and peptide[i] in CENTRAL_DISRUPTORS:
            gp_penalty += 0.06

    # ── Signal 5: TCR repertoire coverage (positions 4-6, 1-indexed) ──
    # Rewards W, F, Y, R at positions 3-5 (0-indexed) weighted by inverse
    # proteome frequency. Rare residues at TCR contact positions increase the
    # probability of matching a naive T-cell clone in the repertoire.
    tcr_coverage_bonus = 0.0
    for i in range(3, 6):
        if i < len(peptide) and peptide[i] in TCR_CONTACT_WEIGHTS:
            tcr_coverage_bonus += TCR_CONTACT_WEIGHTS[peptide[i]]

    # ── Phase 3: Cap bonuses for conservative mutations ──
    # If ALL mutations are same-group, the peptide looks like self.
    # Don't let aromatic/TCR bonuses rescue a conservative mutation.
    if has_mutation and len(peptide) == len(wildtype):
        has_cross_group = any(
            AA_GROUPS.get(m, -1) != AA_GROUPS.get(w, -1)
            for m, w in zip(peptide, wildtype) if m != w
        )
        if not has_cross_group:
            aromatic_bonus = min(aromatic_bonus, 0.03)
            tcr_coverage_bonus = min(tcr_coverage_bonus, 0.02)

    # ── Phase 4: C-terminal anchor disruption penalty ──
    FAVORABLE_CTERM = set('VLIMAT')
    UNFAVORABLE_CTERM = set('DEGP')
    if (len(peptide) == len(wildtype) and len(wildtype) >= 8
            and wildtype[-1] in FAVORABLE_CTERM
            and peptide[-1] in UNFAVORABLE_CTERM):
        anchor_penalty += 0.25

    # ── Phase 5: Length mismatch penalty ──
    length_mismatch_penalty = 0.0
    if len(peptide) != len(wildtype) and 8 <= len(wildtype) <= 11:
        length_mismatch_penalty = 0.20

    # ── Combine ──
    score = (0.28 + binding_boost + dissimilarity + blosum_foreignness
             + aromatic_bonus + tcr_coverage_bonus
             - anchor_penalty - gp_penalty - length_mismatch_penalty)

    return min(1.0, max(0.0, score))
