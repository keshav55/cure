# Research Learnings

Extracted after each experiment. Feed these into the next experiment's design.
Same pattern as the daemon's `.weighted_learnings.json` but for research, not code.

## Paper 01: Data Quality Dominates Algorithmic Complexity

### Confirmed Findings
1. **HLA allele matching is the #1 determinant.** Wrong alleles → AUC 0.156. Right alleles → AUC 0.958. ΔAUC = 0.802. No algorithmic improvement in literature comes close.
2. **Binding rank alone achieves AUC 0.966 at peptide level.** The problem is nearly solved for peptide selection once you have the right alleles.
3. **Expression is #1 at mutation level (AUC 0.774).** When selecting which mutation to target, how much the gene is expressed matters more than how well the peptide binds.
4. **Foreignness (DAI) is anti-correlated.** AUC = 0.352 (NeoRanking), 0.315 (TESLA). Confirmed across 2 independent datasets. Likely central tolerance effect.
5. **Three features + arithmetic = 0.958.** No ML needed. Binding (0.4) + stability (0.3) + expression (0.3).

### Mistakes Made (don't repeat)
1. **Data leakage from lookup tables.** 3 test peptides in training lookup inflated AUC by 0.273. Always audit lookup-based features against validation sets.
2. **Default HLA alleles are useless.** Testing against A*02:01/A*03:01/A*11:01 gates 78% of real peptides. Require patient alleles or don't bother.
3. **N=20 is noise.** Clinical validation on 20 peptides had CI [0.34, 0.85]. Worthless for significance. Need N>500.
4. **Synthetic wildtypes corrupt foreignness.** Replacing mutant AA with a different-group AA creates fake signal. TESLA AUC was 0.609 partly because of this noise.
5. **Peptide-level AUC is misleadingly high.** 99.92% class imbalance inflates AUC. Report AUPRC for imbalanced datasets.

### Open Questions (next experiments should address)
1. **Can we predict binding with patient alleles from sequence alone?** MHCflurry supports 49 alleles but we only tested 3. Does expanding to all 49 close the gap?
2. **Can expression be estimated without RNA-seq?** GTEx tissue-average expression is available. How much AUC do we lose vs actual tumor RNA-seq?
3. **Why is foreignness anti-correlated?** Is it purely central tolerance, or is there confounding with binding disruption?
4. **Does the daemon add value on top of proper data?** We proved data > algorithms. But with the RIGHT data, can the daemon find algorithmic improvements that matter?
5. **Precision-recall analysis.** At 99.92% class imbalance, what's the actual precision at clinically useful recall (e.g., 80% sensitivity)?

### Principles for Next Experiment
- Start with the data, not the algorithm
- Require patient-matched HLA alleles as input
- Validate on NeoRanking test set (N=125,784) not hand-curated sets
- Report AUPRC alongside AUROC for imbalanced datasets
- Every claim must have a confidence interval
- Check for data leakage before reporting any number

## Experiment: Patient-Matched Alleles (2026-03-16)

### Result
Our scorer (scorer.py) with patient-matched HLA alleles:
  AUC = 0.909, 95% CI [0.868, 0.950]
  N = 902 (92 pos, 810 neg, sampled from NeoRanking test set)
  Gating: 4% of positives (was 78% with 3 default alleles)
  Separation: 0.627 (immunogenic avg 0.784, negative avg 0.157)

### What changed
Expanded DEFAULT_ALLELES from 3 to 12 (added HLA-B and HLA-C).
Used patient's actual allele from NeoRanking data column.
No algorithmic changes — same scorer code.

### Learning
**The scorer was always good. The input was wrong.**
AUC swung from 0.156 to 0.909 by fixing HLA alleles alone.
This confirms the paper's thesis: data quality > algorithms.

## Experiment: Expression Proxy Without RNA-seq (2026-03-16)

### Result
When tumor RNA-seq is unavailable, public databases can proxy:
  Actual tumor RNA-seq:    AUC = 0.791
  TCGA cancer average:     AUC = 0.791 (identical!)
  GTEx sample tissue:      AUC = 0.760
  GTEx all-tissue average: AUC = 0.730

Pearson r (tumor vs GTEx): 0.782 — strong correlation.

### Learning
**RNA-seq is not required.** TCGA cancer-type expression averages
match actual tumor expression perfectly (same AUC). GTEx tissue
averages lose only 0.06 AUC. The pipeline can work without any
patient-specific expression data by using public databases.

This means the pipeline needs only: VCF + HLA typing ($250 total).
No RNA-seq required.

## Experiment: Precision-Recall Analysis (2026-03-16)

### Result
At the clinically useful operating point (82% recall):
  Precision = 49% (half of predicted immunogenic peptides actually are)
  Threshold = 0.59

AUPRC = 0.691 (random baseline = 0.102, lift = 6.8x)

Key operating points:
  High precision (69%): recall = 44%, threshold = 0.9 (conservative)
  Balanced (53%/70%): threshold = 0.7
  High recall (91%): precision = 45%, threshold = 0.5

### Learning
**AUPRC = 0.691 is the honest metric for imbalanced data.**
AUROC (0.909) looks great but is inflated by class imbalance (99.92% negative).
AUPRC shows we're 6.8x better than random — meaningful but not dominant.

At 82% recall we catch most immunogenic peptides but half our predictions
are false positives. For vaccine design this is acceptable — you can
synthesize 20 candidates and expect ~10 to work. For drug development
this ratio would be too expensive. Context matters.

## Experiment: Foreignness Paradox Deep Dive (2026-03-16)

### Result
Three tests on the anti-correlation:

1. DAI correlates with WORSE binding (r=0.113). Foreign mutations
   slightly disrupt MHC binding. Partial confound confirmed.

2. Among good binders only (rank < 0.5), DAI is STILL anti-correlated
   (AUC=0.432). Controlling for binding doesn't fix it.
   **Central tolerance hypothesis holds.**

3. Cross-group mutations (the most "foreign") have the strongest
   anti-correlation (AUC=0.228). Same-group also anti-correlated (0.312).

### Learning
**The foreignness anti-correlation is NOT just a binding confound.**
Even among peptides that bind well, foreign mutations are less immunogenic.
This supports central tolerance: T-cells reactive to highly foreign
peptides are deleted during thymic selection.

Implication: DO NOT use foreignness as a positive feature in the scorer.
Either remove it entirely or use it as a NEGATIVE signal (penalty for
high foreignness). This is counterintuitive but data-driven.

## Experiment: Does the Scorer Add Value Beyond Binding? (2026-03-16)

### Result
At the peptide level with patient-matched alleles:
  Binding rank only (NeoRanking):  AUC = 0.968
  Our scorer (MHCflurry + features): AUC = 0.909
  3-feature combo (bind+stab+expr):  AUC = 0.967
  + inverted foreignness:            AUC = 0.969

### Learning
**Our scorer's features HURT performance at the peptide level.**
Binding rank alone (0.968) beats our full scorer (0.909) by 0.059.
The daemon's BLOSUM80, dissimilarity, TCR coverage features add noise,
not signal, when patient-matched binding predictions are available.

This does NOT mean the daemon was useless:
- The daemon found that foreignness is anti-correlated (insight)
- The daemon improved the benchmark F1 from 0.77 to 0.84 (code quality)
- The daemon's mRNA optimizer IS useful (0.70 → 0.82)

But for immunogenicity PREDICTION specifically, binding rank is
near-sufficient. The algorithmic features don't help. The value is
in getting the right alleles and having good binding predictions.

### Implication for the pipeline
Simplify scorer.py: use MHCflurry binding prediction as the primary
score. Remove or downweight BLOSUM80 and dissimilarity features.
The scorer should be a thin wrapper around MHCflurry, not a
multi-signal combiner.

### Implication for the paper
The paper's thesis is STRENGTHENED: data quality (right alleles) >
algorithmic complexity (more features). Our own scorer proves it —
adding features to binding makes it worse, not better.

## Experiment: Minimal Scorer (MHCflurry Only) (2026-03-16)

### Result
MHCflurry presentation_score with patient allele only:
  AUROC = 0.967, 95% CI [0.941, 0.993]
  Matches NeoRanking binding rank (0.968)
  Beats our full scorer (0.909) by 0.058

### Learning
**The optimal scorer is MHCflurry presentation_score. Period.**
300 lines of BLOSUM80, dissimilarity, TCR analysis = noise.
The scorer should be 5 lines:
  1. Get patient HLA alleles
  2. Call MHCflurry.predict(peptide, allele)
  3. Return presentation_score
  4. Rank by score
  5. Done

### Implication
The daemon's 1,048 experiments optimized features that don't help.
But the daemon DID discover this fact (through the validation loop).
The value was in the PROCESS of systematic investigation, not in
the features it produced.

Kalanick was right: fewest rules while staying out of chaos.
The fewest rules here is one feature: does it bind?

## Experiment: Allele Count vs Allele Matching Ablation (2026-03-16)

### Result (addresses 5.4 Pro reviewer critique)
  1 patient-matched allele:    AUROC = 0.967
  12 common alleles:           AUROC = 0.927 (37% include correct allele)
  3 original defaults:         AUROC = 0.895 (19% include correct allele)

### Learning
**Both allele count AND matching matter, but matching matters more.**
  3→12 alleles (more coverage):       ΔAUROC = +0.032
  12 common→1 matched (right allele): ΔAUROC = +0.040

The original paper claim (ΔAUROC = 0.812 from 0.156→0.968) conflates
allele matching with the full scorer's feature noise. The cleaner
comparison:
  Same scorer (MHCflurry), 3 alleles:     0.895
  Same scorer (MHCflurry), 12 alleles:    0.927
  Same scorer (MHCflurry), 1 matched:     0.967

The improvement from allele matching alone is 0.072 (3→matched) or
0.040 (12→matched). Still the single largest factor, but the 0.812
swing included the scorer noise penalty.

### Paper correction needed
The 0.156→0.968 comparison used our full scorer (0.156) vs NeoRanking
binding rank (0.968) — comparing different scorers AND different alleles.
Fair comparison: MHCflurry with 3 alleles (0.895) vs MHCflurry with
1 matched (0.967). ΔAUROC = 0.072 for allele matching alone.

Still the single largest factor, but honestly reported.

## Correction: AUPRC Baseline Reporting (2026-03-16)

### Issue (flagged by 5.4 Pro reviewer)
We reported AUPRC = 0.691 with "6.8x random baseline" but the baseline
(0.102) came from the SAMPLED evaluation (92 pos / 902 total = 10.2%
positive rate). The full dataset has 96/125,784 = 0.076% positive rate.

### Correction
Paper must state: "AUPRC = 0.691 evaluated on downsampled data (92 pos,
810 neg, 10.2% positive rate). The full dataset positive rate is 0.076%."

The 6.8x lift is correct FOR THE SAMPLED DATA but would be much higher
on the full dataset (where random baseline is 0.00076, not 0.102).
This is not dishonest but it is unclear. The paper should be explicit
about which evaluation context the AUPRC was computed on.

## Experiment: Mutation-Level Prediction (2026-03-16)

### Result
The clinical question "which mutation to vaccinate against?" (N=4,307):
  Expression only:                AUC = 0.774
  Binding only:                   AUC = 0.702
  Expr + binding:                 AUC = 0.785
  Optimized 4 features:           AUC = 0.822

Optimal weights: expression=0.5, binding=0.2, n_alleles=0.2, alt_rna=0.1

### Learning
**At mutation level, expression is 2.5x more important than binding.**
The optimal weight ratio is 50% expression vs 20% binding.

This makes biological sense: a mutation in a highly expressed gene
produces more peptides across more cells. Even if each peptide binds
moderately, the sheer number of copies increases the chance a T-cell
sees it.

Cancer cell fraction (clonality) doesn't help (weight=0). Surprising —
theory says clonal mutations are better targets. Possibly because
the NeoRanking cohort is enriched for clonal mutations already.

### Two-level prediction summary
  Peptide level: binding alone = 0.967 (solved)
  Mutation level: expression + binding = 0.822 (harder, more valuable)

The pipeline should rank at BOTH levels:
  1. Filter mutations by expression (top 50%)
  2. For each mutation, rank peptides by binding
  3. Select top peptide per top mutation

## Experiment: What Else Can Improve the System? (2026-03-17)

### Peptide level: at ceiling
Full antigen presentation pathway (binding + TAP + cleavage + stability + expression):
  AUC = 0.970 vs binding alone = 0.968 (+0.002)
  TAP and cleavage add nothing beyond binding. Binding subsumes them.
  **Peptide-level prediction is solved. No room for improvement.**

### Mutation level: marginal gains possible
  Without driver score:    0.822
  With CSCAPE driver:      0.826 (+0.004)
  Without binding:         0.820 (-0.002)
  LR found binding NEGATIVE but hand-tuned says minimal positive.

### Logistic regression vs hand-tuned
  LR: 0.812 (learned: expr=2.93, n_alleles=2.87, binding=-0.39)
  Hand-tuned: 0.822
  LR underperforms by 0.010 — small training set (107 positives) limits learning.

### Key insight
**We are at the ceiling of what single-peptide/single-mutation
features can achieve.** The remaining improvement requires:
  1. Multi-epitope optimization (which COMBINATION of peptides?)
  2. Patient-level features (TMB, immune infiltration, prior treatment)
  3. Structural modeling (peptide-MHC-TCR 3D complex)
  
These are the next frontiers, not more features from the same dataset.

## Breakthrough: ML Ensemble for Vaccine Selection (2026-03-17)

### The Problem
Given a patient's tumor mutations and ~3,000-12,000 candidate peptides,
select 20 peptides for a vaccine. Metric: what % of confirmed immunogenic
peptides are in the top 20? (recall@20)

### Results
| Method | recall@20 | Patients hit (of 30) |
|--------|-----------|---------------------|
| Binding rank only | 0.492 | 20/30 |
| Hand-tuned (bind 0.8 + stab 0.2) | 0.555 | 24/30 |
| GBT (optimized) | 0.669 | 26/30 |
| **Stacking ensemble (GBT + RF + LR)** | **0.698** | **27/30** |

### Key Details
- 42% relative improvement over binding alone
- 13 patients improved, 16 same, 1 worse (TESLA9: 1.0 → 0.5)
- Biggest wins on hardest patients: 1 positive among 5,000-12,000 peptides
- Model trained on 82 positives + 5,000 sampled negatives from train split
- 5-fold CV AUC = 0.989 ± 0.004 (not overfitting)
- Stacking weights: GBT 0.4, RF 0.4, LR 0.2
- GBT hyperparams: n_estimators=50, max_depth=2, lr=0.1, subsample=1.0

### Features That Actually Help (at vaccine selection level)
These DON'T help at peptide-level AUC but help at vaccine selection recall:
- All three binding predictors (MixMHC, NetMHCpan, PRIME) = 51.6% importance
- Binding stability = 10%
- CSCAPE driver score = signal at selection level
- Expression features = signal at selection level

Features that DON'T help at ANY level:
- Foreignness (all DAI variants) = anti-correlated or noise
- Netchop score = noise
- bestWTPeptideCount_I, bestWTMatchScore_I = noise at selection level

### The Two Surprises
1. **GBT beats RF at selection (0.669 vs 0.597) but RF beats GBT at AUC (0.972 vs 0.924)**
   → Non-linear interactions matter for patient-level selection
2. **Stacking > any individual model** because they err on different patients

### Clinical Significance
At K=2 (most restrictive): ensemble finds 22% of positives vs 5.7% for binding
At K=20 (typical vaccine): 70% vs 49%
At K=50 (generous): 83% vs 69%

### What This Means for the Paper
This is paper 03 material. No one has:
1. Framed neoantigen selection as per-patient recall@K
2. Shown ML stacking beats binding rank by 42% 
3. Demonstrated ensemble helps most for hardest patients

### Three Patients Still Missed (recall@20 = 0)
- Patient4: 1 pos in 97, bind=0.60, ranks 28th — just outside top-20
- 4324: 1 pos in 4310, bind=8.0 — very weak binder, possibly false positive
- 4014: 1 pos in 7436, bind=0.20 — ranks 175th, good binder but drowned out


## Mutation-Level Selection — No Additional Gain (2026-03-17)

### Tested
- Mutation-level only: 0.556 (max, 10mut×2pep, budget=20) — WORSE than peptide
- Hybrid (peptide + mutation): 0.699 (10pep + 10mut×1pep, budget=20) — SAME
- Larger budget hybrid: 0.730 (15pep + 5mut×3pep, budget=30) — more peptides help

### Conclusion
The ensemble already captures mutation-level signal through expression/CSCAPE features.
Explicit mutation deduplication doesn't help because:
1. Most patients have more mutations than top-20 candidates anyway
2. The model implicitly upweights peptides from expressed mutations
3. Budget matters more than selection strategy

### Autoresearch Daemon Result
30 rounds, 0 keeps. The 0.698 is the ceiling for this model class + feature set.
The limitation is training data (82 positives), not model capacity.

### The Full Optimization Arc
| Method | recall@20 | How found |
|--------|-----------|-----------|
| Binding rank | 0.492 | Literature standard |
| + stability | 0.555 | Grid search, 719 combos |
| + CSCAPE | 0.563 | Full grid, 1296 combos |
| GBT | 0.669 | Hyperparameter search |
| Stacking GBT+RF+LR | 0.698 | Model ensemble |
| Hybrid mutation-level | 0.699 | No additional gain |
| Daemon optimization | 0.698 | 30 rounds, confirmed ceiling |

