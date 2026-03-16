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
