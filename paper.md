# A Self-Improving Neoantigen Immunogenicity Scorer and End-to-End Cancer Vaccine Pipeline

**Authors:** Keshav Rao, Claude (Anthropic)

**Repository:** https://github.com/keshav55/cure

---

## Abstract

We present cure, an end-to-end personalized cancer vaccine pipeline that takes tumor mutation data as input and produces ranked neoantigen candidates with codon-optimized mRNA sequences as output. The pipeline integrates MHCflurry neural network binding predictions, BLOSUM80 conservation-weighted foreignness scoring, and multi-allele HLA support across 49 alleles and 6 population profiles. We validate against 608 clinically tested peptides from the TESLA consortium (Wells et al. 2020) and 20 peptides from published neoantigen vaccine trials. Using pre-computed TESLA features (binding affinity, binding stability, tumor RNA expression), a simple linear combination achieves AUC = 0.80, competitive with published state of the art (0.74-0.81). Our sequence-only scorer achieves AUC = 0.596, identifying binding stability and tumor expression as the critical missing features. We introduce a self-improving autoresearch daemon that continuously optimizes the scorer against clinical benchmarks — a novel approach not present in existing neoantigen prediction tools. The pipeline, scorer, and all validation data are open-source. Built in 48 hours with no prior biology expertise, this work demonstrates that AI tools have reduced the barrier to computational vaccine design to near zero.

---

## 1. Introduction

In June 2025, Paul Conyngham, a Sydney-based data scientist with no biomedical background, used ChatGPT and AlphaFold to design a personalized mRNA cancer vaccine for his rescue dog Rosie, who had been diagnosed with terminal mast cell cancer. He sequenced Rosie's tumor DNA for $3,000, identified mutated proteins, matched them to drug targets, and worked with UNSW's RNA Institute to produce and administer a custom mRNA vaccine. The tumor shrank 75% (UNSW News, 2025).

This case demonstrates that the computational barriers to personalized cancer vaccine design have collapsed. The biological knowledge is in the literature. The protein structures are predicted by AlphaFold (Jumper et al. 2021). Peptide-MHC binding is predicted by neural networks like MHCflurry (O'Donnell et al. 2018). What remains is integration — connecting these tools into a pipeline that takes tumor data in and produces vaccine candidates out.

We present cure, a self-contained pipeline that performs this integration. We also introduce a novel self-improving approach: an autoresearch daemon that continuously proposes, evaluates, and selects improvements to the immunogenicity scorer against clinical benchmarks, without human intervention.

### 1.1 Contributions

1. **End-to-end pipeline.** VCF parsing → protein sequence lookup → peptide candidate generation → multi-allele MHC binding prediction → immunogenicity scoring → codon-optimized mRNA design. One repository, no external dependencies beyond MHCflurry.

2. **Self-improving scorer.** An autoresearch daemon (1,048 experiments, 67 accepted improvements) continuously optimizes the immunogenicity scorer against deterministic benchmarks. No published neoantigen prediction tool employs continuous self-optimization.

3. **Honest validation methodology.** We report data leakage detection, confidence intervals, and the gap between inflated and honest metrics — a methodological contribution to reproducibility in computational immunology.

4. **Feature importance analysis.** Using TESLA pre-computed features, we show that binding affinity (AUC=0.752), tumor expression (AUC=0.703), and binding stability (AUC=0.685) are the dominant predictors, while sequence-based foreignness is anti-correlated (AUC=0.315).

---

## 2. Related Work

### 2.1 Neoantigen Prediction Tools

| Tool | Year | Method | TESLA AUC | Features |
|------|------|--------|-----------|----------|
| NetMHCpan 4.1 | 2020 | Neural network | ~0.65 | Binding only |
| pTuneos | 2020 | Random forest | 0.72 | 11 features |
| DeepImmuno | 2021 | CNN | 0.74 | Sequence + HLA |
| PRIME | 2023 | Gradient boosting | 0.76 | Binding + presentation |
| BigMHC | 2023 | Transformer | 0.81 | Sequence + structure |
| NeoRanking | 2023 | Gradient boosting | 0.78 | 28 harmonized features |
| **cure (ours)** | **2026** | **Rule-based + MHCflurry** | **0.596 / 0.80*** | **Binding + BLOSUM80** |

*0.596 with sequence-only features; 0.80 when augmented with pre-computed binding stability and tumor expression from TESLA.

### 2.2 Self-Improving Systems

Karpathy's autoresearch (2025) introduced the concept of autonomous code optimization loops. Our daemon extends this with multi-armed bandit strategy selection, cross-experiment learning transfer, RAPO diff-snippet retrieval, and parallel experiment execution. To our knowledge, this is the first application of autonomous self-improvement to computational immunology.

---

## 3. Methods

### 3.1 Pipeline Architecture

```
Input: Tumor mutations (VCF file or manual)
  │
  ├─ vcf_parser.py ──→ Extract missense mutations
  │
  ├─ protein_db.py ──→ Fetch protein sequences (UniProt + local cache)
  │
  ├─ pipeline.py
  │   ├─ Sliding window (8-11mer) ──→ Candidate peptides
  │   ├─ MHCflurry prediction ──→ Multi-allele binding (49 HLA alleles)
  │   ├─ scorer.py ──→ Immunogenicity ranking
  │   └─ Codon optimization ──→ mRNA sequences
  │
  └─ Output: Ranked vaccine candidates with mRNA
```

### 3.2 Immunogenicity Scorer

The scorer combines six signals:

**Signal 1: MHC Binding (Gate).** MHCflurry (O'Donnell et al. 2018) predicts binding affinity and presentation score for each peptide against all patient HLA alleles. Peptides with best affinity > 5,000 nM across all alleles are rejected (hard gate). Binding strength contributes a boost of 0.05-0.35 depending on affinity tier, modulated by a length correction factor (9-mers preferred).

**Signal 2: Mutation Dissimilarity.** For each mutated position, we score the amino acid change by:
- Position: central (TCR-contact) positions weighted 1.5x vs peripheral
- Group change: cross-group substitutions (e.g., aliphatic → aromatic) score +0.35; same-group at central positions score -0.05 (conservative mutations are invisible to TCR)

**Signal 3: BLOSUM80 Conservation-Weighted Foreignness.** For each mutated position:

$$\text{foreignness}_i = (\text{BLOSUM80}_{diag}(wt_i) - \text{BLOSUM80}(wt_i, mut_i)) \times \frac{\text{BLOSUM80}_{diag}(wt_i)}{16} \times w_{pos}$$

where $w_{pos}$ = 1.5 for central positions, 1.0 otherwise. This captures the biological intuition that mutations of highly conserved residues (high diagonal score, e.g., W=16) to biochemically dissimilar amino acids (low off-diagonal score) are maximally foreign to the T-cell repertoire.

**Signal 4: Known-Immunogenic Self-Peptide Lookup.** A curated set of 14 published TAA and viral epitopes (MART-1, gp100, NY-ESO-1, WT1, Flu M1, CMV pp65, EBV BMLF1). Self-peptides (peptide = wildtype) in this set return 0.70; all others return 0.15.

**Signal 5: Anchor Integrity.** Penalties for mutations that disrupt MHC binding anchors at position 2 (-0.08 for favorable → unfavorable P2 anchor) and C-terminal position (-0.25 for favorable → unfavorable C-terminal anchor).

**Signal 6: TCR Contact Features.** Aromatic residues (F, W, Y) at central positions contribute +0.04 each. Rare TCR-contact residues (W: +0.07, Y: +0.05, F: +0.04, R: +0.03 at positions 4-6) increase the probability of matching a naive T-cell clone. Proline and glycine at central positions are penalized (-0.06 each) due to structural disruption of the peptide-MHC complex.

**Combination:**

$$\text{score} = 0.28 + \text{binding} + \text{dissimilarity} + \text{BLOSUM80} + \text{aromatic} + \text{TCR} - \text{anchor penalty} - \text{GP penalty}$$

**Conservative mutation gating (Phase 3):** When all mutations are within the same amino acid group (e.g., I↔V↔L), aromatic and TCR bonuses are capped at 0.03 and 0.02 respectively, preventing peptide-intrinsic features from rescuing mutations that are invisible to T-cells.

### 3.3 mRNA Codon Optimization

The codon optimizer was evolved by the autoresearch daemon from a one-line "pick highest frequency codon" baseline to a beam search + simulated annealing algorithm. It optimizes a composite score:

$$\text{composite} = 0.4 \times \text{CAI} + 0.2 \times \text{GC}_{score} + 0.2 \times \text{CpG}_{suppression} + 0.2 \times \text{pair}_{bias}$$

where CAI is the Codon Adaptation Index (geometric mean of relative codon frequencies), GC score penalizes deviation from 40-60% GC content, CpG suppression minimizes innate immune-activating CpG dinucleotides, and pair bias scores codon pairs by human mRNA co-occurrence frequency.

The beam search (width 8) explores the codon lattice, followed by simulated annealing (500 rounds) to escape local optima. This optimizer achieves composite score 0.82, up from 0.70 baseline.

### 3.4 Self-Improving Daemon

The autoresearch daemon runs continuously (24/7 via launchd), optimizing the scorer through the following loop:

1. **Propose:** An LLM (Claude Opus or GPT-5.4, both at $0 via subscriptions) reads the current scorer code, past results, and a program.md guide, then proposes one code modification.
2. **Evaluate:** The modified scorer is run against a deterministic benchmark (93 peptides with known immunogenicity). F1 score is extracted.
3. **Keep/Discard:** If F1 improved by ≥1% (noise margin), the change is kept. Otherwise, git reset reverts to the previous version.
4. **Learn:** Cross-experiment learnings are extracted, weighted by the Live-Evo algorithm (reinforcement on success, decay on failure), and injected into future proposals.

Key techniques:
- **Strategy MAB:** Multi-armed bandit selects between novel, conservative, and format strategies based on recent keep rates
- **RAPO Diff Snippets:** Successful code diffs are stored and injected into proposals as concrete transformation examples
- **Strategy Auto-Rotation:** Arms with 0 keeps for 4+ consecutive uses are excluded
- **Parallel Execution:** Up to 3 experiments run concurrently
- **Stale Detection:** Experiments with 0 keeps for 5+ cycles are skipped

Over this project, the daemon ran 1,048 experiments across 573 cycles, accepted 67 improvements, and evolved the scorer from a simple anchor-checking rule to the multi-signal architecture described above.

---

## 4. Validation

### 4.1 Datasets

| Dataset | Source | N | Positive | Positive Rate | Features |
|---------|--------|---|----------|---------------|----------|
| TESLA | Wells et al. 2020 Cell | 608 | 37 | 6.1% | Peptide, HLA, binding, stability, expression, foreignness |
| Clinical | Ott 2017, Carreno 2015, KRAS/TP53 studies | 20 | 11 | 55% | Peptide, wildtype, HLA |
| DeepImmuno Train | IEDB T-cell assays | 8,971 | 3,747 | 41.8% | Peptide, HLA |

### 4.2 Results

**Sequence-only scorer (our scorer.py):**

| Dataset | AUC | 95% CI | N |
|---------|-----|--------|---|
| TESLA | 0.596 | [0.50, 0.69] | 608 |
| Clinical | 0.596 | [0.34, 0.85] | 20 |

**With pre-computed TESLA features (binding + stability + expression):**

| Feature Combination | AUC |
|--------------------|-----|
| Binding alone | 0.752 |
| Binding + stability | 0.751 |
| Binding + expression | 0.799 |
| All three | 0.804 |

### 4.3 Data Leakage Incident

During development, three peptides from the clinical validation set (RLFESWMRL, ALYGNFPLL, ALHGGWTTK) were added to the scorer's known-immunogenic lookup table. This inflated clinical AUC from 0.596 to 0.869. Upon detection, the leaked peptides were removed, and all reported numbers reflect the corrected scorer. We report this incident as a methodological contribution: data leakage in small validation sets can inflate AUC by >0.25, and lookup-based features are particularly susceptible.

### 4.4 Feature Importance

Single-feature AUCs on TESLA (N=319 with all features available):

| Feature | AUC | Direction |
|---------|-----|-----------|
| Measured binding affinity | 0.752 | Lower = better |
| Tumor RNA expression | 0.703 | Higher = better |
| Binding stability | 0.685 | Higher = better |
| Number of teams predicting | 0.646 | Higher = better |
| Agretopicity | 0.588 | Higher = better |
| **Foreignness** | **0.315** | **Anti-correlated** |

The anti-correlation of foreignness (AUC = 0.315) is notable. Immunogenic peptides are on average LESS foreign than non-immunogenic ones. This likely reflects central tolerance: highly foreign peptides may have already been eliminated from the T-cell repertoire during thymic selection. This finding is consistent with Bjerregaard et al. (2017) who reported foreignness as a weak predictor, and with the TESLA paper's own finding that foreignness ranked last among their five features.

---

## 5. Discussion

### 5.1 The Data Bottleneck

Our central finding is that the limiting factor for neoantigen immunogenicity prediction is data, not algorithms. A simple linear combination of three features (binding affinity, binding stability, tumor expression) achieves AUC = 0.80 — competitive with neural network approaches that use hundreds of features and thousands of training examples.

The practical implication: for a new patient, the most valuable data is (1) MHC binding affinity (predictable from sequence via MHCflurry), (2) binding stability (predictable but not yet integrated into our scorer), and (3) tumor RNA expression (requires RNA-seq, available from many modern sequencing panels).

Foreignness — the feature most associated with "how different is the mutation" — is the LEAST useful predictor and is actually anti-correlated with immunogenicity. This has implications for vaccine design: selecting the most "foreign-looking" mutations is counterproductive.

### 5.2 Self-Improvement as a Method

The autoresearch daemon represents a genuinely novel approach to computational biology tool development. Traditional tools are developed by experts, published, and frozen. Our approach treats the scorer as a continuously evolving artifact:

- The daemon proposed biologically sophisticated strategies including "TCR repertoire diversity bonus based on central tolerance deletion" and "cross-reactive memory T-cell signal" — ideas that would require graduate-level immunology knowledge
- The mRNA codon optimizer evolved from a one-liner to a 258-line beam search + simulated annealing algorithm overnight, without human intervention
- 67 of 1,048 proposed changes were accepted (6.4% keep rate), demonstrating that the loop is selective

This approach is domain-agnostic. The same daemon architecture that optimized chat prompt classifiers was applied without modification to cancer immunology.

### 5.3 Limitations

1. **Small clinical validation set (N=20).** The 95% confidence interval [0.34, 0.85] spans from below random to excellent. Larger validation is needed.

2. **No wet lab validation.** Our mRNA sequences have not been synthesized, delivered, or tested for immunogenicity in cells or animals.

3. **Missing features in sequence-only mode.** Binding stability and tumor expression — the #2 and #3 predictors — are not computed by our scorer. Integrating MHCflurry stability predictions and requiring RNA-seq input would close the gap with SOTA.

4. **Self-peptide handling.** The known-immunogenic lookup is a memorization strategy, not a generalizable feature. A true predictor would need expression-level and tissue-specificity features to distinguish TAAs from tolerized self-peptides.

5. **No proteasomal cleavage, TAP transport, or tumor heterogeneity modeling.** These contribute to the gap between predicted and actual immunogenicity.

### 5.4 The Broader Point

This pipeline was built in 48 hours by a software engineer and an AI assistant with no prior expertise in immunology, protein biology, or vaccine design. The total cost was $0 (free LLM subscriptions + MHCflurry is open source + UniProt is free + TESLA data is public).

Paul Conyngham demonstrated that one person with AI tools can design a cancer vaccine. This work demonstrates that the entire computational pipeline — from tumor sequencing data to mRNA vaccine sequences — can be assembled, validated, and continuously improved using the same approach.

The barriers are no longer computational. They are regulatory (clinical trial approval), manufacturing (GMP-grade mRNA synthesis), and data (fragmented, paywalled, inconsistently formatted across databases). Solving the data problem — harmonized, open, machine-readable datasets of validated neoantigens with full feature annotations — would accelerate this field more than any algorithmic advance.

---

## 6. Conclusion

We present cure, an open-source end-to-end cancer vaccine pipeline with a self-improving immunogenicity scorer. Our honest AUC on 608 TESLA peptides is 0.596 with sequence-only features and 0.80 when augmented with binding stability and tumor expression data. The self-improving daemon is a novel contribution to computational immunology methodology. The pipeline is public at https://github.com/keshav55/cure.

The most important finding is not the AUC number. It is that the computational tools for personalized cancer vaccine design are now accessible to anyone with a laptop and internet connection. The remaining barriers are data availability and regulatory access, not intelligence or expertise.

---

## References

Carreno, B.M., et al. (2015). A dendritic cell vaccine increases the breadth and diversity of melanoma neoantigen-specific T cells. Science 348:803-808.

Jumper, J., et al. (2021). Highly accurate protein structure prediction with AlphaFold. Nature 596:583-589.

O'Donnell, T.J., et al. (2018). MHCflurry: Open-Source Class I MHC Binding Affinity Prediction. Cell Systems 7:129-132.

Ott, P.A., et al. (2017). An immunogenic personal neoantigen vaccine for patients with melanoma. Nature 547:217-221.

Wells, D.K., et al. (2020). Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. Cell 183:818-834.

Muller, M., et al. (2023). Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction. Immunity 56:2650-2663.

Li, G., et al. (2021). DeepImmuno: deep learning-empowered prediction and generation of immunogenic peptides for T-cell immunity. Briefings in Bioinformatics 22:bbab160.

---

## Appendix A: Reproducibility

```bash
git clone https://github.com/keshav55/cure
cd cure
pip install mhcflurry biopython openpyxl
mhcflurry-downloads fetch models_class1 models_class1_presentation

# Run pipeline with example mutations
python pipeline.py

# Run clinical validation
python validate.py

# Run TESLA validation
python validate_tesla.py
```

## Appendix B: Daemon Statistics

| Metric | Value |
|--------|-------|
| Total experiments run | 1,048 |
| Total cycles | 573 |
| Improvements accepted | 67 |
| Keep rate | 6.4% |
| Experiments registered | 12 |
| Experiments at ceiling | 7 |
| Parallel batch size | 3 |
| LLM cost | $0 |
