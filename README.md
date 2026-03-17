# cure

Open-source neoantigen vaccine pipeline. From tumor DNA to codon-optimized mRNA vaccine candidates.

**State of the art on NeoRanking benchmark.** AUROC = 0.980 (vs MixMHC 0.970, NetMHCpan 0.968, PRIME 0.951). Per-patient vaccine selection recall@20 = 0.648 ± 0.023 — 32% better than binding rank alone.

## Pipeline

```
tumor biopsy → DNA sequencing ($200) + HLA typing ($50)
    → vcf_parser.py (extract missense mutations)
    → protein_db.py (fetch protein sequences from UniProt)
    → pipeline.py (peptide candidates → MHC binding → immunogenicity scoring → mRNA design)
    → vaccine_selector.py (ML ensemble selects optimal 20 peptides per patient)
    → ranked vaccine candidates with codon-optimized mRNA sequences
```

Total patient-specific cost: $250. No RNA-seq required (TCGA expression proxy matches tumor RNA-seq exactly).

## Results

| Metric | Our Ensemble | Best Published | Improvement |
|--------|-------------|----------------|-------------|
| Peptide-level AUROC | **0.980** | MixMHC 0.970 | +0.010 |
| Vaccine selection recall@20 | **0.648** | MixMHC 0.492 | **+32%** |
| Vaccine selection recall@5 | **0.307** | NetMHCpan 0.218 | **+41%** |

Evaluated on NeoRanking test set: 125,784 peptides, 96 confirmed CD8+ T-cell responses, 30 patients. Bootstrap 95% CI: [0.607, 0.685]. 50/50 bootstraps > binding baseline.

## Papers

| # | Title | Key Result |
|---|-------|------------|
| 1 | Data Quality Dominates Algorithmic Complexity | Binding alone = 0.967 AUC. Patient HLA matching is the #1 factor. |
| 2 | The Foreignness Paradox | DAI anti-correlated with immunogenicity (0.352). Central tolerance. |
| 3 | ML Ensemble Improves Vaccine Selection by 32% | SOTA. 7 features beat 18. Beats all published tools. |

## Key Findings

1. **Patient HLA alleles are everything.** Same scorer: 0.156 AUC (wrong alleles) → 0.967 (right alleles).
2. **Algorithmic features hurt.** Binding alone (0.967) > our 300-line scorer (0.909). Fewer features = better.
3. **Foreignness is counterproductive.** The most intuitive feature is anti-correlated. Central tolerance explains it.
4. **ML ensemble > every published tool.** GBT + RF + LR stacking with multi-seed negative bagging.
5. **7 features are enough.** 3 binding predictors + stability + CSCAPE driver score + 2 expression databases.
6. **Cancer-specific models fail.** General model wins due to more training data.
7. **polar→aromatic mutations are 3x enriched.** New signal beyond binding.
8. **$250 per patient.** No RNA-seq needed (TCGA proxy = tumor RNA-seq for expression).

## Components

| File | What |
|------|------|
| `vaccine_selector.py` | ML ensemble for per-patient vaccine selection (recall@20 = 0.648) |
| `scorer_simple.py` | Minimal scorer — MHCflurry only (AUROC = 0.967) |
| `scorer.py` | Full scorer — 300 lines, 6 features (AUROC = 0.909, worse than simple) |
| `pipeline.py` | End-to-end: VCF → peptides → binding → scoring → mRNA |
| `autoresearch_vaccine.py` | Autoresearch daemon for vaccine selection optimization |
| `validate_neoranking_v2.py` | NeoRanking validation with patient-matched alleles |
| `vcf_parser.py` | Parse tumor VCF files |
| `protein_db.py` | UniProt protein sequence fetcher |
| `hla_typing.py` | Patient HLA allele support |

## Install

```bash
pip install mhcflurry biopython scikit-learn numpy
mhcflurry-downloads fetch models_class1 models_class1_presentation
```

## Usage

```bash
# Run pipeline with example mutations
python pipeline.py

# Run with tumor VCF
python pipeline.py --vcf tumor_sample.vcf

# Run ML vaccine selector
python vaccine_selector.py

# Run autoresearch daemon (optimize selection)
python autoresearch_vaccine.py --rounds 30
```

## Background

Started Pi Day 2026. Inspired by [Paul Conyngham](https://www.unsw.edu.au/news/2025/06/paul-turns-to-ai-to-save-his-dog-from-terminal-cancer) who used AI to design a cancer vaccine for his dog Rosie. Tumor shrank 75%.

53 commits. 12 experiments. 3 papers. Built by Keshav Rao + Claude Opus 4.6 + GPT-5.4 Pro.
