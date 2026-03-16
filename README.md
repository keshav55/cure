# cure

A repo that leverages the best AI tools, models, and research to solve problems in health. We will find a cure, no matter what.

## What this does

```
tumor biopsy → DNA sequencing ($200-3K) → VCF file
    → vcf_parser.py (extract missense mutations)
    → protein_db.py (fetch protein sequences from UniProt)
    → pipeline.py (peptide candidates → MHC binding → immunogenicity scoring → mRNA design)
    → ranked vaccine candidates with codon-optimized mRNA sequences
```

## Validation

Tested against peptides from published clinical trials (Ott 2017, TESLA 2020, Carreno 2015).

| Dataset | N | AUC | 95% CI | Significant? |
|---------|---|-----|--------|-------------|
| TESLA | 608 | 0.595 | [0.50, 0.69] | Yes (p < 0.05) |
| Clinical | 20 | 0.596 | [0.34, 0.85] | No (N too small) |

Published state of the art: 0.74-0.81 AUC (DeepImmuno, PRIME, BigMHC).

**Gap analysis:** We're missing tumor RNA expression (the #2 predictor after binding) and real wildtype sequences for TESLA. Adding these is the next step.

## Components

| File | What |
|------|------|
| `pipeline.py` | End-to-end: mutations → peptide candidates → MHC binding → immunogenicity → mRNA |
| `scorer.py` | Neoantigen immunogenicity scorer (MHCflurry + BLOSUM80 + dissimilarity) |
| `vcf_parser.py` | Parse tumor VCF files from DNA sequencing |
| `protein_db.py` | Fetch protein sequences from UniProt + local cache |
| `hla_typing.py` | Patient HLA allele support (6 populations, 49 alleles) |
| `validate.py` | Clinical validation against published trial data |
| `validate_tesla.py` | TESLA validation (608 peptides from Wells 2020 Cell) |

## Install

```bash
pip install mhcflurry biopython openpyxl
mhcflurry-downloads fetch models_class1 models_class1_presentation
```

## Usage

```bash
# Run with example cancer mutations (KRAS, TP53, BRAF, EGFR, PIK3CA)
python pipeline.py

# Run with your tumor VCF file
python pipeline.py --vcf tumor_sample.vcf

# Run clinical validation
python validate.py
```

## How the scorer works

1. **MHCflurry** (neural network) — predicts peptide-MHC binding across multiple HLA alleles
2. **BLOSUM80 foreignness** — evolutionary distance of the mutation, weighted by conservation
3. **Mutation dissimilarity** — cross-group amino acid changes at TCR-contact positions score highest
4. **Known-immunogenic lookup** — 14 published TAA/viral epitopes (no test-set overlap)
5. **Anchor integrity** — penalizes mutations that destroy MHC binding anchors
6. **TCR repertoire coverage** — rare residues at TCR contact positions increase immunogenicity

## What's novel

1. **Self-improving scorer.** Continuously optimized by an autoresearch daemon (1,000+ experiments, 67 improvements). No published neoantigen tool does this.
2. **End-to-end pipeline.** VCF → protein → peptides → binding → scoring → mRNA. One repo.
3. **BLOSUM80 conservation-weighted foreignness.** Biologically grounded mutation scoring.
4. **mRNA codon optimizer.** Beam search + simulated annealing for translation efficiency, GC balance, CpG suppression.

## Background

Started Pi Day 2026. Inspired by [Paul Conyngham](https://www.unsw.edu.au/news/2025/06/paul-turns-to-ai-to-save-his-dog-from-terminal-cancer) who used ChatGPT + AlphaFold to design a cancer vaccine for his dog Rosie. Tumor shrank 75%.

## Limitations

- AUC = 0.595 (below SOTA 0.74-0.81). Missing tumor expression data and real wildtype sequences.
- Validated on 608 TESLA peptides (statistically significant) + 20 clinical peptides (too small for significance).
- mRNA sequences have not been synthesized or tested.
- No wet lab validation.
- Does not model proteasomal cleavage, TAP transport, or tumor heterogeneity.

## Next steps

1. Download NeoRanking dataset (figshare) — has real wildtype sequences + expression data
2. Add `TUMOR_ABUNDANCE` as a feature — the #2 predictor
3. Train on NeoRanking, test on TESLA — proper train/test split
4. Compare to DeepImmuno, PRIME, BigMHC baselines
5. Write it up
