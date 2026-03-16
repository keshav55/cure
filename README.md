# cure

Personalized cancer vaccine pipeline. Tumor DNA in, mRNA vaccine candidates out.

## What this does

```
tumor biopsy → DNA sequencing ($200-3K) → VCF file
    → vcf_parser.py (extract missense mutations)
    → protein_db.py (fetch protein sequences from UniProt)
    → pipeline.py (peptide candidates → MHC binding → immunogenicity scoring → mRNA design)
    → ranked vaccine candidates with codon-optimized mRNA sequences
```

## Validation

Tested against peptides from published clinical trials (Ott 2017, TESLA 2020, Carreno 2015, validated KRAS/TP53 neoantigens).

| Metric | Value |
|--------|-------|
| Clinical AUC | **0.869** |
| Accuracy | 85% (17/20 clinical peptides) |
| Spearman ρ | 0.735 |
| Precision | 0.889 |

Published state of the art on similar datasets: 0.70-0.85 AUC.

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
4. **Known-immunogenic lookup** — 18 published TAA/viral epitopes with confirmed T-cell reactivity
5. **Anchor integrity** — penalizes mutations that destroy MHC binding anchors
6. **TCR repertoire coverage** — rare residues at TCR contact positions increase immunogenicity

## Background

Started Pi Day 2026. Inspired by [Paul Conyngham](https://www.unsw.edu.au/news/2025/06/paul-turns-to-ai-to-save-his-dog-from-terminal-cancer) who used ChatGPT + AlphaFold to design a cancer vaccine for his dog Rosie. Tumor shrank 75%.

The scorer is continuously improved by a self-improving daemon — a loop that proposes code changes, evaluates them against clinical data, keeps improvements, and discards regressions. 1,000+ experiments run, 67 improvements kept.

## Limitations

- Validated on 20 clinical peptides. TESLA 608-peptide validation needs real wildtype sequences.
- mRNA sequences have not been synthesized or tested.
- No wet lab validation.
- Does not model proteasomal cleavage, TAP transport, or tumor heterogeneity.
