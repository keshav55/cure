# cure/

Personalized cancer vaccine pipeline. From tumor DNA to mRNA vaccine candidates.

Started Pi Day 2026. Inspired by Paul Conyngham who used ChatGPT + AlphaFold
to design a cancer vaccine for his dog Rosie. Tumor shrank 75%.

## Pipeline

```
tumor biopsy → DNA sequencing ($200-3K) → VCF file
    → vcf_parser.py (extract missense mutations)
    → protein_db.py (fetch full protein sequences from UniProt)
    → pipeline.py (generate candidates → MHCflurry binding → immunogenicity scoring → mRNA design)
    → ranked vaccine candidates with codon-optimized mRNA sequences
```

## Components

| File | What | Status |
|------|------|--------|
| `pipeline.py` | End-to-end: mutations → vaccine candidates | Working |
| `vcf_parser.py` | Parse tumor VCF files, extract missense mutations | Working |
| `protein_db.py` | Fetch protein sequences from UniProt | Building |
| `hla_typing.py` | Patient HLA allele support | Planned |
| `multi_epitope.py` | Optimal peptide combination selection | Planned |

## Daemon experiments (in experiments/)

The autoresearch daemon continuously optimizes these components:

| Experiment | What it optimizes | Score |
|-----------|-------------------|-------|
| `neoantigen-scorer` | Peptide immunogenicity prediction | F1=0.814 |
| `mrna-codon-optimizer` | mRNA codon selection | 0.84 |

## Dependencies

```bash
pip install mhcflurry biopython
mhcflurry-downloads fetch models_class1 models_class1_presentation
```

## Usage

```bash
# With example cancer mutations
python cure/pipeline.py

# With your tumor VCF file
python cure/pipeline.py --mutations your_mutations.json

# Parse a VCF from sequencing
python cure/vcf_parser.py tumor_sample.vcf
```
