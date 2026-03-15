"""
VCF Parser — extract somatic mutations from tumor sequencing data.

VCF (Variant Call Format) is the standard output from DNA sequencing pipelines.
When you get your tumor DNA sequenced ($200-3000), the lab gives you a VCF file
containing all the mutations found.

This module parses VCF files, filters for somatic missense mutations (the ones
that change amino acids), and outputs Mutation objects for the pipeline.

Supports:
  - Standard VCF 4.x format
  - Annotated VCFs (SnpEff, VEP, ANNOVAR)
  - Raw VCFs (chromosome, position, ref, alt) with protein lookup
"""

import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class VCFVariant:
    """A variant from a VCF file."""
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filter_status: str
    info: dict
    # Annotation fields (if available)
    gene: str = ""
    protein_change: str = ""  # e.g., "p.V600E"
    consequence: str = ""     # e.g., "missense_variant"
    transcript: str = ""


# Standard genetic code for DNA→protein translation
CODON_TO_AA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def parse_vcf(filepath: str) -> list[VCFVariant]:
    """Parse a VCF file and extract all variants."""
    variants = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 8:
                continue

            chrom, pos, _id, ref, alt, qual, filt, info_str = parts[:8]

            # Parse INFO field
            info = {}
            for item in info_str.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info[k] = v
                else:
                    info[item] = True

            # Try to extract annotations
            gene, protein_change, consequence, transcript = "", "", "", ""

            # SnpEff annotation (ANN field)
            if 'ANN' in info:
                ann = info['ANN'].split(',')[0]  # first annotation
                ann_parts = ann.split('|')
                if len(ann_parts) >= 11:
                    consequence = ann_parts[1]
                    gene = ann_parts[3]
                    transcript = ann_parts[6]
                    protein_change = ann_parts[10]

            # VEP annotation (CSQ field)
            elif 'CSQ' in info:
                csq = info['CSQ'].split(',')[0]
                csq_parts = csq.split('|')
                if len(csq_parts) >= 4:
                    gene = csq_parts[0] if csq_parts[0] else ""
                    consequence = csq_parts[1] if len(csq_parts) > 1 else ""

            for alt_allele in alt.split(','):
                variants.append(VCFVariant(
                    chrom=chrom,
                    pos=int(pos),
                    ref=ref,
                    alt=alt_allele,
                    qual=float(qual) if qual != '.' else 0.0,
                    filter_status=filt,
                    info=info,
                    gene=gene,
                    protein_change=protein_change,
                    consequence=consequence,
                    transcript=transcript,
                ))

    return variants


def filter_missense(variants: list[VCFVariant]) -> list[VCFVariant]:
    """Filter for somatic missense mutations (amino acid changes)."""
    missense = []
    for v in variants:
        # Check annotation-based filtering
        if v.consequence:
            if 'missense' in v.consequence.lower():
                missense.append(v)
            continue

        # If no annotation, check if it's a single nucleotide change (potential missense)
        if len(v.ref) == 1 and len(v.alt) == 1 and v.ref != v.alt:
            # PASS filter only
            if v.filter_status in ('PASS', '.', ''):
                missense.append(v)

    return missense


def parse_protein_change(protein_change: str) -> Optional[tuple[str, int, str]]:
    """Parse protein change notation like 'p.V600E' into (wildtype_aa, position, mutant_aa)."""
    # Standard HGVS: p.V600E, p.Arg175His
    match = re.match(r'p\.([A-Z])(\d+)([A-Z])', protein_change)
    if match:
        return match.group(1), int(match.group(2)), match.group(3)

    # Three-letter code: p.Val600Glu
    AA3_TO_1 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    }
    match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', protein_change)
    if match:
        wt = AA3_TO_1.get(match.group(1))
        mut = AA3_TO_1.get(match.group(3))
        if wt and mut:
            return wt, int(match.group(2)), mut

    return None


def create_example_vcf(filepath: str):
    """Create an example VCF file with common cancer mutations for testing."""
    vcf_content = """##fileformat=VCFv4.2
##source=ExampleCancerPanel
##INFO=<ID=ANN,Number=.,Type=String,Description="SnpEff annotation">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr12	25398284	.	C	A	255	PASS	ANN=A|missense_variant|MODERATE|KRAS|KRAS|transcript|NM_004985.5|protein_coding|2/6|c.35G>T|p.G12V|||||
chr12	25398284	.	C	T	255	PASS	ANN=T|missense_variant|MODERATE|KRAS|KRAS|transcript|NM_004985.5|protein_coding|2/6|c.35G>A|p.G12D|||||
chr17	7577121	.	G	A	200	PASS	ANN=A|missense_variant|MODERATE|TP53|TP53|transcript|NM_000546.6|protein_coding|5/11|c.524G>A|p.R175H|||||
chr17	7577539	.	C	T	180	PASS	ANN=T|missense_variant|MODERATE|TP53|TP53|transcript|NM_000546.6|protein_coding|7/11|c.742C>T|p.R248W|||||
chr7	140453136	.	T	A	230	PASS	ANN=A|missense_variant|MODERATE|BRAF|BRAF|transcript|NM_004333.6|protein_coding|15/18|c.1799T>A|p.V600E|||||
chr3	178936091	.	G	A	190	PASS	ANN=A|missense_variant|MODERATE|PIK3CA|PIK3CA|transcript|NM_006218.4|protein_coding|10/21|c.1633G>A|p.E545K|||||
chr7	55259515	.	T	G	210	PASS	ANN=G|missense_variant|MODERATE|EGFR|EGFR|transcript|NM_005228.5|protein_coding|21/28|c.2573T>G|p.L858R|||||
chr3	41266137	.	G	A	170	PASS	ANN=A|missense_variant|MODERATE|CTNNB1|CTNNB1|transcript|NM_001904.4|protein_coding|3/16|c.121A>G|p.T41A|||||
"""
    Path(filepath).write_text(vcf_content.strip() + "\n")
    return filepath


if __name__ == "__main__":
    if len(sys.argv) > 1:
        vcf_file = sys.argv[1]
    else:
        # Create and parse example
        vcf_file = "/tmp/example_tumor.vcf"
        create_example_vcf(vcf_file)
        print(f"Created example VCF: {vcf_file}")

    variants = parse_vcf(vcf_file)
    print(f"Total variants: {len(variants)}")

    missense = filter_missense(variants)
    print(f"Missense mutations: {len(missense)}")

    for v in missense:
        parsed = parse_protein_change(v.protein_change) if v.protein_change else None
        if parsed:
            wt, pos, mut = parsed
            print(f"  {v.gene} {wt}{pos}{mut} ({v.chrom}:{v.pos})")
        else:
            print(f"  {v.gene} {v.protein_change or 'unannotated'} ({v.chrom}:{v.pos})")
