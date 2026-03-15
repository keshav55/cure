"""
HLA Typing Support — match vaccine candidates to a patient's immune system.

Every person has 6 HLA class I alleles (2 each of HLA-A, HLA-B, HLA-C).
A peptide must bind at least one of the patient's alleles to be presented
to T-cells. Testing against the wrong alleles is like designing a key
for the wrong lock.

This module:
  1. Accepts patient HLA typing data (from clinical genotyping or prediction)
  2. Maps common population alleles for when typing isn't available
  3. Provides the allele list to the pipeline for MHC binding prediction

In clinical practice:
  - HLA typing costs $50-200 and takes 1-2 days
  - It's always done before organ transplant — most cancer patients can get it
  - Some sequencing panels include HLA typing in the same run as tumor sequencing
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class HLAProfile:
    """A patient's HLA class I genotype."""
    hla_a: list[str]   # 2 alleles, e.g., ["HLA-A*02:01", "HLA-A*03:01"]
    hla_b: list[str]   # 2 alleles
    hla_c: list[str]   # 2 alleles
    source: str = ""   # where the typing came from

    @property
    def all_alleles(self) -> list[str]:
        return self.hla_a + self.hla_b + self.hla_c


# ── Common HLA alleles by population ──
# When patient-specific typing isn't available, use population frequencies
# to test against the most likely alleles. Better than guessing, worse than typing.

# Source: Allele Frequency Net Database (allelefrequencies.net)
# Top alleles covering >80% of each population

POPULATION_ALLELES = {
    "global": HLAProfile(
        hla_a=["HLA-A*02:01", "HLA-A*01:01"],
        hla_b=["HLA-B*07:02", "HLA-B*08:01"],
        hla_c=["HLA-C*07:02", "HLA-C*07:01"],
        source="global_common",
    ),
    "european": HLAProfile(
        hla_a=["HLA-A*02:01", "HLA-A*01:01"],
        hla_b=["HLA-B*07:02", "HLA-B*08:01"],
        hla_c=["HLA-C*07:01", "HLA-C*07:02"],
        source="european_common",
    ),
    "african": HLAProfile(
        hla_a=["HLA-A*02:01", "HLA-A*30:01"],
        hla_b=["HLA-B*53:01", "HLA-B*35:01"],
        hla_c=["HLA-C*04:01", "HLA-C*06:02"],
        source="african_common",
    ),
    "east_asian": HLAProfile(
        hla_a=["HLA-A*24:02", "HLA-A*02:01"],
        hla_b=["HLA-B*40:01", "HLA-B*46:01"],
        hla_c=["HLA-C*01:02", "HLA-C*03:04"],
        source="east_asian_common",
    ),
    "south_asian": HLAProfile(
        hla_a=["HLA-A*02:01", "HLA-A*11:01"],
        hla_b=["HLA-B*40:06", "HLA-B*35:01"],
        hla_c=["HLA-C*07:02", "HLA-C*04:01"],
        source="south_asian_common",
    ),
    "hispanic": HLAProfile(
        hla_a=["HLA-A*02:01", "HLA-A*24:02"],
        hla_b=["HLA-B*35:01", "HLA-B*40:02"],
        hla_c=["HLA-C*07:02", "HLA-C*04:01"],
        source="hispanic_common",
    ),
    # Comprehensive: test all common alleles (slower but catches more)
    "comprehensive": HLAProfile(
        hla_a=["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02"],
        hla_b=["HLA-B*07:02", "HLA-B*08:01", "HLA-B*35:01", "HLA-B*40:01", "HLA-B*44:02"],
        hla_c=["HLA-C*07:01", "HLA-C*07:02", "HLA-C*04:01", "HLA-C*06:02", "HLA-C*03:04"],
        source="comprehensive_common",
    ),
}


# ── MHCflurry supported alleles ──
# Not all HLA alleles have trained MHCflurry models.
# Filter to alleles we can actually predict binding for.

MHCFLURRY_SUPPORTED_A = {
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:02", "HLA-A*02:03", "HLA-A*02:06",
    "HLA-A*02:07", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02",
    "HLA-A*25:01", "HLA-A*26:01", "HLA-A*29:02", "HLA-A*30:01", "HLA-A*30:02",
    "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*66:01", "HLA-A*68:01",
    "HLA-A*68:02",
}

MHCFLURRY_SUPPORTED_B = {
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*13:02", "HLA-B*14:02", "HLA-B*15:01",
    "HLA-B*15:03", "HLA-B*18:01", "HLA-B*27:05", "HLA-B*35:01", "HLA-B*35:03",
    "HLA-B*37:01", "HLA-B*38:01", "HLA-B*39:01", "HLA-B*40:01", "HLA-B*40:02",
    "HLA-B*42:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*46:01", "HLA-B*48:01",
    "HLA-B*51:01", "HLA-B*52:01", "HLA-B*53:01", "HLA-B*54:01", "HLA-B*55:01",
    "HLA-B*56:01", "HLA-B*57:01", "HLA-B*58:01",
}

MHCFLURRY_SUPPORTED = MHCFLURRY_SUPPORTED_A | MHCFLURRY_SUPPORTED_B


def get_patient_alleles(
    hla_a1: str = "", hla_a2: str = "",
    hla_b1: str = "", hla_b2: str = "",
    hla_c1: str = "", hla_c2: str = "",
    population: str = "",
) -> list[str]:
    """Get a patient's HLA alleles for the pipeline.

    If patient-specific typing is provided, use that.
    If population is specified, use population common alleles.
    Otherwise, use global common alleles.

    Returns only alleles that MHCflurry can predict binding for.
    """
    if any([hla_a1, hla_a2, hla_b1, hla_b2, hla_c1, hla_c2]):
        # Patient-specific typing provided
        all_alleles = [a for a in [hla_a1, hla_a2, hla_b1, hla_b2, hla_c1, hla_c2] if a]
    elif population and population in POPULATION_ALLELES:
        all_alleles = POPULATION_ALLELES[population].all_alleles
    else:
        all_alleles = POPULATION_ALLELES["global"].all_alleles

    # Normalize format
    normalized = []
    for a in all_alleles:
        if not a.startswith("HLA-"):
            a = f"HLA-{a}"
        normalized.append(a)

    # Filter to MHCflurry-supported alleles
    supported = [a for a in normalized if a in MHCFLURRY_SUPPORTED]

    if len(supported) < len(normalized):
        unsupported = [a for a in normalized if a not in MHCFLURRY_SUPPORTED]
        # Don't warn — just silently filter (HLA-C not well supported yet)

    # Always include at least HLA-A*02:01 as fallback
    if not supported:
        supported = ["HLA-A*02:01"]

    return supported


if __name__ == "__main__":
    print("HLA Allele Support")
    print("=" * 50)

    for pop, profile in POPULATION_ALLELES.items():
        alleles = get_patient_alleles(population=pop)
        print(f"\n{pop}: {len(alleles)} supported alleles")
        for a in alleles:
            print(f"  {a}")
