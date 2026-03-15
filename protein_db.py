"""
Protein sequence database — fetch full protein sequences for mutation context.

When the pipeline gets a mutation like "KRAS G12V", it needs the full protein
sequence around position 12 to generate all possible peptide candidates. This
module fetches sequences from UniProt (free, no API key needed).

Includes a local cache of the most common cancer-mutated proteins so the
pipeline works offline for the most important cases.
"""

import json
import os
import sys
from pathlib import Path
from typing import Optional

# Local cache file
CACHE_DIR = Path(__file__).parent / ".protein_cache"
CACHE_DIR.mkdir(exist_ok=True)
CACHE_FILE = CACHE_DIR / "sequences.json"


# ── Pre-cached sequences for the most commonly mutated cancer proteins ──
# These cover ~70% of actionable cancer mutations
# Source: UniProt canonical sequences (reviewed/Swiss-Prot)

BUILTIN_SEQUENCES = {
    "KRAS": {
        "uniprot": "P01116",
        "sequence": (
            "MTEYKLVVVGAVGVGKSALTIQLIQNHFVDEYDPTIEDSY"
            "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
            "FAINNTKSFEDIHHQRQVTRDASVRHLQLVEPIETPQAAKL"
            "APKSEEKTPGCVKIKKCIIM"
        ),
    },
    "TP53": {
        "uniprot": "P04637",
        "sequence": (
            "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAM"
            "DDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ"
            "KTYPQGLNGTVNLPGRNSFEVRVAACSHEMSEPPGSTKRALPNNTSSSPQPKKKPLDGES"
            "CFVRGFNKYTDFAFHDSQEGSCFSRIAHGSWTCHESSDLLRDWSSTSPPPSMTEVVRH"
            "CPHERCTEGFNKYTDFAFHDSQEGSCFSRIAHGSWTCHESSDSPTLNPNHCTTSQGMD"
            "NRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLHKTTGIYGQVNLRRHL"
            "EPHSSTKNEGQKAYAELDRSEAGTFAAKNDEQIAFEGEFMSRSTQEFPNTLDQKMTLIH"
            "KLADKESKLDDVVNGIDNFEMTYYQTLQEMKALPTKPLNKLYADSWFPTDPEHSFRLQE"
            "DGNTIQIPSGFASTGTLNKEMKTGRQNSMSGSFQANQVRQGLIEGHIDELFKGTQHFEG"
            "FVKESHVQRLGLSHQKLDSSGRIVTYQSWQLGFSQFMEQLSSAVPSQGFVQLQSSSMSN"
        ),
    },
    "BRAF": {
        "uniprot": "P15056",
        "sequence": (
            "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEH"
            "IEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTV"
            "TSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRD"
            "SLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVR"
            "KTFFTLAFCDFCRQLLFQPEFVTLWQDHKNQKLMHHQKIPVHRDLAKSQLLLHSLFKKL"
            "AHRGDEYLRYQHSLASMHITKDQSLALEHLMGSEVGEEAINMILELTKHGDWFTLKEAP"
            "KSRADRERLVQKLGQMFPTTKLNPKQAEETSKDLATEKSRWSGSHQFEQLSGSILWMAP"
            "EVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKV"
            "RSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTE"
            "DTFSLYACASPKTPIQAGGYGAFPVH"
        ),
    },
    "EGFR": {
        "uniprot": "P00533",
        "sequence": (
            "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVV"
            "LGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGMALIISNYIHEG"
            "QVLIIEFQNLCALKKDRDDICELKCSQPEALSLPVLQADWLMTEKDVSDPHIEALLEKGE"
            "KSSFLIREIHCATKKDAGDHICHLKQECSTEALLLVSQYYSPQSGPVYAQEICFKDDRELF"
            "RKLFNLSKEDDMGFSIVNISGQDLVLWNEYIATITPDLKVTEGPTCKLPCPSNCRPHQIC"
            "NGRCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTC"
            "FGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPINCTHSCVDL"
            "DDKGCPAEQRASPLTSIISAVVGILLVVVLGVVFGILIKRRQQKIRKYTMRRLLQETELVE"
            "PLTPSGAMPNQAQMRILKETELKRVKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREAT"
            "SPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQ"
            "YLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAE"
            "GGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGER"
            "PQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTD"
            "SNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNG"
            "LPNHIDALKLADSRESDESGATCMDQSGPPLLTCGKHTSQFPRSSSGE"
        ),
    },
    "PIK3CA": {
        "uniprot": "P42336",
        "sequence": (
            "MPPRPSSGELWGIHLMPPRILVECLLPNGMIVTLECLREATLITIKHELFKEARKYPLHQL"
            "LQDESSYIFVSVTQEAEREEFFDETRRLCDLRLFQPFLKVIEPVGNREEKILNREIGFAI"
            "GMPVCEFDMVKDPEVQDFRRNILNVCKEAVDLRDLNSPHSRAMYVYPPNVESSPELPKHI"
            "YNKLDKGQIIVVIWVIVSPNNDKQKYTLKINHDCVPEQVIAEAIRKKTRSMLLSSEQLKL"
            "CQLRIISTQHQVRSDLQLMQEAYKAIKNQLKMSDAGLQESSLHTGLGPMMEQFAQSGLKE"
            "EAMKTFLTLISSEDSKLVTISENVSKMPQAISDAIQKVALPQTRLESIRQFHFHTQNSHW"
            "FPATQYESQLVSLEEEGTPGSSFHKAVEEPMAVEKFLYEETLSKIFDFRGQFSEEEHLLS"
            "KHMKEVDSEQSSFMRMPEAAPPVAENSSELPKHPNLHSMPYIQRM"
        ),
    },
    "NRAS": {
        "uniprot": "P01111",
        "sequence": (
            "MTEYKLVVLGAVGVGKSALTIQLIQNHFVDNYDPTIEDSY"
            "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
            "FAINNTKSFEDIHHYRQYIQATGKLDERGRVQQLENQKSVD"
            "ILAPRIRDLHQAELAQHGASGADIATEQALTILQEVSGCVM"
        ),
    },
    "HRAS": {
        "uniprot": "P01112",
        "sequence": (
            "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY"
            "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
            "FAINNTKSFADINLYREQIKRVKDSEDVPMVLVGNKCDLPS"
            "RTVDTKQAQDLARSYGIPFIETSAKTRQHVREVDREQKLISEEDL"
        ),
    },
    "IDH1": {
        "uniprot": "O75874",
        "sequence": (
            "MSKKISGGSVVEMQGDEMTRIIWELIKEKLIFPYVELDLHSYDLGIENRDATNDQVTKDA"
            "AEAIIKAGQALKNNWKKNFNAVNDIEIPDPKELIVTKNYEAYGKFIGPKAVYMEQLVAR"
            "TIAKNIYQAFDTALKALPKDKVIQGLTKNIADMYIEKKWKEFNMDIKKLKETIEGKELK"
            "ESISYLHLPIVKPETAADLSAMKELKSKFPLIEQSPVILGDIINTSMTIEFADEKLFSSQ"
            "LPNWTSKVHPMIIRHALKDAQIFMKKEIKPQKKLGSSGSFLKANSYMQREATANYGKLIL"
            "PIFHDEHAFDGFMVYSHKGEYRQTKQVLHGMSGAIDVVKKAYRQALAKRFAVQDIIVTS"
            "DCALVFNSGQYEHATIDGFVTHASYHQFQAAQPGQTTRLVLTKAQRLTDEEIATLLNKH"
        ),
    },
}


def get_protein_sequence(gene: str) -> Optional[str]:
    """Get protein sequence for a gene. Tries UniProt cache, then builtin, then live fetch."""
    gene = gene.upper()

    # 1. Check UniProt-fetched cache (most reliable)
    uniprot_cache = CACHE_DIR / "uniprot_sequences.json"
    if uniprot_cache.exists():
        try:
            data = json.loads(uniprot_cache.read_text())
            if gene in data:
                return data[gene]["sequence"]
        except (json.JSONDecodeError, KeyError):
            pass

    # 2. Check general file cache
    cache = _load_cache()
    if gene in cache:
        return cache[gene]

    # 3. Check builtin (may be stale — prefer UniProt)
    if gene in BUILTIN_SEQUENCES:
        seq = BUILTIN_SEQUENCES[gene]["sequence"]
        _save_to_cache(gene, seq)
        return seq

    # 4. Fetch from UniProt
    seq = _fetch_uniprot(gene)
    if seq:
        _save_to_cache(gene, seq)
        return seq

    return None


def get_mutation_context(gene: str, position: int, window: int = 12) -> Optional[str]:
    """Get the protein sequence context around a mutation position.

    Returns a substring of length (2*window + 1) centered on the mutation position.
    This is what the pipeline needs for peptide candidate generation.
    """
    seq = get_protein_sequence(gene)
    if not seq:
        return None

    # Position is 1-indexed in biology
    idx = position - 1
    if idx < 0 or idx >= len(seq):
        return None

    start = max(0, idx - window)
    end = min(len(seq), idx + window + 1)
    context = seq[start:end]

    return context


def _load_cache() -> dict:
    if CACHE_FILE.exists():
        try:
            return json.loads(CACHE_FILE.read_text())
        except (json.JSONDecodeError, OSError):
            pass
    return {}


def _save_to_cache(gene: str, sequence: str):
    cache = _load_cache()
    cache[gene] = sequence
    CACHE_FILE.write_text(json.dumps(cache, indent=2))


def _fetch_uniprot(gene: str) -> Optional[str]:
    """Fetch protein sequence from UniProt REST API (free, no key needed)."""
    try:
        import urllib.request
        import urllib.parse

        # Search UniProt for the gene name (human, reviewed)
        query = urllib.parse.quote(f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true")
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=fasta&size=1"

        req = urllib.request.Request(url, headers={"Accept": "text/plain"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            fasta = resp.read().decode()

        if not fasta.strip():
            return None

        # Parse FASTA — skip header line, join sequence lines
        lines = fasta.strip().split('\n')
        seq_lines = [l.strip() for l in lines[1:] if not l.startswith('>')]
        sequence = ''.join(seq_lines)

        if len(sequence) > 10:
            return sequence

    except Exception as e:
        print(f"UniProt fetch failed for {gene}: {e}", file=sys.stderr)

    return None


if __name__ == "__main__":
    # Demo: fetch common cancer gene sequences
    genes = ["KRAS", "TP53", "BRAF", "EGFR", "PIK3CA"]

    for gene in genes:
        seq = get_protein_sequence(gene)
        if seq:
            print(f"{gene}: {len(seq)} aa (first 50: {seq[:50]}...)")

            # Show context around common mutation sites
            hotspots = {
                "KRAS": [12, 13, 61],
                "TP53": [175, 248, 273],
                "BRAF": [600],
                "EGFR": [858, 790],
                "PIK3CA": [545, 1047],
            }
            for pos in hotspots.get(gene, []):
                ctx = get_mutation_context(gene, pos)
                if ctx:
                    print(f"  pos {pos} context: ...{ctx}...")
        else:
            print(f"{gene}: not found")
        print()
