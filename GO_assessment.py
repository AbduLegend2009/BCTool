# filesystem + download helpers
from pathlib import Path
import urllib.request
import gzip
import shutil

# GO ontology & annotations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go

# enrichment engine
from goatools.go_enrichment import GOEnrichmentStudy


def ensure_file(path: Path, url: str):
    if not path.exists():
        urllib.request.urlretrieve(url, str(path))



def gunzip_if_needed(gz_path: Path, out_path: Path):
    if not out_path.exists():
        with gzip.open(str(gz_path), 'rb') as fin, out_path.open('wb') as fout:
            shutil.copyfileobj(fin, fout)


def init_goea(taxid: int, universe: set[int]) -> GOEnrichmentStudy:
    """
    1. Parse the GO ontology (only once per run).
    2. Load NCBI gene2go for the given taxid.
    3. Build and return a GOEnrichmentStudy.
    """
    # 1. Ontology graph
    dag = GODag(str(OBO_PATH))

    # 2. Associations: gene→GO terms for our organism
    assoc = read_ncbi_gene2go(str(G2G_TXT_PATH), taxids=[taxid])
    if not assoc:
        raise ValueError(f"No GO annotations found for taxid={taxid}")

    # 3. Enrichment engine: universe + assoc + graph
    goea = GOEnrichmentStudy(
        list(universe),      # background gene list
        assoc,               # geneid → goids
        dag,                 # ontology
        methods=["fdr_bh"],  # apply Benjamini–Hochberg FDR
    )
    return goea




OBO_URL  = "http://purl.obolibrary.org/obo/go/go-basic.obo"
G2G_URL  = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"



OBO_PATH        = Path("data/go-basic.obo")
G2G_GZ_PATH     = Path("data/gene2go.gz")
G2G_TXT_PATH    = Path("data/gene2go.txt")




# 1. download if needed
ensure_file(OBO_PATH, OBO_URL)
ensure_file(G2G_GZ_PATH, G2G_URL)


# 2. decompress if needed
gunzip_if_needed(G2G_GZ_PATH, G2G_TXT_PATH)