from pathlib import Path
import urllib.request
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go

# Canonical GO data sources
OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
G2G_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"

def ensure_file(url: str, local_name: str) -> Path:
    """
    Make sure *local_name* exists in the current folder.
    If it doesn’t, download it from *url*.
    """
    path = Path(local_name)
    if path.exists():
        return path
    print(f"Downloading {local_name} …")
    urllib.request.urlretrieve(url, path)
    return path
OBO = ensure_file(OBO_URL, "go-basic.obo")
G2G = ensure_file(G2G_URL, "gene2go.gz")

oboDAG = Godag(str(OBO))

def get_associations(taxid:int):
    assoc=read_ncbi_gene2go(str(G2G), taxids=[taxid])
    if not assoc:
        raise ValueError("No GO annotations—check taxid or ID type")
    return assoc
human_assoc=get_associations(9606)
print(len(human_assoc))
print(human_assoc[7157])

