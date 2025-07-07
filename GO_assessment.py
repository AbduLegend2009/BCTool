from pathlib import Path
import urllib.request
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy
from collections import OrderedDict
import research_gui


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

oboDAG = GODag(OBO)

def get_associations(taxid:int):
    assoc=read_ncbi_gene2go(str(G2G), taxids=[taxid])
    if not assoc:
        raise ValueError("No GO annotations—check taxid or ID type")
    return assoc

def init_goea(taxid, universe):
    return GOEnrichmentStudy(
        list(universe),
        get_associations(taxid),
        oboDAG,
        methods=["fdr_bh"],
        alpha=1.0,
        propagate_counts=True
    )
def is_enriched(genes, goea, p_cut):
    if len(genes)<2:
        return False
    results = goea.run_study(genes)
    return any(r.p_fdr_bh < p_cut for r in results)
def go_assessment(taxid: int,
                   biclusters: list[list[int]], 
                   universe: set[int], 
                   p_vals=(0.05,0.01,0.001),
                   OBO,
                   G2G):
    goea = init_goea(taxid,universe)
    enriched_flags = {th: 0 for th in p_vals}
    for bic in biclusters:
        for th in p_vals:
            if is_enriched(bic, goea, th):
                enriched_flags[th] += 1
    total = len(biclusters)
    return OrderedDict((th, enriched_flags[th] / total) for th in pvals)