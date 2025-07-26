from pathlib import Path
import urllib.request
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
import gzip
import shutil
from goatools.go_enrichment import GOEnrichmentStudy
from collections import OrderedDict


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
G2G_COMPRESSED = ensure_file(G2G_URL, "gene2go.gz")

def ensure_uncompressed_g2g(path: Path) -> Path:
    """Return uncompressed gene2go file, extracting if needed."""
    target = path.with_suffix("")
    if target.exists():
        return target
    try:
        with gzip.open(path, "rt") as fin, open(target, "w") as fout:
            shutil.copyfileobj(fin, fout)
    except OSError as exc:
        raise RuntimeError(f"Failed to decompress {path}: {exc}")
    return target

G2G = ensure_uncompressed_g2g(G2G_COMPRESSED)

oboDAG = GODag(str(OBO))

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
def is_enriched(genes: list[str], goea, p_cut: float) -> bool:
    if len(genes)<2:
        return False
    results = goea.run_study(genes)
    return any(r.p_fdr_bh < p_cut for r in results)
def go_assessment(
    taxid: int,
    biclusters: list[list[str]],
    universe: set[int],
    p_vals=(0.05, 0.01, 0.001),
) -> OrderedDict:
    goea = init_goea(taxid,universe)
    enriched_flags = {th: 0 for th in p_vals}
    for bic in biclusters:
        for th in p_vals:
            if is_enriched(bic, goea, th):
                enriched_flags[th] += 1
    total = len(biclusters)
    return OrderedDict((th, enriched_flags[th] / total) for th in p_vals)
