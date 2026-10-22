"""
Generate a small, connected knowledge-graph sample for a species.

Implements the "ego-network closure" algorithm from
docs/knowledge-graph/sample-generation.md: start from a curated set of anchor genes,
expand by real cross-references (STRING PPI/coexpression, coxpresdb, TFLink)
up to a node budget, then filter every real source file down to rows that
touch the resulting ID set. Output mirrors the directory layout used by
config/<species>/<species>_adapters_config.yaml so the generated files can be
pointed at directly from a sample adapters config with input_dir: samples/<species>.

Runs standalone against any pre-existing local mirror of the real data
(--input-dir), or directly against the output of
biocypher_dataset_downloader/download_data.py (same --input-dir, since that
tool writes one subdirectory per source_id, matching the layout expected
here). Source file names embed a release version (e.g. "_fb_2026_02",
"BDGP6.54.62", "gencode.v49") that changes over time as each species'
<species>_data_source_config.yaml gets bumped to a newer release — so every
real input file is located via a glob pattern (resolve_source_file) rather
than a hardcoded exact name, and every output file is written under a
stable, version-free name so the generated sample adapters config never
needs to change when upstream releases move on.

Two backbone strategies for closing gene -> transcript/protein/uniprot/entrez:
- dmel: FlyBase's own precomputed mapping tables (fbgn_fbtr_fbpp_expanded,
  fbgn_uniprot) — build_global_maps(). dmel also gets extra FlyBase-specific
  adapters (gene_group, disease_model, genotype_phenotype, allele,
  physical_interaction, gene_genetic_interaction, orthology/paralogy) that
  don't exist for any other species.
- hsa/mmu/rno/cel: no such species-specific mapping table exists, so the
  backbone is derived from generic sources every species has —
  build_generic_backbone_maps(): gene<->transcript from the GTF's own
  "transcript" rows, gene<->protein<->uniprot from the UniProt .dat file's
  "DR   Ensembl; ENST...; ENSP...; ENSG...;" cross-reference lines, and
  gene<->entrez via biocypher_metta.processors.EntrezEnsemblProcessor (the
  same cache the real adapters use at KG-build time, so the closure stays
  consistent with what TFLink/coxpresdb resolve to later).

GAF annotation files key their DB_Object_ID differently per species/GO
consortium member, not always by the primary gene ID: hsa's goa_human.gaf.gz
uses UniProt accessions, mmu's mgi.gaf.gz uses MGI IDs, rno's rgd.gaf.gz uses
RGD IDs (bare, no prefix), cel's wb.gaf.gz and dmel's own GAF use the native
gene ID directly. `gaf_id_space` per species picks which closure ID set to
filter GAF rows against.

Usage:
    python scripts/generate_connected_sample.py --species dmel
    python scripts/generate_connected_sample.py --species hsa --input-dir data/hsa
"""

import argparse
import csv
import functools
import gzip
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

# Make the repo root importable regardless of how this script is invoked
# (directly as `python scripts/generate_connected_sample.py`, or imported by
# biocypher_dataset_downloader/download_data.py after a download).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from biocypher_metta.adapters import Adapter
from biocypher_metta.processors import EntrezEnsemblProcessor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

SPECIES = {
    "dmel": {
        "taxon_id": 7227,
        "is_flybase": True,
        "anchor_genes_file": "config/dmel/dmel_anchor_genes.yaml",
        "anchor_id_key": "fbgn",
        "default_output_dir": "samples/dmel",
        "gene_info_file": "aux_files/dmel/Drosophila_melanogaster.gene_info.gz",
        "string_taxon_prefix": "7227.",
        "gaf_id_space": "gene",
    },
    "hsa": {
        "taxon_id": 9606,
        "is_flybase": False,
        "anchor_genes_file": "config/hsa/hsa_anchor_genes.yaml",
        "anchor_id_key": "ensg",
        "default_output_dir": "samples/hsa",
        "string_taxon_prefix": "9606.",
        "gaf_id_space": "uniprot",
        "uniprot_dr_db": "Ensembl",
        "uniprot_gene_id_marker": None,
    },
    "mmu": {
        "taxon_id": 10090,
        "is_flybase": False,
        "anchor_genes_file": "config/mmu/mmu_anchor_genes.yaml",
        "anchor_id_key": "ensmusg",
        "default_output_dir": "samples/mmu",
        "string_taxon_prefix": "10090.",
        "gaf_id_space": "mgi",
        "uniprot_dr_db": "Ensembl",
        "uniprot_gene_id_marker": None,
    },
    "rno": {
        "taxon_id": 10116,
        "is_flybase": False,
        "anchor_genes_file": "config/rno/rno_anchor_genes.yaml",
        "anchor_id_key": "ensrnog",
        "default_output_dir": "samples/rno",
        "string_taxon_prefix": "10116.",
        "gaf_id_space": "rgd",
        "uniprot_dr_db": "Ensembl",
        "uniprot_gene_id_marker": None,
    },
    "cel": {
        "taxon_id": 6239,
        "is_flybase": False,
        "anchor_genes_file": "config/cel/cel_anchor_genes.yaml",
        "anchor_id_key": "wbgene",
        "default_output_dir": "samples/cel",
        "string_taxon_prefix": "6239.",
        "gaf_id_space": "gene",
        # uniprot_sprot_invertebrates.dat.gz is shared across many invertebrate
        # species (Ciona, insects, mollusks, ...), and cel's own DR cross-refs
        # use "EnsemblMetazoa" (not "Ensembl") with a WormBase locus ID as both
        # transcript and protein field, WBGene as the gene field — same DR
        # shape the real UniprotProteinAdapter matches via
        # TRANSLATION_CONDITION_MAP[6239]. The WBGene marker also guards
        # against picking up other species' EnsemblMetazoa DR lines from the
        # same shared file.
        "uniprot_dr_db": "EnsemblMetazoa",
        "uniprot_gene_id_marker": "WBGene",
        # WormBase locus-style transcript/protein IDs (e.g. "T11F9.4a.1") use
        # dots as part of the identifier itself, not as a GENCODE/Ensembl
        # version separator — stripping "everything after the first dot"
        # would truncate them to garbage ("T11F9"). WBGene gene IDs have no
        # dots either way, so this only affects transcript/protein parsing.
        "strip_dot_version": False,
    },
}

FBGN_RE = re.compile(r"FBgn\d+")
FBAL_RE = re.compile(r"FBal\d+")
GTF_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GTF_TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def data_source_config_path(species):
    return f"config/{species}/{species}_data_source_config.yaml"


@functools.lru_cache(maxsize=None)
def load_data_source_config(species):
    with open(data_source_config_path(species)) as f:
        return yaml.safe_load(f)


def species_input_dir(species):
    """The real, full-size source data dir this species' sample generation reads
    from by default — the `input_dir:` top-level key in that species'
    data_source_config.yaml (also where biocypher_dataset_downloader/download_data.py
    downloads to), overridable per call via --input-dir.
    """
    return load_data_source_config(species)["input_dir"]


@functools.lru_cache(maxsize=None)
def local_files_index(species):
    """Flatten every source_id's `local_files:` sub-dict in this species'
    data_source_config.yaml into name -> (pattern, stable_name, source_id).
    `pattern` is a glob relative to input_dir (absorbs release-version
    differences: fb_2026_02 vs fb_2026_07, GENCODE/Ensembl release bumps,
    STRING/TFLink point releases, ...). `stable_name` is the version-free
    relative path used when writing the filtered sample output (None if this
    source is backbone-only and never copied to samples/). `source_id` is the
    top-level data_source_config.yaml key that owns this file — used to
    auto-download it if missing (see ensure_required_files) — with no url:
    for local_files that aren't download-managed at all (e.g. dmel's
    fca2/afca atlases, hsa's PEREGRINE/motif_diff/roadmap forgedb files).
    """
    cfg = load_data_source_config(species)
    index = {}
    for source_id, source_config in cfg.items():
        if source_id in ("name", "input_dir") or not isinstance(source_config, dict):
            continue
        for name, spec in (source_config.get("local_files") or {}).items():
            index[name] = (spec["pattern"], spec.get("stable_name"), source_id)
    return index


# (species, name) -> fallback directory to search when the primary glob comes up
# empty under input_dir. "INPUT_ROOT" means "input_dir itself, no subfolder"
# (rno's uniprot_sprot_rodents.dat.gz currently sits loose at the input root
# instead of inside uniprot/). config/dmel/dmel_data_source_config.yaml's
# flybase `move_to` relocates fbal_to_fbgn out of <output_dir>/flybase/ and into
# a fixed repo-relative path independent of --input-dir.
AUX_FALLBACK_ROOTS = {
    ("dmel", "fbal_to_fbgn"): Path("aux_files/dmel"),
    ("rno", "uniprot_dat"): "INPUT_ROOT",
}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def open_maybe_gzip(path, mode="rt", **kwargs):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, mode, **kwargs)


def resolve_source_file(input_dir, species, name):
    """Glob-resolve a real input file by logical name (see local_files: in
    config/<species>/<species>_data_source_config.yaml), so release-versioned
    filenames (fb_2026_02, BDGP6.54.62, gencode.v49, v12.0, ...) don't need to
    be hardcoded. Excludes GTF's '.chr.-only-annotation' false positives are
    not an issue here since chr_patch_hapl_scaff is the wanted file; dmel's
    '.chr.' exclusion stays dmel-specific via the pattern itself. Picks the
    lexicographically last match when several releases coexist (version
    strings sort chronologically).
    """
    pattern, _, _ = local_files_index(species)[name]
    exclude = (lambda p: ".chr." not in p.name) if species == "dmel" and name == "gtf" else (lambda p: True)
    matches = sorted(p for p in input_dir.glob(pattern) if exclude(p))
    search_root = input_dir
    fallback = AUX_FALLBACK_ROOTS.get((species, name))
    if not matches and fallback is not None:
        fallback_dir = input_dir if fallback == "INPUT_ROOT" else fallback
        bare_pattern = pattern.split("/", 1)[-1] if fallback != "INPUT_ROOT" else Path(pattern).name
        fallback_matches = sorted(p for p in fallback_dir.glob(bare_pattern) if exclude(p))
        if fallback_matches:
            logger.info(
                "'%s' not found under %s — using fallback %s instead.",
                name, input_dir, fallback_dir,
            )
            matches, search_root = fallback_matches, fallback_dir
    if not matches:
        raise FileNotFoundError(f"No file matching '{pattern}' found under {input_dir} (source: {name})")
    if len(matches) > 1:
        logger.warning(
            "Multiple files match '%s' under %s: %s — using %s",
            pattern, search_root, [m.name for m in matches], matches[-1].name,
        )
    return matches[-1]


def output_path_for(output_dir, species, name):
    _, stable_name, _ = local_files_index(species)[name]
    if stable_name is None:
        raise ValueError(f"Source '{name}' for species '{species}' is backbone-only, has no output path")
    return output_dir / stable_name


def _file_present(input_dir, species, name):
    try:
        resolve_source_file(input_dir, species, name)
        return True
    except FileNotFoundError:
        return False


_DOWNLOAD_KEYS = ("url", "directories", "zip_extract", "bucket", "bed_url", "ep_base_url")


def _is_download_managed(species, source_id):
    """Whether this source_id's data_source_config.yaml entry actually has a
    download recipe (any of DownloadManager.process_download's recognised
    keys), as opposed to being a local_files:-only placeholder for a file
    that isn't download-managed at all (dmel's fca2/afca atlases, hsa's
    PEREGRINE/motif_diff/roadmap forgedb files — real, pre-existing gaps).
    """
    source_config = load_data_source_config(species).get(source_id, {})
    return any(k in source_config for k in _DOWNLOAD_KEYS)


def ensure_required_files(species, input_dir, interactive=True):
    """Check, up front, that every real source file this species' sample
    generation will read is present under input_dir. If any are missing, list
    them (flagging which ones can be fetched automatically via
    biocypher_dataset_downloader vs. which need to be obtained manually) and
    ask the user whether to download the missing ones now. Returns False if
    the user declines, so the caller can abort before doing any processing.

    interactive=False (used by download_data.py's post-download bridge, which
    always runs right after its own download_all() already ran) skips the
    prompt entirely: known-missing files that aren't download-managed at all
    (dmel's fca2/afca, hsa's PEREGRINE/motif_diff/roadmap forgedb files, ...)
    would otherwise block that automated call on stdin every single run.
    """
    data_source_config_path_ = data_source_config_path(species)

    logger.info("Checking required source files for '%s' under %s...", species, input_dir)
    missing = []
    for name, (pattern, _, source_id) in local_files_index(species).items():
        if not _file_present(input_dir, species, name):
            missing.append((name, pattern, source_id if _is_download_managed(species, source_id) else None))

    if not missing:
        logger.info("All required source files are present.")
        return True

    if not interactive:
        logger.warning(
            "Missing (non-interactive call, continuing): %s",
            ", ".join(sorted(name for name, _, _ in missing)),
        )
        return True

    print(f"\n{len(missing)} required source file(s) missing under {input_dir}:")
    downloadable_ids = set()
    for name, pattern, source_id in sorted(missing):
        if source_id:
            print(f"  - [{name}] {pattern}  (downloadable: source '{source_id}')")
            downloadable_ids.add(source_id)
        else:
            print(f"  - [{name}] {pattern}  (no known download source — obtain manually)")

    if downloadable_ids:
        prompt = (
            f"\n{len(downloadable_ids)} download source(s) ({', '.join(sorted(downloadable_ids))}) "
            f"can be fetched via {data_source_config_path_} to cover {sum(1 for *_, s in missing if s)} "
            f"of the {len(missing)} missing file(s). Download them now and continue? [y/N]: "
        )
    else:
        prompt = "\nNone of the missing files can be downloaded automatically. Continue anyway? [y/N]: "

    reply = input(prompt).strip().lower()
    if reply not in ("y", "yes"):
        print("Aborting: required source files are missing.")
        return False

    failed_ids = []
    if downloadable_ids:
        from biocypher_dataset_downloader.download_manager import DownloadManager
        manager = DownloadManager(data_source_config_path_, input_dir)
        for source_id in sorted(downloadable_ids):
            logger.info("Downloading missing source '%s' -> %s/%s ...", source_id, input_dir, source_id)
            try:
                manager.download_source(source_id)
            except Exception as e:
                logger.error("Download of source '%s' failed (%s: %s) — skipping it, continuing with the rest.", source_id, type(e).__name__, e)
                failed_ids.append(source_id)

    still_missing = [name for name in local_files_index(species) if not _file_present(input_dir, species, name)]
    if still_missing:
        logger.warning(
            "Still missing after download: %s. Continuing anyway — the run will fail if any of these are actually needed.",
            ", ".join(sorted(still_missing)),
        )
    if failed_ids:
        logger.warning("Source(s) that failed to download: %s", ", ".join(failed_ids))
    return True


def copy_header_comments_and_filter(input_path, output_path, keep_fn, min_cols=1):
    """Stream a tab-delimited FlyBase-style file, keeping comment/blank lines
    verbatim and data lines for which keep_fn(parts) is True.

    Returns the number of data lines kept.
    """
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= min_cols and keep_fn(parts):
                fout.write(line)
                kept += 1
    return kept


# ---------------------------------------------------------------------------
# Phase 1: anchor genes
# ---------------------------------------------------------------------------

def load_anchor_genes(path, id_key):
    with open(path) as f:
        data = yaml.safe_load(f)
    return {entry[id_key] for entry in data["anchor_genes"]}


# ---------------------------------------------------------------------------
# Phase 2a: backbone closure — dmel (FlyBase precomputed tables)
# ---------------------------------------------------------------------------

def build_global_maps(input_dir, species, gene_info_path):
    """Read the real, complete FlyBase mapping files once and build forward +
    reverse lookups used throughout closure and expansion. dmel-only — other
    species use build_generic_backbone_maps() instead (no such FlyBase-style
    dedicated mapping table exists for them).
    """
    gene_to_transcript = defaultdict(set)
    gene_to_protein = defaultdict(set)
    protein_to_gene = {}
    with open_maybe_gzip(resolve_source_file(input_dir, species, "fbgn_fbtr_fbpp_expanded")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 9 and parts[2]:
                gene = parts[2]
                if parts[7]:
                    gene_to_transcript[gene].add(parts[7])
                if parts[9]:
                    gene_to_protein[gene].add(parts[9])
                    protein_to_gene[parts[9]] = gene

    gene_to_uniprot = defaultdict(set)
    uniprot_to_gene = {}
    with open_maybe_gzip(resolve_source_file(input_dir, species, "fbgn_uniprot")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 3 and parts[0] and parts[3]:
                gene_to_uniprot[parts[0]].add(parts[3])
                uniprot_to_gene[parts[3]] = parts[0]

    gene_to_allele = defaultdict(set)
    allele_to_gene = {}
    with open_maybe_gzip(resolve_source_file(input_dir, species, "fbal_to_fbgn")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2 and parts[0] and parts[2]:
                gene_to_allele[parts[2]].add(parts[0])
                allele_to_gene[parts[0]] = parts[2]

    gene_to_entrez = defaultdict(set)
    entrez_to_gene = {}
    with open_maybe_gzip(gene_info_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 5:
                for xref in parts[5].split("|"):
                    if xref.startswith("FLYBASE:"):
                        fbgn = xref.split(":", 1)[1]
                        gene_to_entrez[fbgn].add(parts[1])
                        entrez_to_gene[parts[1]] = fbgn

    return {
        "gene_to_transcript": gene_to_transcript,
        "gene_to_protein": gene_to_protein,
        "protein_to_gene": protein_to_gene,
        "gene_to_uniprot": gene_to_uniprot,
        "uniprot_to_gene": uniprot_to_gene,
        "gene_to_allele": gene_to_allele,
        "allele_to_gene": allele_to_gene,
        "gene_to_entrez": gene_to_entrez,
        "entrez_to_gene": entrez_to_gene,
    }


# ---------------------------------------------------------------------------
# Phase 2b: backbone closure — generic (hsa/mmu/rno/cel)
# ---------------------------------------------------------------------------

def _strip_dot_version(id_str, enabled=True):
    """Strip a trailing GENCODE/Ensembl-style '.<version>' suffix. Only safe
    for IDs where a dot is exclusively a version separator (ENSG...N.V). Not
    safe for WormBase locus-style IDs (e.g. "T11F9.4a.1"), where dots are part
    of the identifier itself — pass enabled=False for those (see cel's
    "strip_dot_version": False in SPECIES).
    """
    return id_str.split(".")[0] if enabled else id_str


def _parse_gtf_gene_to_transcript(gtf_path, strip_version=True):
    gene_to_transcript = defaultdict(set)
    with open_maybe_gzip(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            gm = GTF_GENE_ID_RE.search(attrs)
            tm = GTF_TRANSCRIPT_ID_RE.search(attrs)
            if gm and tm:
                gene_to_transcript[_strip_dot_version(gm.group(1), strip_version)].add(
                    _strip_dot_version(tm.group(1), strip_version)
                )
    return gene_to_transcript


def _parse_uniprot_dat_gene_protein_uniprot(uniprot_dat_path, dr_db="Ensembl", gene_id_marker=None, strip_version=True):
    """One pass over uniprot.dat: track the current record's primary accession
    (first AC line) and, for every 'DR   <dr_db>; <transcript>; <protein>;
    <gene>;' line, record gene<->protein and gene<->uniprot-accession.

    dr_db varies by species ("Ensembl" for hsa/mmu/rno, "EnsemblMetazoa" for
    cel — same DR database the real UniprotProteinAdapter matches per taxon).
    gene_id_marker, when set, requires that substring in the gene field —
    needed for uniprot_sprot_invertebrates.dat.gz, which is shared across many
    invertebrate species, to reject other species' DR lines for the same db.
    """
    gene_to_protein = defaultdict(set)
    protein_to_gene = {}
    gene_to_uniprot = defaultdict(set)
    uniprot_to_gene = {}

    dr_prefix = f"DR   {dr_db};"
    primary_accession = None
    with open_maybe_gzip(uniprot_dat_path) as f:
        for line in f:
            if line.startswith("AC   ") and primary_accession is None:
                primary_accession = line[5:].strip().rstrip(";").split(";")[0].strip()
            elif line.startswith(dr_prefix):
                fields = [x.strip().rstrip(".") for x in line[len(dr_prefix):].split(";")]
                if len(fields) >= 3 and fields[1] and fields[2]:
                    gene_field = fields[2].split()[0]  # drop trailing " [isoform]" if present
                    if gene_id_marker and gene_id_marker not in gene_field:
                        continue
                    ensp = _strip_dot_version(fields[1], strip_version)
                    ensg = _strip_dot_version(gene_field, strip_version)
                    gene_to_protein[ensg].add(ensp)
                    protein_to_gene[ensp] = ensg
                    if primary_accession:
                        gene_to_uniprot[ensg].add(primary_accession)
                        uniprot_to_gene[primary_accession] = ensg
            elif line.strip() == "//":
                primary_accession = None

    return gene_to_protein, protein_to_gene, gene_to_uniprot, uniprot_to_gene


def _build_entrez_maps_via_processor(taxon_id):
    """Load (or build, via network — same mechanism the real adapters use)
    the entrez<->ensembl cache for this species, so the closure stays
    consistent with what TFLink/coxpresdb resolve to at KG-build time.
    """
    species_info = Adapter.SPECIES_INFO[taxon_id]
    processor = EntrezEnsemblProcessor(
        ncbi_gene_info_url=species_info["ncbi_gene_info_url"],
        gencode_url=species_info["features_data_url"],
        tax_id=str(taxon_id),
        cache_dir=species_info["entrez_ensembl_cache_directory"],
        update_interval_hours=species_info["update_interval_hours"],
    )
    processor.load_or_update()

    gene_to_entrez = defaultdict(set)
    entrez_to_gene = {}
    for entrez_id, ensembl_id in processor.entrez_to_ensembl.items():
        base = ensembl_id.split(".")[0]
        gene_to_entrez[base].add(entrez_id)
        entrez_to_gene[entrez_id] = base
    return gene_to_entrez, entrez_to_gene


# NCBI gene_info dbXrefs prefix -> (prefix to match, chars to strip) for
# building the GAF-matching ID space. MGI's dbXrefs column double-wraps the
# namespace ("MGI:MGI:87854"); stripping one "MGI:" leaves the form GAF's
# DB_Object_ID actually uses ("MGI:87854"). RGD's GAF DB_Object_ID is bare
# (no "RGD:" prefix at all), so strip it entirely.
_GAF_XREF_PREFIXES = {
    "mgi": ("MGI:MGI:", "MGI:"),
    "rgd": ("RGD:", "RGD:"),
}


def build_generic_backbone_maps(input_dir, species):
    cfg = SPECIES[species]
    strip_version = cfg.get("strip_dot_version", True)
    gene_to_transcript = _parse_gtf_gene_to_transcript(
        resolve_source_file(input_dir, species, "gtf"), strip_version=strip_version
    )
    gene_to_protein, protein_to_gene, gene_to_uniprot, uniprot_to_gene = _parse_uniprot_dat_gene_protein_uniprot(
        resolve_source_file(input_dir, species, "uniprot_dat"),
        dr_db=cfg.get("uniprot_dr_db", "Ensembl"),
        gene_id_marker=cfg.get("uniprot_gene_id_marker"),
        strip_version=strip_version,
    )
    gene_to_entrez, entrez_to_gene = _build_entrez_maps_via_processor(cfg["taxon_id"])

    maps = {
        "gene_to_transcript": gene_to_transcript,
        "gene_to_protein": gene_to_protein,
        "protein_to_gene": protein_to_gene,
        "gene_to_uniprot": gene_to_uniprot,
        "uniprot_to_gene": uniprot_to_gene,
        "gene_to_entrez": gene_to_entrez,
        "entrez_to_gene": entrez_to_gene,
    }

    gaf_id_space = cfg.get("gaf_id_space", "gene")
    if gaf_id_space in _GAF_XREF_PREFIXES:
        match_prefix, strip_prefix = _GAF_XREF_PREFIXES[gaf_id_space]
        gene_to_gafid = defaultdict(set)
        with open_maybe_gzip(resolve_source_file(input_dir, species, "gene_info")) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= 5:
                    continue
                ensembl_id = gafid = None
                for xref in parts[5].split("|"):
                    if xref.startswith("Ensembl:"):
                        ensembl_id = xref.split(":", 1)[1]
                    if xref.startswith(match_prefix):
                        gafid = xref[len(strip_prefix):]
                if ensembl_id and gafid:
                    gene_to_gafid[ensembl_id].add(gafid)
        maps["gene_to_gafid"] = gene_to_gafid

    return maps


# ---------------------------------------------------------------------------
# Phase 2c: closure (species-agnostic, consumes the canonical map schema)
# ---------------------------------------------------------------------------

def close_gene_set(gene_ids, maps):
    """Derive transcript/protein/uniprot/entrez (+ allele/gafid when present
    in maps) ID sets for a gene set. Works for both build_global_maps (dmel)
    and build_generic_backbone_maps (others) output, since both use the same
    canonical key names.
    """
    transcript_ids, protein_ids, uniprot_ids, entrez_ids = set(), set(), set(), set()
    for gid in gene_ids:
        transcript_ids |= maps["gene_to_transcript"].get(gid, set())
        protein_ids |= maps["gene_to_protein"].get(gid, set())
        uniprot_ids |= maps["gene_to_uniprot"].get(gid, set())
        entrez_ids |= maps["gene_to_entrez"].get(gid, set())
    id_sets = {
        "gene": set(gene_ids),
        "transcript": transcript_ids,
        "protein": protein_ids,
        "uniprot": uniprot_ids,
        "entrez": entrez_ids,
    }
    if "gene_to_allele" in maps:
        allele_ids = set()
        for gid in gene_ids:
            allele_ids |= maps["gene_to_allele"].get(gid, set())
        id_sets["allele"] = allele_ids
    if "gene_to_gafid" in maps:
        gafid_ids = set()
        for gid in gene_ids:
            gafid_ids |= maps["gene_to_gafid"].get(gid, set())
        id_sets["gafid"] = gafid_ids
    return id_sets


# ---------------------------------------------------------------------------
# Phase 3: expansion by cross-source coverage (species-agnostic)
# ---------------------------------------------------------------------------

def find_expansion_candidates(input_dir, species, id_sets, maps, coexpression_threshold=400):
    """Return {gene_id: {source_name, ...}} for genes connected to the current
    id_sets via STRING PPI/coexpression, coxpresdb, or TFLink, but not
    already in id_sets['gene'].
    """
    partner_sources = defaultdict(set)
    protein_ids = id_sets["protein"]
    entrez_ids = id_sets["entrez"]
    uniprot_ids = id_sets["uniprot"]
    gene_ids = id_sets["gene"]
    taxon_prefix = SPECIES[species]["string_taxon_prefix"]

    def strip_taxon_prefix(protein_id):
        return protein_id[len(taxon_prefix):] if protein_id.startswith(taxon_prefix) else protein_id

    def add_partner(gene_id, source):
        if gene_id and gene_id not in gene_ids:
            partner_sources[gene_id].add(source)

    # STRING PPI
    with open_maybe_gzip(resolve_source_file(input_dir, species, "string_ppi")) as f:
        next(f, None)
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            p1 = strip_taxon_prefix(parts[0])
            p2 = strip_taxon_prefix(parts[1])
            p1_in, p2_in = p1 in protein_ids, p2 in protein_ids
            if p1_in and not p2_in:
                add_partner(maps["protein_to_gene"].get(p2), "string_ppi")
            if p2_in and not p1_in:
                add_partner(maps["protein_to_gene"].get(p1), "string_ppi")

    # STRING coexpression (detailed file, "coexpression" column thresholded)
    with open_maybe_gzip(resolve_source_file(input_dir, species, "string_coexpression")) as f:
        header = next(f).split()
        col_idx = {name: i for i, name in enumerate(header)}
        coex_i = col_idx.get("coexpression")
        for line in f:
            parts = line.split()
            if coex_i is None or len(parts) <= coex_i:
                continue
            try:
                score = int(parts[coex_i])
            except ValueError:
                continue
            if score < coexpression_threshold:
                continue
            p1 = strip_taxon_prefix(parts[0])
            p2 = strip_taxon_prefix(parts[1])
            p1_in, p2_in = p1 in protein_ids, p2 in protein_ids
            if p1_in and not p2_in:
                add_partner(maps["protein_to_gene"].get(p2), "string_coexpression")
            if p2_in and not p1_in:
                add_partner(maps["protein_to_gene"].get(p1), "string_coexpression")

    # coxpresdb: one file per entrez id already in the set
    coxdir = input_dir / "coxpressdb"
    for entrez_id in entrez_ids:
        fp = coxdir / entrez_id
        if not fp.exists():
            continue
        with open(fp) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2 or parts[0] in entrez_ids:
                    continue
                add_partner(maps["entrez_to_gene"].get(parts[0]), "coxpresdb")

    # TFLink
    with open_maybe_gzip(resolve_source_file(input_dir, species, "tflink")) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tf_e = row.get("NCBI.GeneID.TF", "")
            tgt_e = row.get("NCBI.GeneID.Target", "")
            tf_u = row.get("UniprotID.TF", "")
            tgt_u = row.get("UniprotID.Target", "")
            tf_in = tf_e in entrez_ids or tf_u in uniprot_ids
            tgt_in = tgt_e in entrez_ids or tgt_u in uniprot_ids
            if tf_in and not tgt_in:
                add_partner(maps["entrez_to_gene"].get(tgt_e) or maps["uniprot_to_gene"].get(tgt_u), "tflink")
            if tgt_in and not tf_in:
                add_partner(maps["entrez_to_gene"].get(tf_e) or maps["uniprot_to_gene"].get(tf_u), "tflink")

    return partner_sources


def build_backbone_maps(input_dir, species, cfg):
    if cfg["is_flybase"]:
        return build_global_maps(input_dir, species, Path(cfg["gene_info_file"]))
    return build_generic_backbone_maps(input_dir, species)


def expand_gene_set(input_dir, species, anchor_gene_ids, maps, size_budget):
    id_sets = close_gene_set(anchor_gene_ids, maps)
    candidates = find_expansion_candidates(input_dir, species, id_sets, maps)
    ranked = sorted(candidates.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    budget_remaining = max(size_budget - len(anchor_gene_ids), 0)
    selected = [gid for gid, _ in ranked[:budget_remaining]]
    logger.info(
        "Expansion: %d candidate genes found via cross-references, selected %d "
        "(budget %d, anchors %d)",
        len(candidates), len(selected), size_budget, len(anchor_gene_ids),
    )
    final_gene_ids = set(anchor_gene_ids) | set(selected)
    return close_gene_set(final_gene_ids, maps)


# ---------------------------------------------------------------------------
# Phase 4: per-format filters
# ---------------------------------------------------------------------------

def filter_gtf(input_path, output_path, gene_ids):
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            m = GTF_GENE_ID_RE.search(line)
            if not m:
                continue
            gene_id = m.group(1).split(".")[0]
            if gene_id in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_uniprot_dat(input_path, output_path, uniprot_ids):
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        stanza = []
        for line in fin:
            stanza.append(line)
            if line.strip() == "//":
                accessions = set()
                for l in stanza:
                    if l.startswith("AC   "):
                        for acc in l[5:].strip().rstrip(";").split(";"):
                            accessions.add(acc.strip())
                if accessions & uniprot_ids:
                    fout.writelines(stanza)
                    kept += 1
                stanza = []
    return kept


def filter_gaf(input_path, output_path, id_set):
    """id_set is whichever ID space this species' GAF keys DB_Object_ID by —
    gene, uniprot accession, MGI, or RGD (see SPECIES[species]['gaf_id_space'])."""
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if line.startswith("!"):
                fout.write(line)
                continue
            parts = line.split("\t")
            if len(parts) > 1 and parts[1] in id_set:
                fout.write(line)
                kept += 1
    return kept


def filter_string_links(input_path, output_path, protein_ids, taxon_prefix):
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        header = next(fin, None)
        if header:
            fout.write(header)
        for line in fin:
            parts = line.split()
            if len(parts) < 2:
                continue
            p1 = parts[0][len(taxon_prefix):] if parts[0].startswith(taxon_prefix) else parts[0]
            p2 = parts[1][len(taxon_prefix):] if parts[1].startswith(taxon_prefix) else parts[1]
            if p1 in protein_ids and p2 in protein_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_coxpressdb(input_dir_path, output_dir_path, entrez_ids):
    output_dir_path.mkdir(parents=True, exist_ok=True)
    kept_files = kept_rows = 0
    for entrez_id in entrez_ids:
        src = input_dir_path / entrez_id
        if not src.exists():
            continue
        lines_kept = []
        with open(src) as fin:
            for line in fin:
                parts = line.strip().split("\t")
                if parts and parts[0] in entrez_ids:
                    lines_kept.append(line)
        if lines_kept:
            with open(output_dir_path / entrez_id, "w") as fout:
                fout.writelines(lines_kept)
            kept_files += 1
            kept_rows += len(lines_kept)
    return kept_files, kept_rows


def filter_tflink(input_path, output_path, uniprot_ids, entrez_ids):
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        header = next(reader)
        writer.writerow(header)
        idx = {name: i for i, name in enumerate(header)}

        def in_set(row, uni_col, entrez_col):
            return row[idx[uni_col]] in uniprot_ids or row[idx[entrez_col]] in entrez_ids

        for row in reader:
            tf_in = in_set(row, "UniprotID.TF", "NCBI.GeneID.TF")
            tgt_in = in_set(row, "UniprotID.Target", "NCBI.GeneID.Target")
            if tf_in and tgt_in:
                writer.writerow(row)
                kept += 1
    return kept


def filter_mitab_by_fbgn(input_path, output_path, gene_ids):
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            a = FBGN_RE.search(parts[0])
            b = FBGN_RE.search(parts[1])
            if a and b and a.group(0) in gene_ids and b.group(0) in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_gene_group(input_path, output_path, gene_ids, member_col=5, group_col=0):
    """Two-pass: keep every row for a group that has >=1 member in gene_ids,
    so the group node's full member list stays intact.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin:
        lines = fin.readlines()

    kept_groups = set()
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) > member_col and parts[member_col] in gene_ids:
            kept_groups.add(parts[group_col])

    kept = 0
    with open_maybe_gzip(output_path, "wt") as fout:
        for line in lines:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > group_col and parts[group_col] in kept_groups:
                fout.write(line)
                kept += 1
    return kept


def filter_by_fbal_intersection(input_path, output_path, allele_ids, fbal_col_regex=FBAL_RE, text_col=1):
    """Keep rows where any FBal id found (via regex) in text_col intersects allele_ids."""
    def keep_fn(parts):
        if len(parts) <= text_col:
            return False
        found = set(fbal_col_regex.findall(parts[text_col]))
        return bool(found & allele_ids)

    return copy_header_comments_and_filter(input_path, output_path, keep_fn, min_cols=text_col + 1)


def filter_gene_genetic_interactions(input_path, output_path, gene_ids):
    def keep_fn(parts):
        if len(parts) < 4:
            return False
        source_ids = set(FBGN_RE.findall(parts[1]))
        target_ids = set(FBGN_RE.findall(parts[3]))
        if not source_ids or not target_ids:
            return False
        return source_ids <= gene_ids and target_ids <= gene_ids

    return copy_header_comments_and_filter(input_path, output_path, keep_fn, min_cols=4)


def filter_afca(input_path, output_path, symbols_to_fbgn_path, gene_ids):
    """afca_afca_annotation_group_by_mean.tsv.gz: a wide matrix (one row per
    gene, hundreds of cell-type columns, ~6.2GB), keyed by gene SYMBOL in
    col0 (not FBgn id) -- e.g. "128up", not "FBgn...". ExpressionValueAdapter
    resolves symbol -> FBgn via fbgn_fbtr_fbpp_expanded.tsv.gz (gene_symbol
    col3 -> gene_ID col2) before matching; replicate that exact resolution
    here so kept rows are precisely the ones the real adapter would keep.
    Only the first tab-delimited field is ever parsed for rows we don't
    keep -- lines are huge, splitting the whole row would be wasteful.
    """
    symbol_to_fbgn = {}
    with open_maybe_gzip(symbols_to_fbgn_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 3:
                symbol_to_fbgn[parts[3]] = parts[2]

    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            tab_idx = line.find("\t")
            symbol = line[:tab_idx] if tab_idx != -1 else line.rstrip("\n")
            fbgn = symbol_to_fbgn.get(symbol)
            if fbgn is not None and fbgn in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_by_column_in_set(input_path, output_path, col, id_set, min_cols=None):
    min_cols = min_cols if min_cols is not None else col + 1

    def keep_fn(parts):
        return parts[col] in id_set

    return copy_header_comments_and_filter(input_path, output_path, keep_fn, min_cols=min_cols)


def filter_by_two_columns_in_set(input_path, output_path, col_a, col_b, id_set):
    def keep_fn(parts):
        return parts[col_a] in id_set and parts[col_b] in id_set

    return copy_header_comments_and_filter(
        input_path, output_path, keep_fn, min_cols=max(col_a, col_b) + 1
    )


def filter_by_any_column_in_set(input_path, output_path, cols, id_set, split_on=None, min_cols=None):
    """Keep rows where ANY of `cols` (optionally split_on-delimited multi-value
    fields, e.g. GWAS's semicolon-separated SNP_GENE_IDS) intersects id_set.
    """
    min_cols = min_cols if min_cols is not None else max(cols) + 1

    def keep_fn(parts):
        for col in cols:
            if col >= len(parts) or not parts[col]:
                continue
            values = parts[col].split(split_on) if split_on else [parts[col]]
            if any(v.strip() in id_set for v in values):
                return True
        return False

    return copy_header_comments_and_filter(input_path, output_path, keep_fn, min_cols=min_cols)


def filter_by_column_in_set_with_header(input_path, output_path, col, id_set, min_cols=None, strip_prefix=None):
    """Like filter_by_column_in_set, but the first line is a real column-name
    header (not a "#"-prefixed comment) that some adapters validate strictly
    (e.g. HPOAdapter/HPOGeneDiseaseAdapter) — write it through unfiltered,
    then filter every row after it. strip_prefix optionally strips a CURIE
    prefix (e.g. "NCBIGene:") from the compared value first, matching what
    the real adapter does internally before comparing IDs.
    """
    min_cols = min_cols if min_cols is not None else col + 1
    output_path.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < min_cols:
                continue
            value = parts[col]
            if strip_prefix and value.startswith(strip_prefix):
                value = value[len(strip_prefix):]
            if value in id_set:
                fout.write(line)
                kept += 1
    return kept


# ---------------------------------------------------------------------------
# Phase 2 (round 2): bgee / alliance / epd / reactome / hsa-only extras
# ---------------------------------------------------------------------------

def filter_bgee(input_path, output_path, gene_ids):
    """bgee's gene column is the bare gene ID (no CURIE prefix) — the adapter
    adds the taxon CURIE prefix itself at KG-build time."""
    return filter_by_column_in_set(input_path, output_path, col=0, id_set=gene_ids, min_cols=9)


def filter_alliance_disease(input_path, output_path, taxon_id, gene_ids):
    """DISEASE-ALLIANCE_COMBINED.tsv.gz is a single file combining every
    species Alliance tracks; every per-species alliance_gene_disease_*
    adapter block re-reads this same file at KG-build time, filtering by
    (taxon, label) itself — so here we only filter by taxon + gene closure,
    keeping all association-type rows for our species' genes. Reuses
    AllianceGeneDiseaseAdapter._resolve_bare_gene_id so the species-specific
    ID resolution (HGNC for hsa, BioMart for mmu/rno, strip-prefix for
    dmel/cel) exactly matches what the real adapter does later.
    """
    from biocypher_metta.adapters.alliance_gene_disease_adapter import AllianceGeneDiseaseAdapter, COLUMNS

    resolver = AllianceGeneDiseaseAdapter(
        filepath=str(input_path), label="is_implicated_in", taxon_id=taxon_id,
        write_properties=False, add_provenance=False,
    )
    taxon_str = str(taxon_id)
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= COLUMNS["db_object_id"] or parts[COLUMNS["taxon"]] == "Taxon":
                fout.write(line)  # header row
                continue
            if parts[COLUMNS["taxon"]].replace("NCBITaxon:", "") != taxon_str:
                continue
            if parts[COLUMNS["db_object_type"]] != "gene":
                continue
            resolved = resolver._resolve_bare_gene_id(parts[COLUMNS["db_object_id"]])
            if resolved is not None and resolved in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_alliance_orthology(input_path, output_path, taxon_id, gene_ids):
    """ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz: keep rows where Gene1 (our
    species' side) resolves into our gene closure. Gene2 (the ortholog, in
    another species) will always be a dangling edge in a single-species
    sample/KG — same accepted pattern as dmel's existing orthology_association
    (Phase 1) — since the target gene node only exists in a combined
    multi-species KG (--species all), not a single-species one.
    """
    from biocypher_metta.adapters.alliance_gene_orthology_adapter import AllianceGeneOrthologyAdapter, COLUMNS

    resolver = AllianceGeneOrthologyAdapter(
        filepath=str(input_path), label="orthologs_genes", taxon_id=taxon_id,
        write_properties=False, add_provenance=False,
    )
    target_taxon = f"NCBITaxon:{taxon_id}"
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= COLUMNS["gene2_species_name"] or parts[COLUMNS["gene1_id"]] == "Gene1ID":
                fout.write(line)  # header row
                continue
            if parts[COLUMNS["gene1_taxon"]] != target_taxon:
                continue
            resolved = resolver._resolve_bare_gene_id(parts[COLUMNS["gene1_id"]], str(taxon_id))
            if resolved is not None and resolved in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_epd(input_path, output_path, taxon_id, gene_ids, hgnc_to_ensembl_map=None):
    """EPD bed files are space-delimited (not tab), and the gene reference is
    a symbol embedded in the promoter name column (e.g. "NOC2L_1" -> "NOC2L"),
    not a direct ID — reuses EPDAdapter._resolve_symbol for the exact same
    per-species resolution the real adapter does (HGNC for hsa, legacy pickle
    for dmel, BioMart-backed GeneSymbolEnsemblProcessor for mmu/rno/cel).
    """
    from biocypher_metta.adapters.epd_adapter import EPDAdapter

    resolver = EPDAdapter(
        filepath=str(input_path), label="promoter", type="promoter", taxon_id=taxon_id,
        hgnc_to_ensembl_map=hgnc_to_ensembl_map, write_properties=False, add_provenance=False,
    )
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line)
                continue
            parts = line.split()
            if len(parts) <= EPDAdapter.INDEX["gene_id"]:
                continue
            symbol = parts[EPDAdapter.INDEX["gene_id"]].rsplit("_", 1)[0]
            ensembl_id, _ = resolver._resolve_symbol(symbol)
            if ensembl_id is not None and ensembl_id in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_reactome_gene_pathway_or_reaction(input_path, output_path, species_full_name, gene_ids,
                                              transcript_ids, uniprot_ids, ensembl_uniprot_map):
    """Shared filter for Ensembl2Reactome_All_Levels.txt (gene/protein/transcript
    -> pathway) and Ensembl2ReactomeReactions.txt (-> reaction): same 6-column
    layout (entity_id, pathway_or_reaction_id, url, name, evidence, species_name).
    Keeps a row if entity_id is a bare gene/transcript ID already in our
    closure, or (for protein-type IDs) resolves via the same
    EnsemblUniProtProcessor.mapping the real ReactomeEdgesAdapter uses to a
    UniProt ID in our closure. Returns (kept_row_count, referenced_ids set)
    so the caller can filter ReactomePathways.txt/ReactomePathwaysRelation.txt
    down to just the pathways actually reachable from our genes.
    """
    referenced_ids = set()
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if not line.strip():
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[5] != species_full_name:
                continue
            entity_id = parts[0].split(".")[0]  # strip Ensembl version suffix if present
            keep = entity_id in gene_ids or entity_id in transcript_ids
            if not keep and ensembl_uniprot_map:
                uniprot_id = ensembl_uniprot_map.get(parts[0]) or ensembl_uniprot_map.get(entity_id)
                keep = uniprot_id is not None and uniprot_id in uniprot_ids
            if keep:
                fout.write(line)
                kept += 1
                referenced_ids.add(parts[1])
    return kept, referenced_ids


def filter_reactome_pathways_relation(input_path, output_path, reactome_prefix, referenced_ids):
    """ReactomePathwaysRelation.txt: parent_id, child_id. Keeps species-prefixed
    rows touching a pathway already referenced by our filtered gene/pathway
    edges, and folds in the opposite endpoint too (one hop of hierarchy
    context). Returns the full retained pathway-id set (referenced_ids plus
    anything pulled in by this pass) for filtering ReactomePathways.txt.
    """
    retained = set(referenced_ids)
    kept_lines = []
    with open_maybe_gzip(input_path) as fin:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0].startswith(reactome_prefix):
                continue
            if parts[0] in referenced_ids or parts[1] in referenced_ids:
                kept_lines.append(line)
                retained.add(parts[0])
                retained.add(parts[1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(output_path, "wt") as fout:
        for line in kept_lines:
            fout.write(line)
    return len(kept_lines), retained


def filter_reactome_pathways_nodes(input_path, output_path, species_full_name, retained_pathway_ids):
    """ReactomePathways.txt: pathway_id, name, species_name — global file
    across every organism Reactome tracks; keep species rows whose ID was
    actually referenced by our filtered edges/hierarchy."""
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if not line.strip():
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == species_full_name and parts[0] in retained_pathway_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_reactome_small_molecule(input_path, output_path, species_full_name, retained_ids):
    """ChEBI2Reactome_All_Levels.txt / ChEBI2ReactomeReactions.txt: same
    6-column layout as Ensembl2Reactome_All_Levels.txt/Ensembl2ReactomeReactions.txt
    (chebi_id, pathway_or_reaction_id, url, name, evidence, species_name) —
    keep rows whose pathway/reaction was already retained by the gene-based
    reactome filter (filter_reactome_gene_pathway_or_reaction), so every
    small_molecule edge lands on a pathway/reaction node that actually exists
    in the sample.
    """
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if not line.strip():
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6 and parts[5] == species_full_name and parts[1] in retained_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_reactome_ppi(input_path, output_path, reactome_prefix, gene_ids):
    """reactome.all_species.interactions.tab-delimited.txt: interactor 1/2
    Ensembl gene id columns are pipe-separated bundles of ENST/ENSP/ENSG ids
    (e.g. "ENSEMBL:ENST...|ENSEMBL:ENSP...|ENSEMBL:ENSG..."); keep a row if
    the 'reactome:R-<PREFIX>-...' context column matches our species AND
    either interactor's bundle contains a bare Ensembl gene id in gene_ids.
    """
    def has_our_gene(bundle_col):
        for token in bundle_col.split("|"):
            token = token.split(":")[-1]
            if token in gene_ids:
                return True
        return False

    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or not parts[7].startswith(f"reactome:{reactome_prefix}"):
                continue
            if has_our_gene(parts[1]) or has_our_gene(parts[4]):
                fout.write(line)
                kept += 1
    return kept


def filter_reactome_reaction_exporter(input_path, output_path, retained_pathway_ids, reaction_ids):
    """reactome_reaction_exporter_All_species.txt: pathway_id, reaction_id,
    reaction_name, uniprot_acc, role_in_reaction (no header). Feeds both
    reactome_reaction_to_pathway (ReactomeEdgesAdapter) and the 5
    protein-role adapters (ReactomeInferenceEdgesAdapter), which emit both a
    protein->pathway and a protein->reaction edge per row. Require BOTH ids
    already in our closure's retained sets so neither side ever dangles.
    """
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if parts[0] in retained_pathway_ids and parts[1] in reaction_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_reactome_go_terms(input_path, output_path, retained_ids):
    """Pathways2GoTerms_human.txt / Reactions2GoTerms_human.txt: header +
    'Identifier\tName\tGO_Term' rows, Identifier = pathway_id or reaction_id.
    Keep rows whose Identifier is already in our closure's retained set.
    """
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 1 and parts[0] in retained_ids:
                fout.write(line)
                kept += 1
    return kept


# ---------------------------------------------------------------------------
# Phase 3: TADmap / TFBS / rna_central / ABC+dbSNP
# ---------------------------------------------------------------------------

def filter_tadmap(input_path, output_path, taxon_id, gene_ids):
    """TADMap_geneset_hs.csv: col0 = "chr|start|end", col1 = ";"-separated
    "SYMBOL|Ensembl:ENSGxxx|HGNC:SYMBOL" entries — same parsing TADMapAdapter
    itself does. Keeps a TAD row if any of its genes is in our closure.
    """
    from biocypher_metta.adapters import Adapter

    prefix = Adapter.CURIE_PREFIX[taxon_id]
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 2:
                continue
            keep = False
            for gene_entry in parts[1].split(";"):
                fields = gene_entry.split("|")
                if len(fields) < 2:
                    continue
                ensembl_id = fields[1].split(":")[-1].upper()
                if f"{prefix}:{ensembl_id}" in gene_ids or ensembl_id in gene_ids:
                    keep = True
                    break
            if keep:
                fout.write(line)
                kept += 1
    return kept


def filter_tfbs(input_path, output_path, gene_ids, hgnc_processor):
    """encRegTfbsClustered.txt.gz: tab-delimited, no header, col4 = HGNC gene
    symbol — resolve via the same HGNCProcessor TfbsAdapter uses, keep the
    row if the resolved gene is in our closure."""
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(input_path) as fin, open_maybe_gzip(output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            ensembl_id = hgnc_processor.get_ensembl_id(parts[4].strip('"'))
            if ensembl_id is not None and ensembl_id in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def _load_gene_regions_from_filtered_gtf(gtf_path):
    """Parse chr/start/end for 'gene' feature rows out of the sample's own
    already-filtered (tiny, ~180-gene) GTF — reused so rna_central filtering
    doesn't need to re-parse the full-size original GTF."""
    regions = defaultdict(list)
    with open_maybe_gzip(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or parts[2] != "gene":
                continue
            chr_ = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            start, end = int(parts[3]), int(parts[4])
            regions[chr_].append((start, end))
    return regions


def filter_rna_central(bed_path, rfam_path, bed_output_path, rfam_output_path, gtf_output_path, taxon_id):
    """Keep RNAcentral BED rows whose genomic interval overlaps any gene
    region in the sample's own filtered GTF (plain linear scan — the gene
    list is tiny, ~180 entries, so no interval tree needed here), then keep
    rfam GO-annotation rows for exactly the RNA ids that survived.
    """
    gene_regions = _load_gene_regions_from_filtered_gtf(gtf_output_path if gtf_output_path.exists() else Path(gtf_output_path))

    def overlaps_any_gene(chr_, start, end):
        for g_start, g_end in gene_regions.get(chr_, ()):
            if start <= g_end and end >= g_start:
                return True
        return False

    kept_rna_ids = set()
    bed_kept = 0
    bed_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(bed_path) as fin, open_maybe_gzip(bed_output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chr_ = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            start = int(parts[1].strip()) + 1  # BED is 0-indexed
            end = int(parts[2].strip())
            if overlaps_any_gene(chr_, start, end):
                fout.write(line)
                bed_kept += 1
                kept_rna_ids.add(parts[3].split("_")[0])

    rfam_kept = 0
    rfam_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(rfam_path) as fin, open_maybe_gzip(rfam_output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            rna_id = parts[0]
            if not rna_id.endswith(f"_{taxon_id}"):
                continue
            if rna_id.split("_")[0] in kept_rna_ids:
                fout.write(line)
                rfam_kept += 1

    return bed_kept, rfam_kept


def filter_abc(abc_path, abc_output_path, gene_ids, hgnc_processor):
    """abc.forgedb.csv.gz: col10 = target_gene (HGNC symbol). Keep rows whose
    resolved gene is in our closure. Returns (kept_rows, triples) where
    triples are (rsid, chr, pos) for the shared dbSNP cache — rsid/
    chromosome/start_position are already columns 0/1/2, no VCF needed.
    """
    kept_rows = []
    with open_maybe_gzip(abc_path) as fin:
        header = next(fin, None)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 11:
                continue
            ensembl_id = hgnc_processor.get_ensembl_id(parts[10])
            if ensembl_id is not None and ensembl_id in gene_ids:
                kept_rows.append(parts)

    abc_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(abc_output_path, "wt") as fout:
        if header is not None:
            fout.write(header)
        for parts in kept_rows:
            fout.write(",".join(parts) + "\n")

    triples = [(p[0], p[1], p[2]) for p in kept_rows]
    return len(kept_rows), triples


def filter_cadd(cadd_path, cadd_output_path, chr_filter="chr16"):
    """cadd.forgedb.csv.gz: rsid=0, chromosome=1, position=2. CADD has no
    gene reference at all (get_edges() is a no-op in CADDAdapter — it only
    emits SNP nodes) and the real full config restricts it to chr_filter, so
    filtering by chromosome alone matches what the adapter will actually
    process. Returns (kept_rows, triples) for the shared dbSNP cache.
    """
    kept_rows = []
    with open_maybe_gzip(cadd_path) as fin:
        header = next(fin, None)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 7:
                continue
            if parts[1] == chr_filter:
                kept_rows.append(parts)

    cadd_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(cadd_output_path, "wt") as fout:
        if header is not None:
            fout.write(header)
        for parts in kept_rows:
            fout.write(",".join(parts) + "\n")

    triples = [(p[0], p[1], p[2]) for p in kept_rows]
    return len(kept_rows), triples


def filter_refseq(refseq_path, refseq_output_path, gene_ids, hgnc_processor):
    """closest_gene.forgedb.csv.gz: rsid=0, chromosome=1, start=2, end=3,
    gene_symbol=7 (HGNC symbol). The real full config has no chr/start/end
    restriction (whole genome, gene-driven), so filter by resolved gene
    membership only. Returns (kept_rows, triples) for the shared dbSNP cache
    (end_position == the SNP's own dbSNP position, matching CADD's `position`
    for the same rsid).
    """
    kept_rows = []
    with open_maybe_gzip(refseq_path) as fin:
        header = next(fin, None)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 8:
                continue
            ensembl_id = hgnc_processor.get_ensembl_id(parts[7])
            if ensembl_id is not None and ensembl_id in gene_ids:
                kept_rows.append(parts)

    refseq_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(refseq_output_path, "wt") as fout:
        if header is not None:
            fout.write(header)
        for parts in kept_rows:
            fout.write(",".join(parts) + "\n")

    triples = [(p[0], p[1], p[3]) for p in kept_rows]
    return len(kept_rows), triples


def filter_topld(topld_path, topld_output_path, start=53000000, end=56000000):
    """<POP>_chr16_no_filter_0.2_1000000_LD.csv.gz: SNP1=0, SNP2=1 (integer
    chr16 positions, no rsid column at all — TopLDAdapter resolves rsids via
    a reverse dbsnp_pos_map lookup). Keep rows where either endpoint falls in
    [start, end], matching the real full config's chr16:53M-56M window.
    """
    kept = 0
    topld_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(topld_path) as fin, open_maybe_gzip(topld_output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 2:
                continue
            try:
                pos1, pos2 = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            if (start <= pos1 <= end) or (start <= pos2 <= end):
                fout.write(line)
                kept += 1
    return kept


def build_dbsnp_cache(triples, dbsnp_cache_dir):
    """Write a scoped rsid_to_pos SQLite cache (DBSNPProcessor's schema) from
    (rsid, chr, pos) triples gathered from ABC/CADD/RefSeq's own rows — no
    genome-wide VCF processing. The same table serves both the forward
    dbsnp_rsid_map (ABC/CADD/RefSeq) and reverse dbsnp_pos_map (TopLD) lookups,
    since DBSNPProcessor.get_dict_wrappers() queries this one table both ways.
    """
    import sqlite3
    import json as _json

    dbsnp_cache_dir = Path(dbsnp_cache_dir)
    dbsnp_cache_dir.mkdir(parents=True, exist_ok=True)
    db_path = dbsnp_cache_dir / "dbsnp_mapping.db"
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE rsid_to_pos (rsid TEXT PRIMARY KEY, chr TEXT, pos INTEGER)")
    conn.execute("CREATE INDEX idx_pos ON rsid_to_pos (chr, pos)")
    seen_rsids = set()
    rows_to_insert = []
    for rsid, chrom, pos in triples:
        if rsid in seen_rsids or not rsid.startswith("rs"):
            continue
        seen_rsids.add(rsid)
        rows_to_insert.append((rsid, chrom, int(pos)))
    conn.executemany("INSERT INTO rsid_to_pos VALUES (?, ?, ?)", rows_to_insert)
    conn.commit()
    conn.close()

    with open(dbsnp_cache_dir / "dbsnp_version.json", "w") as f:
        _json.dump({
            "sample_from_common": True,
            "filter": "scoped to connected-sample gene/region closure (built from abc/cadd/refseq rows directly, no VCF processed)",
            "entries": len(rows_to_insert),
        }, f, indent=2)

    return len(rows_to_insert)


def filter_gtex(gtex_path, gtex_output_path, gene_ids):
    """gtex.forgedb.csv.gz: comma-delimited (despite living next to mostly
    tab-delimited siblings) — col2 = ensembl_gene_id, already bare/unversioned,
    no HGNC resolution or dbsnp_rsid_map needed (chr/pos are inline columns
    18/19). Shared by both gtex_eqtl and gtex_expression adapters."""
    kept = 0
    gtex_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(gtex_path) as fin, open_maybe_gzip(gtex_output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) >= 22 and parts[2] in gene_ids:
                fout.write(line)
                kept += 1
    return kept


def filter_hocomoco(annotation_path, annotation_output_path, pwm_input_dir, pwm_output_dir, gene_ids, hgnc_processor):
    """HOCOMOCOv11_core_annotation_HUMAN_mono.tsv (401 models, ~1.8MB total incl.
    pwm/) is copied wholesale, unfiltered. HOCOMOCOAdapter emits the node id as
    the TF's own ENSEMBL:<id>, and motif_diff's edges reference the TF's Model
    column across its ~400-wide header regardless of our gene closure -- if we
    filtered HOCOMOCO down to closure-only TFs (as an earlier version of this
    function did), motif_diff's SNP->motif edges would mostly dangle (target
    'motif' node never created). The data is small enough that filtering buys
    nothing; copy every model/pwm file so motif_diff's full TF universe always
    resolves.
    """
    annotation_output_path.parent.mkdir(parents=True, exist_ok=True)
    annotation_output_path.write_bytes(annotation_path.read_bytes())

    pwm_output_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    for src_pwm in pwm_input_dir.glob("*.pwm"):
        (pwm_output_dir / src_pwm.name).write_bytes(src_pwm.read_bytes())
        copied += 1

    with open_maybe_gzip(annotation_output_path) as fin:
        next(fin, None)
        n_models = sum(1 for _ in fin)
    return n_models, copied


def filter_dbsuper(dbsuper_path, dbsuper_output_path, gene_ids, hgnc_processor):
    """dbSUPER_SuperEnhancers_hg19.tsv.gz: col4 = gene_symbol (HGNC symbol);
    hg19->hg38 liftover is handled inside the adapter itself, not here."""
    def keep_fn(parts):
        if len(parts) < 6:
            return False
        ensembl_id = hgnc_processor.get_ensembl_id(parts[4])
        return ensembl_id is not None and ensembl_id in gene_ids

    kept = 0
    dbsuper_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(dbsuper_path) as fin, open_maybe_gzip(dbsuper_output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if keep_fn(parts):
                fout.write(line)
                kept += 1
    return kept


def filter_peregrine(enhancers_path, sources_path, gene_link_path,
                      enhancers_output_path, sources_output_path, gene_link_output_path,
                      gene_ids, hgnc_processor):
    """enhancer_gene_link_18.tsv.gz: col1 = "HUMAN|HGNC=<id>|..." (99.994% of
    rows) — same HGNC-ID resolution PeregrineAdapter.handle_gene() does. Keep
    rows whose resolved gene is in our closure, collect the referenced
    enhancer ids (col0), then filter PEREGRINEenhancershg38.gz (col3=enhancer_id)
    and PEREGRINEenhancersources.gz (col0=enhancer_id) down to just those ids.
    """
    referenced_enhancer_ids = set()
    kept_link_rows = 0
    gene_link_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(gene_link_path) as fin, open_maybe_gzip(gene_link_output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene_field = parts[1]
            if not gene_field.startswith("HUMAN|HGNC="):
                continue
            hgnc_id = "HGNC:" + gene_field.split("HGNC=", 1)[1].split("|")[0]
            ensembl_id = hgnc_processor.get_ensembl_id(hgnc_id)
            if ensembl_id is not None and ensembl_id in gene_ids:
                fout.write(line)
                kept_link_rows += 1
                referenced_enhancer_ids.add(parts[0])

    kept_enhancers = 0
    enhancers_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(enhancers_path) as fin, open_maybe_gzip(enhancers_output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[3] in referenced_enhancer_ids:
                fout.write(line)
                kept_enhancers += 1

    kept_sources = 0
    sources_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(sources_path) as fin, open_maybe_gzip(sources_output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 1 and parts[0] in referenced_enhancer_ids:
                fout.write(line)
                kept_sources += 1

    return kept_link_rows, kept_enhancers, kept_sources


def collect_column_values(path, col, min_cols=None, split_on=None, delimiter="\t", skip_header=False):
    """Read an already-filtered sample file and collect the distinct values
    in one column — used to gather the exact rsid set referenced by GWAS/GTEx
    sample rows, to know which dbsnp_snps nodes are actually needed."""
    min_cols = min_cols if min_cols is not None else col + 1
    values = set()
    with open_maybe_gzip(path) as fin:
        if skip_header:
            next(fin, None)
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(delimiter)
            if len(parts) < min_cols or not parts[col]:
                continue
            if split_on:
                values.update(v.strip() for v in parts[col].split(split_on) if v.strip())
            else:
                values.add(parts[col])
    return values


def filter_enhancer_atlas(hs_bed_path, hs_bed_output_path, enhancer_gene_input_dir, enhancer_gene_output_dir, gene_ids):
    """enhancer_atlas has no gene column in its node file (hs.bed.gz) — the
    gene link only exists in the per-tissue enhancer_gene/<tissue>_EP.txt
    files, each line shaped "chr:start-end_ENSGxxxxx$SYMBOL$chr$tss$strand\\tscore".
    Filter every tissue file by gene membership, collect the referenced
    (chr, start+1, end) regions (matching EnhancerAtlasAdapter's own +1
    0-based-to-1-based offset), then filter hs.bed.gz down to just those
    regions so every enhancer node kept has a real edge.
    """
    referenced_regions = set()
    enhancer_gene_output_dir.mkdir(parents=True, exist_ok=True)
    kept_edge_rows = 0
    for tissue_file in sorted(enhancer_gene_input_dir.iterdir()):
        if not tissue_file.is_file():
            continue
        out_path = enhancer_gene_output_dir / tissue_file.name
        with open(tissue_file, "r") as fin, open(out_path, "w") as fout:
            for line in fin:
                info = line.strip().split("\t")
                if not info or "_" not in info[0]:
                    continue
                enhancer_info, gene_part = info[0].split("_", 1)
                gene = gene_part.split("$")[0]
                if gene not in gene_ids:
                    continue
                chrom, coords = enhancer_info.split(":")
                start, end = coords.split("-")
                referenced_regions.add((chrom, int(start) + 1, int(end)))
                fout.write(line)
                kept_edge_rows += 1

    kept_nodes = 0
    hs_bed_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(hs_bed_path) as fin, open_maybe_gzip(hs_bed_output_path, "wt") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start, end = int(parts[1]) + 1, int(parts[2])
            except ValueError:
                continue
            if (parts[0], start, end) in referenced_regions:
                fout.write(line)
                kept_nodes += 1
    return kept_nodes, kept_edge_rows


def filter_ccre_closest_genes(input_path, output_path, gene_ids):
    """GRCh38-Closest-Genes-{All,PC}.tsv.gz: no header, col12 = versioned
    Ensembl gene id (e.g. "ENSG00000186092.7") — strip the version before
    comparing. Shared by every enhancer/promoter cCRE block that reads this
    file (element_type/edge_type just select which rows those adapters use
    at KG-build time), so filter once per unique file and reuse the output.
    """
    def keep_fn(parts):
        return parts[12].split(".")[0] in gene_ids

    return copy_header_comments_and_filter(input_path, output_path, keep_fn, min_cols=13)


def filter_ccre_eqtl(input_path, output_path, gene_ids):
    """V4-hg38.Gene-Links.eQTLs.txt.gz: no header, col1 = bare (unversioned)
    Ensembl gene id already."""
    return filter_by_column_in_set(input_path, output_path, col=1, id_set=gene_ids, min_cols=9)


def filter_catlas(abc_scores_input_dir, abc_scores_output_dir, ccre_master_path, ccre_master_output_path,
                   ccres_input_dir, ccres_output_dir, ccre_label_pkl_path, gene_ids, hgnc_processor):
    """Catlas is 3 linked pieces: per-cell-type ABC_scores/*.tsv.gz (gene
    resolution: "Gene Name" col = "SYMBOL:transcript_id"), the cCRE_hg38.tsv.gz
    master catalog (coordinates + Class, no gene column), and per-cell-type
    cCREs/*.bed (coordinates only). Filter ABC_scores by gene first, collect
    the referenced cCRE coordinates, then filter the master catalog and every
    per-cell-type .bed down to just those coordinates so everything stays
    consistent. Cell_ontology.tsv / catlas_abc_cell_type_aliases.tsv are tiny
    filename-keyed lookups (not coordinate-dependent) — copy them unfiltered,
    reused as-is by the caller.

    Also builds catlas_ccre_label_map.pkl directly (mirrors
    scripts/create_catlas_ccre_label_map.py's exact logic) at a
    connected-sample-specific path — this pkl is coordinate-dependent, so it
    must NOT be built at the shared aux_files/hsa/catlas/ location the real
    full/legacy-sample runs use, or it would silently corrupt their cache with
    our small subset.
    """
    referenced_ccre_keys = set()
    abc_scores_output_dir.mkdir(parents=True, exist_ok=True)
    kept_abc_files = 0
    kept_abc_rows = 0
    for f in sorted(abc_scores_input_dir.iterdir()):
        if not (f.name.endswith(".tsv") or f.name.endswith(".tsv.gz")):
            continue
        out_path = abc_scores_output_dir / f.name
        kept_here = 0
        with open_maybe_gzip(f) as fin, open_maybe_gzip(out_path, "wt") as fout:
            header = next(fin, None)
            if header is not None:
                fout.write(header)
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                symbol = parts[3].split(":")[0]
                ensembl_id = hgnc_processor.get_ensembl_id(symbol)
                if ensembl_id is None or ensembl_id not in gene_ids:
                    continue
                fout.write(line)
                kept_here += 1
                for coord_field in (parts[0], parts[1]):  # cCRE, Promoter columns
                    if ":" not in coord_field or "-" not in coord_field:
                        continue
                    chrom, rng = coord_field.split(":", 1)
                    start, end = rng.split("-", 1)
                    try:
                        referenced_ccre_keys.add((chrom, int(start), int(end)))
                    except ValueError:
                        continue
        if kept_here > 0:
            kept_abc_files += 1
            kept_abc_rows += kept_here
        else:
            out_path.unlink(missing_ok=True)

    ccre_label_map = {}
    ccre_master_output_path.parent.mkdir(parents=True, exist_ok=True)
    kept_master = 0
    with open_maybe_gzip(ccre_master_path) as fin, open_maybe_gzip(ccre_master_output_path, "wt") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start0, end = parts[0], int(parts[1]), int(parts[2])
            if (chrom, start0, end) not in referenced_ccre_keys:
                continue
            fout.write(line)
            kept_master += 1
            cls = parts[3].strip().lower()
            label = "enhancer" if cls == "distal" else ("promoter" if cls in ("promoter", "promoter proximal") else None)
            if label is not None:
                ccre_label_map[(chrom, start0 + 1, end)] = label

    ccres_output_dir.mkdir(parents=True, exist_ok=True)
    kept_bed_files = 0
    for f in sorted(ccres_input_dir.iterdir()):
        if not f.is_file():
            continue
        out_path = ccres_output_dir / f.name
        kept_here = 0
        with open(f) as fin, open(out_path, "w") as fout:
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    key = (parts[0], int(parts[1]), int(parts[2]))
                except ValueError:
                    continue
                if key in referenced_ccre_keys:
                    fout.write(line)
                    kept_here += 1
        if kept_here > 0:
            kept_bed_files += 1
        else:
            out_path.unlink()

    import pickle
    ccre_label_pkl_path.parent.mkdir(parents=True, exist_ok=True)
    with open(ccre_label_pkl_path, "wb") as f:
        pickle.dump(ccre_label_map, f)

    return kept_abc_files, kept_abc_rows, kept_master, kept_bed_files


def filter_motif_diff(motif_diff_path, motif_diff_output_path, rsid_set, max_rows=None):
    """_mono_probNorm_average.diff: ~10GB after filtering to real_rows alone
    (770 TF-model score columns/row is the actual size driver, not row
    count), tab-delimited, col0=bare rsid (no CURIE prefix —
    MotifDiffAdapter.get_edges() prefixes with DBSNP: itself). MotifDiffAdapter
    has no dbsnp_rsid_map/region filter at all (confirmed — it's the only
    rsid-keyed adapter without one), so the only workable filter is rsid
    membership against an externally supplied set (GWAS + GTEx-eqtl + the
    ABC/CADD/RefSeq dbSNP cache rsids). Confirmed (2026-08-01, via a
    connectivity check with this edge type entirely excluded) that
    alters_binding contributes zero essential gene-closure connectivity —
    it's pure schema/coverage breadth, so max_rows caps it to a small
    representative sample rather than every matching row; the scan stops
    as soon as the cap is hit, so this is also usually much faster than a
    full linear pass once max_rows is reached early.
    """
    kept = 0
    motif_diff_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(motif_diff_path, "r") as fin, open(motif_diff_output_path, "w") as fout:
        header = next(fin, None)
        if header is not None:
            fout.write(header)
        for line in fin:
            if max_rows is not None and kept >= max_rows:
                break
            tab_idx = line.find("\t")
            rsid = line[:tab_idx] if tab_idx != -1 else line.rstrip("\n")
            if rsid in rsid_set:
                fout.write(line)
                kept += 1
    return kept


def filter_roadmap_dir(input_dir, output_dir, rsid_set):
    """Roadmap Epigenomics forgedb CSVs (chromatin_state/h3_marks, each split
    across 10 parts, ~480M+ rows each): header
    'rsid,dataset,cell,tissue,datatype'. Comment previously said this needs
    csv.reader throughout because cell/tissue/datatype can carry quoted
    commas -- verified against the real files (all 10 chromatin_state parts,
    482M+ rows each, zero quote characters anywhere) that this never
    actually happens. The overwhelming majority of rows are discarded (rsid
    not in the closure's window), so csv.reader's per-row parsing overhead
    across every field was pure waste for them. Only the rsid (first field,
    never contains a comma) is needed to decide membership -- read raw
    lines and cheaply slice out just that field; only rows that actually
    pass the filter get properly parsed (via csv.reader, still, for
    correctness/safety) before being rewritten. Every adapter
    (RoadMapChromatinStateAdapter/RoadMapH3MarkAdapter) resolves rsid ->
    chr/pos via dbsnp_rsid_map at runtime and then applies the config's
    chr16:53M-56M window -- so filtering rows here to rsids already known
    (via that same cache) to fall in that window reproduces exactly what
    the adapter would keep from the full files, without downloading/keeping
    the full ~20GB.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    total_kept = 0
    for part_path in sorted(input_dir.glob("*.csv.gz")):
        kept = 0
        out_path = output_dir / part_path.name
        with gzip.open(part_path, "rt", newline="") as fin, gzip.open(out_path, "wt", newline="") as fout:
            header = next(fin, None)
            writer = None
            if header is not None:
                writer = csv.writer(fout, delimiter=",")
                writer.writerow(next(csv.reader([header])))
            for line in fin:
                comma_idx = line.find(",")
                rsid = line[:comma_idx] if comma_idx != -1 else line.rstrip("\n")
                if rsid not in rsid_set:
                    continue
                row = next(csv.reader([line]))
                if writer is None:
                    writer = csv.writer(fout, delimiter=",")
                writer.writerow(row)
                kept += 1
        total_kept += kept
    return total_kept


def filter_roadmap_file(input_path, output_path, rsid_set):
    """forge2.erc2-DHS.forgedb.csv.gz: single file (not a directory), same
    'rsid,dataset,cell,tissue,datatype' format as chromatin_state/h3_marks.
    Same rsid-first fast path as filter_roadmap_dir -- see its docstring.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with gzip.open(input_path, "rt", newline="") as fin, gzip.open(output_path, "wt", newline="") as fout:
        header = next(fin, None)
        writer = None
        if header is not None:
            writer = csv.writer(fout, delimiter=",")
            writer.writerow(next(csv.reader([header])))
        for line in fin:
            comma_idx = line.find(",")
            rsid = line[:comma_idx] if comma_idx != -1 else line.rstrip("\n")
            if rsid in rsid_set:
                row = next(csv.reader([line]))
                if writer is None:
                    writer = csv.writer(fout, delimiter=",")
                writer.writerow(row)
                kept += 1
    return kept


def filter_dbsnp_snps(vcf_path, vcf_output_path, rsid_set):
    """00-common_all.vcf.gz: col2 = rsID (e.g. "rs367896724"). Keep only the
    rsids already referenced by our filtered GWAS/GTEx-eqtl rows, so
    dbsnp_snps produces exactly the SNP nodes those edges need (closing their
    dangling DBSNP: targets) without touching the 1.6GB file's full contents.
    """
    kept = 0
    vcf_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(vcf_path) as fin, open_maybe_gzip(vcf_output_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] in rsid_set:
                fout.write(line)
                kept += 1
    return kept


def filter_eva_snps(vcf_path, vcf_output_path, gtf_output_path):
    """EVA VCF (dmel only for now): keep variant rows whose position falls
    inside a gene body from the sample's own already-filtered GTF — same
    position-overlap pattern as filter_rna_central, but simpler (chrom+pos
    only, no interval). dmel's chromosome names (2L/2R/3L/3R/4/X/MT) already
    match between the EVA VCF and the GTF with no 'chr' prefix, so this
    doesn't need _load_gene_regions_from_filtered_gtf's gencode-style
    normalization.
    """
    gene_regions = defaultdict(list)
    with open_maybe_gzip(gtf_output_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or parts[2] != "gene":
                continue
            gene_regions[parts[0]].append((int(parts[3]), int(parts[4])))

    def overlaps_any_gene(chrom, pos):
        return any(start <= pos <= end for start, end in gene_regions.get(chrom, ()))

    kept = 0
    vcf_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(vcf_path) as fin, open_maybe_gzip(vcf_output_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or not parts[1].isdigit():
                continue
            if overlaps_any_gene(parts[0], int(parts[1])):
                fout.write(line)
                kept += 1
    return kept


# ---------------------------------------------------------------------------
# Phase 5: synthetic fallback
# ---------------------------------------------------------------------------

def write_synthetic_tsv(output_path, num_columns, id_columns_values, num_rows=3):
    """Write a minimal synthetic .SYNTHETIC. file: num_rows rows of num_columns
    columns, with the given {col_index: [values...]} cycled in, "SYNTHETIC"
    elsewhere. Tagged filename, never overwrites the real (empty) file.
    """
    suffix = "".join(output_path.suffixes)
    stem = output_path.name[: -len(suffix)] if suffix else output_path.name
    synthetic_path = output_path.with_name(f"{stem}.SYNTHETIC{suffix}")
    synthetic_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(synthetic_path, "wt") as fout:
        for r in range(num_rows):
            row = ["SYNTHETIC"] * num_columns
            for col, values in id_columns_values.items():
                if values:
                    row[col] = values[r % len(values)]
            fout.write("\t".join(row) + "\n")
    return synthetic_path, num_rows


def write_synthetic_uniprot_chebi_part_of(output_path, part_chebi_id, ligand_chebi_id, taxon_id):
    """UniprotProteinAdapter's 'chemical_substance_part_of_chemical_substance'
    label reads BINDING features' /ligand_id and /ligand_part_id qualifiers
    from uniprot_sprot_human.dat.gz -- a real but rare annotation (a binding
    site contacting only part of a larger ligand). Our closure's 179-protein
    filtered file legitimately has zero such qualifiers (confirmed: zero
    "DR   ChEBI" lines at all), so there is no real row to keep, only a
    minimal record to fabricate.

    Written as a STANDALONE synthetic .dat.gz (never appended to the real,
    shared uniprot_sprot_human.dat.gz — every other uniprot_* adapter reads
    that same file and would otherwise pick up a fake protein). This
    adapter's own config block gets its filepath repointed at this file
    instead; since this label only ever extracts BINDING-feature ligand
    qualifiers (see uniprot_protein_adapter.py's dbxref == "CHEBI" branch)
    and (nodes: False in this block) never emits a node, nothing else about
    the run is affected. part_chebi_id/ligand_chebi_id must be real small
    molecule CHEBI ids already present in this build's ChEBI ontology.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def ft_line(name, location):
        return f"FT   {name:<8}{'':<8}{location}\n"

    def ft_qual(key, value):
        return f"FT{'':19}/{key}=\"{value}\"\n"

    record_text = (
        "ID   SYNTH_ORG               Reviewed;          10 AA.\n"
        "AC   Q00SYN0;\n"
        f"OX   NCBI_TaxID={taxon_id};\n"
        + ft_line("BINDING", "1")
        + ft_qual("ligand", "synthetic placeholder ligand")
        + ft_qual("ligand_id", f"ChEBI:CHEBI:{ligand_chebi_id}")
        + ft_qual("ligand_part", "synthetic placeholder ligand part")
        + ft_qual("ligand_part_id", f"ChEBI:CHEBI:{part_chebi_id}")
        + "//\n"
    )
    with open_maybe_gzip(output_path, "wt") as fout:
        fout.write(record_text)
    return output_path, 1


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def _generate_core(species, input_dir, output_dir, cfg, id_sets, record, src, dst):
    """The adapter set every species has: gencode/uniprot/GAF/STRING/coxpresdb/TFLink."""
    n = filter_gtf(src("gtf"), dst("gtf"), id_sets["gene"])
    record("gencode_gtf", n)

    n = filter_uniprot_dat(src("uniprot_dat"), dst("uniprot_dat"), id_sets["uniprot"])
    record("uniprot_sprot", n)

    gaf_id_space = cfg.get("gaf_id_space", "gene")
    gaf_id_set = {
        "gene": id_sets["gene"],
        "uniprot": id_sets["uniprot"],
    }.get(gaf_id_space, id_sets.get("gafid", set()))
    n = filter_gaf(src("gaf"), dst("gaf"), gaf_id_set)
    record("gaf", n)

    taxon_prefix = cfg["string_taxon_prefix"]
    n = filter_string_links(src("string_ppi"), dst("string_ppi"), id_sets["protein"], taxon_prefix)
    record("string_ppi", n)

    n = filter_string_links(
        src("string_coexpression"), dst("string_coexpression"), id_sets["protein"], taxon_prefix,
    )
    record("string_coexpression", n)

    kept_files, kept_rows = filter_coxpressdb(
        input_dir / "coxpressdb", output_dir / "coxpressdb", id_sets["entrez"]
    )
    record("coxpresdb", kept_rows, empty_ok_reason="no coxpresdb coverage for genes in set (e.g. histones)")

    n = filter_tflink(src("tflink"), dst("tflink"), id_sets["uniprot"], id_sets["entrez"])
    record("tflink", n)


def _generate_dmel_extension(species, input_dir, output_dir, id_sets, record, src, dst):
    """FlyBase-specific adapters that only exist for dmel: gene_group,
    disease_model, genotype_phenotype, allele, physical_interaction,
    gene_genetic_interaction, orthology/paralogy, RNASeq expression family.
    """
    n = filter_by_column_in_set(
        src("fbgn_fbtr_fbpp_expanded"), dst("fbgn_fbtr_fbpp_expanded"),
        col=2, id_set=id_sets["gene"], min_cols=10,
    )
    record("fbgn_fbtr_fbpp_expanded", n)

    n = filter_by_column_in_set(
        src("fbgn_uniprot"), dst("fbgn_uniprot"), col=0, id_set=id_sets["gene"], min_cols=4,
    )
    record("fbgn_uniprot", n)

    n = filter_by_column_in_set(
        src("fbal_to_fbgn"), dst("fbal_to_fbgn"), col=2, id_set=id_sets["gene"], min_cols=4,
    )
    record("fbal_to_fbgn", n)

    n = filter_mitab_by_fbgn(src("physical_interactions_mitab"), dst("physical_interactions_mitab"), id_sets["gene"])
    if n == 0:
        info = write_synthetic_tsv(
            dst("physical_interactions_mitab"),
            num_columns=42, id_columns_values={0: sorted(id_sets["gene"]), 1: sorted(id_sets["gene"])},
        )
        record("physical_interactions_mitab", n, synthetic_info=info)
    else:
        record("physical_interactions_mitab", n)

    n = filter_by_column_in_set(
        src("dmel_human_orthologs_disease"), dst("dmel_human_orthologs_disease"),
        col=0, id_set=id_sets["gene"], min_cols=6,
    )
    record("orthology_association", n)

    n = filter_by_two_columns_in_set(
        src("dmel_paralogs"), dst("dmel_paralogs"), col_a=0, col_b=5, id_set=id_sets["gene"],
    )
    if n == 0:
        info = write_synthetic_tsv(
            dst("dmel_paralogs"),
            num_columns=11, id_columns_values={0: sorted(id_sets["gene"]), 5: sorted(id_sets["gene"])[::-1]},
        )
        record("paralogy_association", n, synthetic_info=info)
    else:
        record("paralogy_association", n)

    n = filter_gene_genetic_interactions(
        src("gene_genetic_interactions"), dst("gene_genetic_interactions"), id_sets["gene"],
    )
    if n == 0:
        info = write_synthetic_tsv(
            dst("gene_genetic_interactions"),
            num_columns=6, id_columns_values={1: sorted(id_sets["gene"]), 3: sorted(id_sets["gene"])[::-1]},
        )
        record("gene_genetic_interactions", n, synthetic_info=info)
    else:
        record("gene_genetic_interactions", n)

    n = filter_by_fbal_intersection(
        src("allele_genetic_interactions"), dst("allele_genetic_interactions"), id_sets["allele"], text_col=1,
    )
    if n == 0:
        info = write_synthetic_tsv(
            dst("allele_genetic_interactions"), num_columns=4, id_columns_values={1: sorted(id_sets["allele"])},
        )
        record("allele_genetic_interactions", n, synthetic_info=info)
    else:
        record("allele_genetic_interactions", n)

    for label, source_name in [
        ("gene_group", "gene_group_data"),
        ("signaling_pathway_group", "signaling_pathway_group_data"),
        ("metabolic_pathway_group", "metabolic_pathway_group_data"),
    ]:
        n = filter_gene_group(src(source_name), dst(source_name), id_sets["gene"])
        record(label, n)

    gene_groups_hgnc_src = src("gene_groups_hgnc")
    gene_groups_hgnc_dst = dst("gene_groups_hgnc")
    gene_groups_hgnc_dst.parent.mkdir(parents=True, exist_ok=True)
    gene_groups_hgnc_dst.write_bytes(gene_groups_hgnc_src.read_bytes())
    record("gene_groups_HGNC_aux_copy", sum(1 for _ in open_maybe_gzip(gene_groups_hgnc_src)))

    n = filter_by_column_in_set(
        src("gene_sequence_ontology"), dst("gene_sequence_ontology"),
        col=0, id_set=id_sets["gene"], min_cols=4,
    )
    record("gene_sequence_ontology", n)

    n = filter_by_column_in_set(
        src("disease_model_annotations"), dst("disease_model_annotations"),
        col=0, id_set=id_sets["gene"], min_cols=7,
    )
    if n == 0:
        info = write_synthetic_tsv(
            dst("disease_model_annotations"),
            num_columns=12,
            id_columns_values={0: sorted(id_sets["gene"]), 6: sorted(id_sets["allele"]) or ["FBal0000000"]},
        )
        record("disease_model", n, synthetic_info=info)
    else:
        record("disease_model", n)

    n = filter_by_fbal_intersection(
        src("genotype_phenotype_data"), dst("genotype_phenotype_data"), id_sets["allele"], text_col=1,
    )
    if n == 0:
        info = write_synthetic_tsv(
            dst("genotype_phenotype_data"),
            num_columns=7,
            id_columns_values={1: sorted(id_sets["allele"]) or ["FBal0000000"]},
        )
        record("genotype_phenotype", n, synthetic_info=info)
    else:
        record("genotype_phenotype", n)

    fbrf_src = src("fbrf_pmid_pmcid_doi")
    fbrf_dst = dst("fbrf_pmid_pmcid_doi")
    fbrf_dst.parent.mkdir(parents=True, exist_ok=True)
    fbrf_dst.write_bytes(fbrf_src.read_bytes())
    record("fbrf_aux_copy", sum(1 for _ in open_maybe_gzip(fbrf_src)))

    n = filter_by_column_in_set(
        src("scrna_seq_gene_expression"), dst("scrna_seq_gene_expression"),
        col=11, id_set=id_sets["gene"], min_cols=12,
    )
    record("scrna_seq_gene_expression", n)

    n = filter_by_column_in_set(
        src("high_throughput_gene_expression"), dst("high_throughput_gene_expression"),
        col=5, id_set=id_sets["gene"], min_cols=6,
    )
    record("high_throughput_gene_expression", n)

    n = filter_by_column_in_set(
        src("gene_rpkm_report"), dst("gene_rpkm_report"), col=1, id_set=id_sets["gene"], min_cols=12,
    )
    record("gene_rpkm_report", n)

    n = filter_by_column_in_set(src("fca2_fbgn_gene"), dst("fca2_fbgn_gene"), col=0, id_set=id_sets["gene"])
    record("fca2_fbgn_gene", n, empty_ok_reason="no fca2 gene-level expression rows for genes in set")

    n = filter_by_column_in_set(
        src("fca2_fbgn_transcript_gene"), dst("fca2_fbgn_transcript_gene"), col=0, id_set=id_sets["gene"],
    )
    record("fca2_fbgn_transcript_gene", n, empty_ok_reason="no fca2 transcript-gene expression rows for genes in set")

    n = filter_by_column_in_set(src("fca2_fbgn_mir_gene"), dst("fca2_fbgn_mir_gene"), col=0, id_set=id_sets["gene"])
    record("fca2_fbgn_mir_gene", n, empty_ok_reason="no fca2 miRNA gene-level expression rows for genes in set")

    n = filter_by_column_in_set(
        src("fca2_fbgn_mir_transcript"), dst("fca2_fbgn_mir_transcript"), col=0, id_set=id_sets["gene"],
    )
    record("fca2_fbgn_mir_transcript", n, empty_ok_reason="no fca2 miRNA transcript expression rows for genes in set")

    n = filter_afca(src("afca_annotation"), dst("afca_annotation"), dst("fbgn_fbtr_fbpp_expanded"), id_sets["gene"])
    record("afca_annotation", n, empty_ok_reason="no afca rows whose gene symbol resolves into genes in set")

    n = filter_eva_snps(src("eva_vcf"), dst("eva_vcf"), dst("gtf"))
    record("eva_snp", n, empty_ok_reason="no EVA variant rows overlap the closure's gene bodies")


def _generate_phase2_shared(species, input_dir, output_dir, cfg, id_sets, record, src, dst):
    """Phase 2 "round 1" sources that apply the same way across every
    species: bgee, alliance (gene-disease + gene-orthology), epd, and
    reactome (pathway/reaction nodes + gene/transcript/protein/small-molecule
    linkage, one hop of pathway hierarchy, PPI, reaction_to_pathway +
    protein-role edges via reactome_reaction_exporter_All_species.txt, and
    the Pathways2GoTerms/Reactions2GoTerms GO cross-links).
    """
    from biocypher_metta.adapters import Adapter
    from biocypher_metta.processors.ensembl_uniprot_processor import EnsemblUniProtProcessor

    taxon_id = cfg["taxon_id"]
    species_info = Adapter.SPECIES_INFO[taxon_id]

    n = filter_bgee(src("bgee"), dst("bgee"), id_sets["gene"])
    record("bgee", n, empty_ok_reason="no bgee expression coverage for genes in set")

    if cfg.get("has_alliance_data", True):
        n = filter_alliance_disease(src("alliance_disease"), dst("alliance_disease"), taxon_id, id_sets["gene"])
        record("alliance_disease", n, empty_ok_reason="no Alliance disease associations for genes in set")

        n = filter_alliance_orthology(src("alliance_orthology"), dst("alliance_orthology"), taxon_id, id_sets["gene"])
        record("alliance_orthology", n, empty_ok_reason="no Alliance orthology rows for genes in set")
    else:
        logger.info("  [alliance] no real data available for %s yet — skipped", species)

    epd_hgnc_map = "./aux_files/dmel/flybase_synonym_mapping.pkl" if species == "dmel" else None
    n = filter_epd(src("epd"), dst("epd"), taxon_id, id_sets["gene"], hgnc_to_ensembl_map=epd_hgnc_map)
    record("epd", n, empty_ok_reason="no EPD promoters resolved to genes in set")

    # --- reactome (scoped core) ---
    processor = EnsemblUniProtProcessor(
        organism=species_info["ensembl_uniprot_organism"],
        cache_dir=species_info["ensembl_uniprot_cache_directory"],
        update_interval_hours=species_info["update_interval_hours"],
    )
    processor.load_or_update()
    ensembl_uniprot_map = processor.mapping

    species_full_name = species_info["full_name"]
    reactome_prefix = species_info["reactome_prefix"]

    n, pathway_ids = filter_reactome_gene_pathway_or_reaction(
        src("reactome_all_levels"), dst("reactome_all_levels"), species_full_name,
        id_sets["gene"], id_sets["transcript"], id_sets["uniprot"], ensembl_uniprot_map,
    )
    record("reactome_all_levels", n, empty_ok_reason="no Reactome pathway membership for genes in set")

    n, reaction_ids = filter_reactome_gene_pathway_or_reaction(
        src("reactome_reactions"), dst("reactome_reactions"), species_full_name,
        id_sets["gene"], id_sets["transcript"], id_sets["uniprot"], ensembl_uniprot_map,
    )
    record("reactome_reactions", n, empty_ok_reason="no Reactome reaction membership for genes in set")

    n, retained_pathway_ids = filter_reactome_pathways_relation(
        src("reactome_pathways_relation"), dst("reactome_pathways_relation"), reactome_prefix, pathway_ids,
    )
    record("reactome_pathways_relation", n, empty_ok_reason="no pathway-hierarchy rows touching referenced pathways")

    n = filter_reactome_pathways_nodes(
        src("reactome_pathways"), dst("reactome_pathways"), species_full_name, retained_pathway_ids,
    )
    record("reactome_pathways", n, empty_ok_reason="no referenced pathways found")

    # Ensembl2ReactomeReactions.txt doubles as the reaction *node* source too
    # (reactome_reaction adapter) — reaction ids referenced only need to
    # exist in the file we already wrote above (dst("reactome_reactions")),
    # no separate node-only file required.
    n = filter_by_column_in_set(
        src("reactome_reaction_pmids"), dst("reactome_reaction_pmids"), col=0, id_set=reaction_ids, min_cols=2,
    )
    record("reactome_reaction_pmids", n, empty_ok_reason="no PubMed refs for referenced reactions")

    n = filter_reactome_small_molecule(
        src("reactome_chebi_pathways"), dst("reactome_chebi_pathways"), species_full_name, retained_pathway_ids,
    )
    record("reactome_chebi_pathways", n, empty_ok_reason="no small-molecule rows for referenced pathways")

    n = filter_reactome_small_molecule(
        src("reactome_chebi_reactions"), dst("reactome_chebi_reactions"), species_full_name, reaction_ids,
    )
    record("reactome_chebi_reactions", n, empty_ok_reason="no small-molecule rows for referenced reactions")

    # --- uniprot_chebi_part_of_chebi: real BINDING-feature ligand_part_id
    # annotations do exist in the full uniprot_sprot_human.dat.gz (378 across
    # 177 accessions), but none of those 177 accessions are in hsa's own
    # 180-gene closure (confirmed: zero overlap) — a scoping gap, not a
    # UniProt data-availability gap. Either way, the only way to exercise
    # this edge type in the sample is a synthetic record, built from two real
    # CHEBI ids our closure's own Reactome small-molecule rows already
    # reference (so it links two real small_molecule nodes).
    chebi_ids = []
    for chebi_source in ("reactome_chebi_pathways", "reactome_chebi_reactions"):
        with open_maybe_gzip(dst(chebi_source)) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if parts and parts[0] and parts[0] not in chebi_ids:
                    chebi_ids.append(parts[0])
                if len(chebi_ids) >= 2:
                    break
        if len(chebi_ids) >= 2:
            break
    if len(chebi_ids) >= 2:
        synthetic_path, n = write_synthetic_uniprot_chebi_part_of(
            output_dir / "uniprot" / "chebi_part_of_synthetic.dat.gz", chebi_ids[0], chebi_ids[1], taxon_id,
        )
        record("uniprot_chebi_part_of_chebi", 0, synthetic_info=(synthetic_path, n))
    else:
        record("uniprot_chebi_part_of_chebi", 0, empty_ok_reason="no referenced CHEBI ids to build even a synthetic pair from")

    n = filter_reactome_ppi(src("reactome_ppi"), dst("reactome_ppi"), reactome_prefix, id_sets["gene"])
    record("reactome_ppi", n, empty_ok_reason="no PPI rows for genes in set")

    # reactome_reaction_exporter_All_species.txt is genuinely multi-species
    # (from an as-yet-unmerged PR's re-export — see reference in project
    # memory), so this is shared across every species, not hsa-only.
    n = filter_reactome_reaction_exporter(
        src("reactome_reaction_exporter"), dst("reactome_reaction_exporter"), retained_pathway_ids, reaction_ids,
    )
    record("reactome_reaction_exporter", n, empty_ok_reason="no reaction-exporter rows for referenced pathways/reactions")

    logger.info("  [reactome] %d pathways, %d reactions referenced by gene closure", len(retained_pathway_ids), len(reaction_ids))

    # --- rna_central (gene-overlap linking, Phase 3) ---
    bed_kept, rfam_kept = filter_rna_central(
        src("rna_central_bed"), src("rna_central_rfam"),
        dst("rna_central_bed"), dst("rna_central_rfam"),
        dst("gtf"), taxon_id,
    )
    record("rna_central_bed", bed_kept, empty_ok_reason="no ncRNA overlapping genes in set")
    record("rna_central_rfam", rfam_kept, empty_ok_reason="no GO annotations for the retained ncRNAs")


def _generate_hsa_phase2_extension(species, input_dir, output_dir, id_sets, record, src, dst):
    """hsa-only Phase 2/3 additions: GWAS (direct gene-ID edges), HPO's two
    gene-linking adapters (gene<->phenotype, gene<->disease, both keyed by
    bare Entrez ID), TADmap, TFBS, and ABC (+ a scoped dbsnp_rsid_map cache
    built alongside it). HPO's bare ontology-node adapter and GWAS's
    cCRE-adjacent context aren't sampled here — no source filtering needed
    for chebi_ontology/hpo_human_phenotype_ontology (small standalone
    ontologies, wired directly into the sample adapters config exactly like
    GO already is).
    """
    from biocypher_metta.processors import HGNCProcessor

    n = filter_by_any_column_in_set(
        src("gwas"), dst("gwas"),
        cols=[15, 16, 17], id_set=id_sets["gene"], split_on=";", min_cols=18,
    )
    record("gwas", n, empty_ok_reason="no GWAS SNP-gene associations for genes in set")

    n = filter_by_column_in_set_with_header(
        src("hpo_gene_phenotype"), dst("hpo_gene_phenotype"), col=0, id_set=id_sets["entrez"], min_cols=6,
    )
    record("hpo_gene_phenotype", n, empty_ok_reason="no HPO phenotype rows for genes in set")

    n = filter_by_column_in_set_with_header(
        src("hpo_gene_disease"), dst("hpo_gene_disease"), col=0, id_set=id_sets["entrez"], min_cols=5,
        strip_prefix="NCBIGene:",
    )
    record("hpo_gene_disease", n, empty_ok_reason="no HPO disease rows for genes in set")

    n = filter_tadmap(src("tadmap"), dst("tadmap"), 9606, id_sets["gene"])
    record("tadmap", n, empty_ok_reason="no TAD regions overlapping genes in set")

    hgnc_processor = HGNCProcessor()
    hgnc_processor.load_or_update()

    n = filter_tfbs(src("tfbs"), dst("tfbs"), id_sets["gene"], hgnc_processor)
    record("tfbs", n, empty_ok_reason="no TFBS regions resolved to genes in set")

    # --- ABC + CADD + RefSeq share one scoped dbsnp cache (forward + reverse
    # lookup, built from their own rsid/chr/pos columns — no VCF needed) ---
    n_abc, abc_triples = filter_abc(src("abc"), dst("abc"), id_sets["gene"], hgnc_processor)
    record("abc", n_abc, empty_ok_reason="no ABC enhancer-gene rows for genes in set")

    n_cadd, cadd_triples = filter_cadd(src("cadd"), dst("cadd"), chr_filter="chr16")
    record("cadd", n_cadd, empty_ok_reason="no chr16 CADD rows")

    n_refseq, refseq_triples = filter_refseq(src("refseq_closest_gene"), dst("refseq_closest_gene"), id_sets["gene"], hgnc_processor)
    record("refseq_closest_gene", n_refseq, empty_ok_reason="no RefSeq closest-gene rows for genes in set")

    dbsnp_cache_rsids = {t[0] for t in abc_triples + cadd_triples + refseq_triples}
    n_rsids = build_dbsnp_cache(abc_triples + cadd_triples + refseq_triples, "aux_files/hsa/sample_dbsnp")
    logger.info("  [dbsnp cache] built from ABC+CADD+RefSeq rows: %d distinct rsids", n_rsids)

    for pop in ("afr", "eas", "eur", "sas"):
        name = f"topld_{pop}"
        n = filter_topld(src(name), dst(name))
        record(name, n, empty_ok_reason="no chr16:53M-56M LD pairs")

    n = filter_gtex(src("gtex_forgedb"), dst("gtex_forgedb"), id_sets["gene"])
    record("gtex_forgedb", n, empty_ok_reason="no GTEx eQTL/expression rows for genes in set")

    hocomoco_pwm_input = src("hocomoco_annotation").parent / "pwm"
    hocomoco_pwm_output = dst("hocomoco_annotation").parent / "pwm"
    n_models, n_pwm = filter_hocomoco(
        src("hocomoco_annotation"), dst("hocomoco_annotation"),
        hocomoco_pwm_input, hocomoco_pwm_output, id_sets["gene"], hgnc_processor,
    )
    record("hocomoco_annotation", n_models, empty_ok_reason="no HOCOMOCO TF motifs for genes in set")
    logger.info("  [hocomoco] copied %d matching .pwm files", n_pwm)

    n = filter_dbsuper(src("dbsuper"), dst("dbsuper"), id_sets["gene"], hgnc_processor)
    record("dbsuper", n, empty_ok_reason="no dbSuper super-enhancers for genes in set")

    n_link, n_enh, n_src = filter_peregrine(
        src("peregrine_enhancers"), src("peregrine_sources"), src("peregrine_gene_link"),
        dst("peregrine_enhancers"), dst("peregrine_sources"), dst("peregrine_gene_link"),
        id_sets["gene"], hgnc_processor,
    )
    record("peregrine_gene_link", n_link, empty_ok_reason="no Peregrine enhancer-gene links for genes in set")
    logger.info("  [peregrine] %d enhancers, %d sources kept", n_enh, n_src)

    gwas_rsids = collect_column_values(dst("gwas"), col=21, min_cols=22, delimiter="\t")
    gtex_rsids = collect_column_values(dst("gtex_forgedb"), col=0, min_cols=22, delimiter=",", skip_header=True)
    # SNP nodes must cover every rsid referenced as a source by ABC/CADD/RefSeq/
    # TopLD's edges too (source: snp in the schema for activity_by_contact,
    # closest_gene, upstream_of/downstream_of, eqtl_association) -- not just
    # GWAS/GTEx -- or the vast majority of those edges dangle on import.
    all_snp_rsids = gwas_rsids | gtex_rsids | dbsnp_cache_rsids
    n = filter_dbsnp_snps(src("dbsnp_common_vcf"), dst("dbsnp_common_vcf"), all_snp_rsids)
    record("dbsnp_common_vcf", n, empty_ok_reason="no dbSNP rows for rsids referenced by genes in set")

    n_nodes, n_edges = filter_enhancer_atlas(
        src("enhancer_atlas_bed"), dst("enhancer_atlas_bed"),
        input_dir / "enhancer_atlas" / "enhancer_gene", output_dir / "enhancer_atlas" / "enhancer_gene",
        id_sets["gene"],
    )
    record("enhancer_atlas_bed", n_nodes, empty_ok_reason="no EnhancerAtlas regions referencing genes in set")
    logger.info("  [enhancer_atlas] %d enhancer-gene edge rows kept across tissue files", n_edges)

    n = filter_ccre_closest_genes(src("ccre_closest_genes_all"), dst("ccre_closest_genes_all"), id_sets["gene"])
    record("ccre_closest_genes_all", n, empty_ok_reason="no cCRE closest-gene rows (All) for genes in set")

    n = filter_ccre_closest_genes(src("ccre_closest_genes_pc"), dst("ccre_closest_genes_pc"), id_sets["gene"])
    record("ccre_closest_genes_pc", n, empty_ok_reason="no cCRE closest-gene rows (PC) for genes in set")

    n = filter_ccre_eqtl(src("ccre_eqtl_gene_links"), dst("ccre_eqtl_gene_links"), id_sets["gene"])
    record("ccre_eqtl_gene_links", n, empty_ok_reason="no cCRE eQTL gene-link rows for genes in set")

    # --- Catlas: filter ABC_scores by gene, then cCRE master + per-cell-type
    # .bed files down to just the coordinates ABC_scores actually referenced.
    # catlas_ccre_label_map.pkl is coordinate-dependent, so it's built at a
    # sample-specific path, NOT the shared aux_files/hsa/catlas/ location.
    catlas_input_dir = input_dir / "catlasv1"
    catlas_output_dir = output_dir / "catlasv1"
    n_abc_files, n_abc_rows, n_master, n_bed_files = filter_catlas(
        catlas_input_dir / "ABC_scores", catlas_output_dir / "ABC_scores",
        catlas_input_dir / "cCRE_hg38.tsv.gz", catlas_output_dir / "cCRE_hg38.tsv.gz",
        catlas_input_dir / "cCREs", catlas_output_dir / "cCREs",
        Path("aux_files/hsa/catlas_connected_sample/catlas_ccre_label_map.pkl"),
        id_sets["gene"], hgnc_processor,
    )
    record("catlas_abc_scores", n_abc_rows, empty_ok_reason="no Catlas ABC-score rows for genes in set")
    logger.info(
        "  [catlas] %d/112 ABC_scores files, %d master cCREs, %d/222 per-cell-type .bed files kept",
        n_abc_files, n_master, n_bed_files,
    )
    cell_ontology_src = catlas_input_dir / "Cell_ontology.tsv"
    cell_ontology_dst = catlas_output_dir / "Cell_ontology.tsv"
    cell_ontology_dst.parent.mkdir(parents=True, exist_ok=True)
    cell_ontology_dst.write_bytes(cell_ontology_src.read_bytes())

    # --- motif_diff: no dbsnp_rsid_map support at all — filter by every rsid
    # already known relevant to our gene closure from other sources. Capped
    # to a small representative sample (confirmed 2026-08-01: alters_binding
    # contributes zero essential gene-closure connectivity — pure schema/
    # coverage breadth, and each row is ~10.8KB wide with 770 TF-model score
    # columns, so row count alone drives most of the ~10GB this source used
    # to occupy uncapped).
    n = filter_motif_diff(src("motif_diff"), dst("motif_diff"), all_snp_rsids, max_rows=2000)
    record("motif_diff", n, empty_ok_reason="no motif-diff rows for rsids referenced by genes in set")

    # --- Roadmap Epigenomics (chromatin_state/h3_marks/dhs): rsid-keyed, no
    # chr/pos of their own -- the real adapters resolve rsid -> chr/pos via
    # dbsnp_rsid_map at runtime and then apply the config's fixed
    # chr16:53M-56M window (same window CADD/TopLD use), so pre-filtering to
    # rsids already known (from the same ABC/CADD/RefSeq cache) to fall in
    # that window reproduces exactly what the adapters would keep, without
    # keeping the full ~20GB of source files.
    roadmap_window_rsids = {
        t[0] for t in abc_triples + cadd_triples + refseq_triples
        if t[1] == "chr16" and 53_000_000 <= int(t[2]) <= 56_000_000
    }
    n_state = filter_roadmap_dir(
        input_dir / "forgedb" / "roadmap" / "chromatin_state",
        output_dir / "forgedb" / "roadmap" / "chromatin_state",
        roadmap_window_rsids,
    )
    record("roadmap_chromatin_state", n_state, empty_ok_reason="no chromatin-state rows for chr16:53M-56M rsids in set")

    n_h3 = filter_roadmap_dir(
        input_dir / "forgedb" / "roadmap" / "h3_marks",
        output_dir / "forgedb" / "roadmap" / "h3_marks",
        roadmap_window_rsids,
    )
    record("roadmap_h3_mark", n_h3, empty_ok_reason="no H3-mark rows for chr16:53M-56M rsids in set")

    n_dhs = filter_roadmap_file(src("roadmap_dhs"), dst("roadmap_dhs"), roadmap_window_rsids)
    record("roadmap_dhs", n_dhs, empty_ok_reason="no DHS rows for chr16:53M-56M rsids in set")

    # --- Reactome GO cross-links (Pathways2GoTerms_human.txt /
    # Reactions2GoTerms_human.txt): hsa-only -- these files have no
    # per-species variant on reactome.org (unlike reaction_exporter, which
    # got a genuinely multi-species re-export), and the other 4 species'
    # full configs don't define reactome_pathway_to_biological_process/
    # reaction_to_molecular_function at all. Reconstruct the pathway/reaction
    # id sets from the already-written filtered files (_generate_phase2_shared
    # computed them but doesn't pass them through to this hsa-only function).
    retained_pathway_ids = set()
    with open_maybe_gzip(dst("reactome_pathways")) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                retained_pathway_ids.add(parts[0])

    reaction_ids = set()
    with open_maybe_gzip(dst("reactome_reactions")) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                reaction_ids.add(parts[1])

    n = filter_reactome_go_terms(src("reactome_pathways_go_bp"), dst("reactome_pathways_go_bp"), retained_pathway_ids)
    record("reactome_pathways_go_bp", n, empty_ok_reason="no GO biological_process terms for referenced pathways")

    n = filter_reactome_go_terms(src("reactome_reactions_go_mf"), dst("reactome_reactions_go_mf"), reaction_ids)
    record("reactome_reactions_go_mf", n, empty_ok_reason="no GO molecular_function terms for referenced reactions")


def generate_sample(species, input_dir, output_dir, size_budget, anchor_genes_file, interactive=True):
    cfg = SPECIES[species]
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    if not ensure_required_files(species, input_dir, interactive=interactive):
        return None

    logger.info("Loading anchor genes from %s", anchor_genes_file)
    anchor_gene_ids = load_anchor_genes(anchor_genes_file, cfg["anchor_id_key"])
    logger.info("Loaded %d anchor genes", len(anchor_gene_ids))

    logger.info("Building backbone maps from real source files under %s...", input_dir)
    maps = build_backbone_maps(input_dir, species, cfg)

    logger.info("Expanding gene set by cross-reference coverage (budget=%d)...", size_budget)
    id_sets = expand_gene_set(input_dir, species, anchor_gene_ids, maps, size_budget)
    for key, ids in id_sets.items():
        logger.info("  id_set[%s] = %d", key, len(ids))

    manifest = {
        "species": species,
        "input_dir": str(input_dir),
        "anchor_genes": sorted(anchor_gene_ids),
        "expanded_gene_count": len(id_sets["gene"]),
        "size_budget": size_budget,
        "sources": {},
    }

    def record(name, real_rows, empty_ok_reason=None, synthetic_info=None):
        entry = {"real_rows": real_rows, "synthetic": synthetic_info is not None}
        if synthetic_info:
            entry["synthetic_rows"] = synthetic_info[1]
            entry["synthetic_file"] = str(synthetic_info[0])
        if real_rows == 0 and synthetic_info is None and empty_ok_reason:
            entry["reason"] = empty_ok_reason
        manifest["sources"][name] = entry
        logger.info(
            "  [%s] real_rows=%d%s", name, real_rows,
            f" (+{synthetic_info[1]} synthetic)" if synthetic_info else "",
        )

    def src(name):
        return resolve_source_file(input_dir, species, name)

    def dst(name):
        return output_path_for(output_dir, species, name)

    _generate_core(species, input_dir, output_dir, cfg, id_sets, record, src, dst)
    if cfg["is_flybase"]:
        _generate_dmel_extension(species, input_dir, output_dir, id_sets, record, src, dst)
    _generate_phase2_shared(species, input_dir, output_dir, cfg, id_sets, record, src, dst)
    if species == "hsa":
        _generate_hsa_phase2_extension(species, input_dir, output_dir, id_sets, record, src, dst)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "sample_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info("Wrote manifest: %s", manifest_path)
    return manifest


def _configured_species_names():
    species_config_path = Path(__file__).resolve().parent.parent / "config" / "species_config.yaml"
    with open(species_config_path) as f:
        return sorted(yaml.safe_load(f).keys())


def main():
    species_names = _configured_species_names()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", choices=species_names + ["all"], default=None)
    parser.add_argument("--input-dir", default=None, help="Real, full-size source data dir (long-lived mirror, or a download_data.py output dir). Ignored with --species all.")
    parser.add_argument("--output-dir", default=None, help="Where to write the filtered sample. Ignored with --species all.")
    parser.add_argument("--size-budget", type=int, default=180)
    parser.add_argument("--anchor-genes-file", default=None, help="Ignored with --species all.")
    args = parser.parse_args()

    if args.species is None:
        print("Usage: generate_connected_sample.py --species {" + ",".join(species_names + ["all"]) + "} [options]")
        print()
        print("Available species (from config/species_config.yaml):")
        for name in species_names:
            print(f"  {name}")
        print("  all  (generate every species above, one after another)")
        return

    if args.species == "all":
        if args.input_dir or args.output_dir or args.anchor_genes_file:
            parser.error("--input-dir/--output-dir/--anchor-genes-file are per-species and can't be combined with --species all")
        for name in species_names:
            cfg = SPECIES[name]
            generate_sample(name, species_input_dir(name), cfg["default_output_dir"], args.size_budget, cfg["anchor_genes_file"])
        return

    cfg = SPECIES[args.species]
    input_dir = args.input_dir or species_input_dir(args.species)
    output_dir = args.output_dir or cfg["default_output_dir"]
    anchor_genes_file = args.anchor_genes_file or cfg["anchor_genes_file"]

    generate_sample(args.species, input_dir, output_dir, args.size_budget, anchor_genes_file)


if __name__ == "__main__":
    main()
