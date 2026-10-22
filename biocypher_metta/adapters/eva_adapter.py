"""
EVA (European Variation Archive) adapter for multi-species variant import.

Supports four non-human species:
  - Drosophila melanogaster (dmel, taxon_id=7227)
  - Caenorhabditis elegans  (cel,  taxon_id=6239)
  - Mus musculus             (mmu,  taxon_id=10090)
  - Rattus norvegicus        (rno,  taxon_id=10116)

Only RS_VALIDATED or SS_VALIDATED (count > 0) variants are imported.

VCF data from:
  https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_9/by_species/

Example VCF record:
  #CHROM  POS   ID            REF  ALT  QUAL  FILTER  INFO
  I       380   rs7089510721  C    A    .     .       SID=PRJEB80000;SS_VALIDATED=0;VC=SO:0001483
  I       392   rs7089703857  C    A    .     .       SID=PRJEB80000;SS_VALIDATED=1;RS_VALIDATED;VC=SO:0001483
"""

import gzip
import numpy as np
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import check_genomic_location
from biocypher._logger import logger

# NCLS is required for interval-overlap queries (gene body detection).

try:
    from ncls import NCLS
    _NCLS_AVAILABLE = True
except ImportError:
    _NCLS_AVAILABLE = False
    logger.warning(
        "ncls package not found — EVA edge generation will fail if enabled. "
        "Install it with: pip install ncls"
    )

# CURIE prefix for EVA variant identifiers
_EVA_PREFIX = "EVA"


def _parse_info(info_string: str) -> dict:
    """Parse a VCF INFO field into a dict.

    Flags (no '=' sign) are stored as True.
    Integer / string values are stored as-is (strings).
    """
    result = {}
    if not info_string or info_string == ".":
        return result

    for entry in info_string.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            result[key] = value
        else:
            result[entry] = True
    return result


def _open_vcf(filepath: str):
    """Open a VCF file whether it is gzip-compressed (.vcf.gz) or plain text (.vcf)."""
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt")
    return open(filepath, "rt")


def _is_validated(info: dict) -> bool:
    """Return True if the variant passes the validation filter.

    A variant is imported if:
      - RS_VALIDATED flag is present (True or != "0"), OR
      - SS_VALIDATED flag is present (True or != "0").
    """
    rs_val = info.get("RS_VALIDATED")
    if rs_val is True or (rs_val is not None and str(rs_val) != "0"):
        return True
        
    ss_val = info.get("SS_VALIDATED")
    if ss_val is True or (ss_val is not None and str(ss_val) != "0"):
        return True
        
    return False


class EVAAdapter(Adapter):
    """Adapter for EVA (European Variation Archive) VCF data.

    Produces SNP nodes and/or 'in gene to variant association' edges.
    """

    # VCF column indices
    _IDX = {"chr": 0, "pos": 1, "id": 2, "ref": 3, "alt": 4, "info": 7}

    # GTF column indices (tab-separated)
    _GTF_IDX = {"chr": 0, "type": 2, "start": 3, "end": 4, "info": 8}

    def __init__(
        self,
        filepath: str,
        taxon_id: int,
        label: str,
        write_properties: bool,
        add_provenance: bool,
        gene_filepath: str = None,
        chr: str = None,
        start: int = None,
        end: int = None,
    ):
        self.filepath = filepath
        self.taxon_id = taxon_id
        self.label = label
        self.gene_filepath = gene_filepath
        self.chr = chr
        self.start = start
        self.end = end

        meta = self.SPECIES_INFO.get(taxon_id, {})
        self.source = f"EVA Release 9 – {meta.get('full_name', str(taxon_id))}"
        self.source_url = meta.get("eva_source_url", "https://www.ebi.ac.uk/eva/")

        # Interval index built lazily in get_edges()
        self._gene_index = None   # dict[chr_str] -> (NCLS, list[gene_id])

        super(EVAAdapter, self).__init__(write_properties, add_provenance)


    def _build_gene_index(self) -> dict:
        """Build a chromosome-keyed interval index from the gene GTF file.

        Returns
        -------
        dict[str, tuple[NCLS, list[str]]]
            Maps chromosome name → (NCLS object, parallel list of gene IDs).
            The NCLS stores 0-based half-open intervals [start-1, end).
        """
        if not _NCLS_AVAILABLE:
            raise RuntimeError(
                "ncls is required for EVA edge generation. "
                "Install it with: pip install ncls"
            )
        if not self.gene_filepath:
            raise ValueError(
                "gene_filepath must be provided for EVA edge generation."
            )

        logger.info(f"EVAAdapter: building gene interval index from {self.gene_filepath}")

        # Accumulate per-chromosome lists before building NCLS
        chr_starts: dict[str, list] = {}
        chr_ends:   dict[str, list] = {}
        chr_ids:    dict[str, list] = {}

        open_fn = gzip.open if str(self.gene_filepath).endswith(".gz") else open

        with open_fn(self.gene_filepath, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[self._GTF_IDX["type"]] != "gene":
                    continue

                chrom = parts[self._GTF_IDX["chr"]]
                try:
                    g_start = int(parts[self._GTF_IDX["start"]]) - 1  
                    g_end = int(parts[self._GTF_IDX["end"]])           
                except ValueError:
                    continue

                gene_id = self._extract_gene_id(parts[self._GTF_IDX["info"]])
                if not gene_id:
                    continue

                chr_starts.setdefault(chrom, []).append(g_start)
                chr_ends.setdefault(chrom, []).append(g_end)
                chr_ids.setdefault(chrom, []).append(gene_id)

        index = {}
        for chrom in chr_starts:
            starts = np.array(chr_starts[chrom], dtype=np.int64)
            ends   = np.array(chr_ends[chrom],   dtype=np.int64)
            ids    = np.arange(len(starts),        dtype=np.int64)
            ncls   = NCLS(starts, ends, ids)
            index[chrom] = (ncls, chr_ids[chrom])

        gene_count = sum(len(v[1]) for v in index.values())
        logger.info(f"EVAAdapter: indexed {gene_count} genes across {len(index)} chromosomes")
        return index

    @staticmethod
    def _extract_gene_id(info_field: str) -> str:
        """Extract the gene_id value from a GTF attributes field.

        Handles both GENCODE/Ensembl style (quoted) and bare values.
        Returns the ID stripped of version suffix (e.g. .1), or empty string.
        """
        for attr in info_field.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                parts = attr.split(None, 1)
                if len(parts) == 2:
                    raw = parts[1].strip().strip('"').strip("'")
                    # Remove version suffix (.N) — keep stable Ensembl/FlyBase ID
                    raw = raw.split(".")[0]
                    return raw
        return ""

    def _query_overlapping_genes(self, chrom: str, pos: int) -> list:
        """Return a list of gene IDs whose body overlaps the given position.

        Parameters
        ----------
        chrom : str
            Chromosome name as it appears in the VCF.
        pos : int
            1-based variant position.
        """
        if self._gene_index is None:
            self._gene_index = self._build_gene_index()

        entry = self._gene_index.get(chrom)
        if entry is None:
            return []

        ncls, gene_ids = entry
        pos0 = pos - 1  # convert to 0-based
        hits = list(ncls.find_overlap(pos0, pos0 + 1))
        return [gene_ids[i] for _, _, i in hits]

    
    def get_nodes(self):
        
        processed = 0
        skipped_validation = 0
        skipped_location = 0

        with _open_vcf(self.filepath) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue

                chrom = parts[self._IDX["chr"]]
                try:
                    pos = int(parts[self._IDX["pos"]])
                except ValueError:
                    continue

                rsid_raw = parts[self._IDX["id"]]
                ref = parts[self._IDX["ref"]]
                alt = parts[self._IDX["alt"]]
                info = _parse_info(parts[self._IDX["info"]])

                if not _is_validated(info):
                    skipped_validation += 1
                    continue

                if not check_genomic_location(
                    self.chr, self.start, self.end, chrom, pos, pos
                ):
                    skipped_location += 1
                    continue

                node_id = f"{_EVA_PREFIX}:{rsid_raw}"

                props = {}
                if self.write_properties:
                    props["chr"] = chrom
                    props["start"] = pos
                    props["end"] = pos
                    props["ref"] = ref
                    props["alt"] = alt
                    props["taxon_id"] = self.taxon_id
                    if self.add_provenance:
                        props["source"] = self.source
                        props["source_url"] = self.source_url

                processed += 1
                yield node_id, self.label, props

        logger.info(
            f"EVAAdapter.get_nodes [{self.taxon_id}]: "
            f"yielded={processed}, "
            f"skipped_validation={skipped_validation}, "
            f"skipped_location={skipped_location}"
        )

    def get_edges(self):
        if not self.gene_filepath:
            raise ValueError(
                "gene_filepath must be set to generate EVA snp-gene edges."
            )

        processed = 0
        skipped_validation = 0
        skipped_no_gene = 0

        with _open_vcf(self.filepath) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue

                chrom = parts[self._IDX["chr"]]
                try:
                    pos = int(parts[self._IDX["pos"]])
                except ValueError:
                    continue

                info = _parse_info(parts[self._IDX["info"]])

                if not _is_validated(info):
                    skipped_validation += 1
                    continue

                if not check_genomic_location(
                    self.chr, self.start, self.end, chrom, pos, pos
                ):
                    continue

                rsid_raw = parts[self._IDX["id"]]
                snp_id = f"{_EVA_PREFIX}:{rsid_raw}"

                overlapping = self._query_overlapping_genes(chrom, pos)
                if not overlapping:
                    skipped_no_gene += 1
                    continue

                for gene_id in overlapping:
                    props = {}
                    if self.write_properties:
                        props["taxon_id"] = self.taxon_id
                        if self.add_provenance:
                            props["source"] = self.source
                            props["source_url"] = self.source_url

                    processed += 1
                    target_id = gene_id
                    if ":" not in gene_id:
                        curie_prefix = self.CURIE_PREFIX.get(self.taxon_id)
                        if curie_prefix:
                            target_id = f"{curie_prefix}:{gene_id}"
                    yield snp_id, target_id, self.label, props

        logger.info(
            f"EVAAdapter.get_edges [{self.taxon_id}]: "
            f"yielded={processed}, "
            f"skipped_validation={skipped_validation}, "
            f"skipped_no_gene={skipped_no_gene}"
        )
