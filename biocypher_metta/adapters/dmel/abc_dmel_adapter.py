# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher_metta.adapters import Adapter
import csv
import gzip
from pathlib import Path
import pickle
from typing import Dict, List, Optional, Union
from biocypher._logger import logger
from biocypher_metta.adapters.helpers import check_genomic_location, to_float

# Expected ABC Dmel column layout:
# rsid, chromosome, position, element_chr, element_start, element_end, element_name, target_gene, target_gene_tss, abc_score, cell_type
DEFAULT_COL_DICT = {
    'rsid': 0,
    'chromosome': 1,
    'position': 2,
    'element_chr': 3,
    'element_start': 4,
    'element_end': 5,
    'element_name': 6,
    'target_gene': 7,
    'target_gene_tss': 8,
    'abc_score': 9,
    'cell_type': 10,
}

# Default tissue / cell type to FlyBase Anatomy Ontology (FBbt) mapping:
DEFAULT_TISSUE_TO_ONTOLOGY_ID_MAP = {
    'Kc167': 'FBbt:00005349',
    'S2': 'FBbt:00005737',
}


class ABCDmelAdapter(Adapter):
    """
    Adapter for Drosophila melanogaster Activity-By-Contact (ABC) prediction data.
    Converts variant-to-gene ABC prediction records into 'activity_by_contact' edges
    connecting snp nodes (e.g. RS...) to target genes (e.g. FBgn...).
    """

    def __init__(
        self,
        filepath: Union[str, Path, List[Union[str, Path]]],
        label: str = 'activity_by_contact',
        taxon_id: int = 7227,
        tissue_to_ontology_id_map: Optional[Union[Dict[str, str], str, Path]] = None,
        threshold: Optional[float] = None,
        gene_prefix: str = '',
        variant_prefix: str = '',
        variant_case: str = 'upper',
        chr: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        write_properties: Optional[bool] = None,
        add_provenance: Optional[bool] = None,
    ):
        if isinstance(filepath, (str, Path)):
            self.filepaths = [Path(filepath)]
        elif isinstance(filepath, (list, tuple)):
            self.filepaths = [Path(p) for p in filepath]
        else:
            self.filepaths = [Path(filepath)]

        self.label = label
        self.taxon_id = int(taxon_id) if taxon_id is not None else 7227
        self.threshold = threshold
        self.gene_prefix = gene_prefix or ''
        self.variant_prefix = variant_prefix or ''
        self.variant_case = variant_case.lower() if isinstance(variant_case, str) else 'upper'
        self.chr = chr
        self.start = start
        self.end = end

        # Resolve tissue mapping (dict, pickle path, or default)
        if tissue_to_ontology_id_map is None:
            self.tissue_to_ontology_id_map = DEFAULT_TISSUE_TO_ONTOLOGY_ID_MAP.copy()
        elif isinstance(tissue_to_ontology_id_map, dict):
            self.tissue_to_ontology_id_map = dict(tissue_to_ontology_id_map)
        elif isinstance(tissue_to_ontology_id_map, (str, Path)):
            with open(tissue_to_ontology_id_map, 'rb') as fp:
                self.tissue_to_ontology_id_map = pickle.load(fp)
        else:
            self.tissue_to_ontology_id_map = dict(tissue_to_ontology_id_map)

        self.source = 'ABC'
        self.source_url = 'https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction'
        self.version = 'r6.69-dm6'
        super(ABCDmelAdapter, self).__init__(write_properties, add_provenance)

    def _normalize_variant_id(self, raw_rsid: str) -> str:
        s = raw_rsid.strip()
        if self.variant_case == 'upper':
            s = s.upper()
        elif self.variant_case == 'lower':
            s = s.lower()
        if self.variant_prefix and not s.startswith(self.variant_prefix):
            s = f"{self.variant_prefix}{s}"
        return s

    def _normalize_gene_id(self, raw_gene: str) -> str:
        s = raw_gene.strip()
        if self.gene_prefix and not s.startswith(self.gene_prefix):
            s = f"{self.gene_prefix}{s}"
        return s

    def _map_cell_type(self, raw_cell_type: str) -> Optional[str]:
        raw_cell_type = raw_cell_type.strip()
        if not raw_cell_type:
            return None

        # Direct dictionary match
        if raw_cell_type in self.tissue_to_ontology_id_map:
            return self.tissue_to_ontology_id_map[raw_cell_type]

        # Comma-separated cell types
        if ',' in raw_cell_type:
            parts = [p.strip() for p in raw_cell_type.split(',') if p.strip()]
            mapped = [self.tissue_to_ontology_id_map.get(p) for p in parts]
            if all(m is not None for m in mapped):
                return ','.join(mapped) if len(mapped) > 1 else mapped[0]
            # If some are mapped, keep mapped values for known ones
            mapped_filtered = [m for m in mapped if m is not None]
            if mapped_filtered:
                return ','.join(mapped_filtered) if len(mapped_filtered) > 1 else mapped_filtered[0]

        return None

    def _matches_chr(self, filter_chr: Optional[str], row_chr: str) -> bool:
        if filter_chr is None:
            return True
        f_chr = str(filter_chr).strip()
        r_chr = str(row_chr).strip()
        if f_chr == r_chr:
            return True
        # Allow comparing e.g. "chr2L" vs "2L"
        f_norm = f_chr[3:] if f_chr.lower().startswith('chr') else f_chr
        r_norm = r_chr[3:] if r_chr.lower().startswith('chr') else r_chr
        return f_norm == r_norm

    def _open_file(self, path: Path):
        path_str = str(path)
        if path_str.endswith('.gz'):
            return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
        return open(path, 'r', encoding='utf-8', errors='replace')

    def get_edges(self):
        """
        Yields edges in the format:
        (variant_id, gene_id, label, {"score": score, "biological_context": context, "taxon_id": taxon_id})
        """
        not_processed = 0
        processed = 0

        for file_path in self.filepaths:
            if not file_path.exists():
                logger.error(f"File not found: {file_path}")
                continue

            try:
                with self._open_file(file_path) as fp:
                    reader = csv.reader(fp, delimiter=',')
                    try:
                        header = next(reader)
                    except StopIteration:
                        # Empty file
                        continue

                    # Build column index mapping from header if possible
                    col_map = DEFAULT_COL_DICT.copy()
                    if header and any(h.strip() in DEFAULT_COL_DICT for h in header):
                        header_stripped = [h.strip() for h in header]
                        for col_name in DEFAULT_COL_DICT:
                            if col_name in header_stripped:
                                col_map[col_name] = header_stripped.index(col_name)

                    min_required_cols = max(col_map['rsid'], col_map['chromosome'],
                                            col_map['position'], col_map['target_gene'],
                                            col_map['abc_score'], col_map['cell_type']) + 1

                    for row in reader:
                        if not row or len(row) < min_required_cols:
                            not_processed += 1
                            continue

                        try:
                            raw_rsid = row[col_map['rsid']].strip()
                            raw_chr = row[col_map['chromosome']].strip()
                            raw_pos = row[col_map['position']].strip()
                            raw_gene = row[col_map['target_gene']].strip()
                            raw_score = row[col_map['abc_score']].strip()
                            raw_cell = row[col_map['cell_type']].strip()

                            if not raw_rsid or not raw_gene:
                                not_processed += 1
                                continue

                            # Genomic location filter on variant position
                            if self.chr is not None or self.start is not None or self.end is not None:
                                if not self._matches_chr(self.chr, raw_chr):
                                    continue
                                try:
                                    pos_int = int(raw_pos)
                                except ValueError:
                                    not_processed += 1
                                    continue
                                if not check_genomic_location(raw_chr, self.start, self.end, raw_chr, pos_int, pos_int):
                                    continue

                            # Parse score
                            try:
                                score = to_float(raw_score)
                            except (ValueError, TypeError):
                                not_processed += 1
                                continue

                            # Threshold filter: keep all when threshold is None
                            if self.threshold is not None and score < self.threshold:
                                continue

                            # Map cell type to biological context
                            ontology_id = self._map_cell_type(raw_cell)
                            if ontology_id is None:
                                not_processed += 1
                                continue

                            source = self._normalize_variant_id(raw_rsid)
                            target = self._normalize_gene_id(raw_gene)

                            props = {
                                'score': score,
                                'biological_context': ontology_id,
                                'taxon_id': self.taxon_id,
                            }
                            if self.write_properties is not None and not self.write_properties:
                                props = {}

                            processed += 1
                            yield source, target, self.label, props

                        except Exception as e:
                            logger.warning(f"Error while parsing row {row}: {e}")
                            not_processed += 1
                            continue

            except Exception as e:
                logger.error(f"Error reading {file_path}: {e}")
                continue

        logger.info(f"ABCDmelAdapter processed {processed} records ({not_processed} skipped).")
