import csv
import os
from biocypher_metta.adapters import Adapter
from biocypher._logger import logger
from biocypher_metta.adapters.helpers import to_float
from biocypher_metta.processors import HGNCProcessor


class MotifDiffAdapter(Adapter):
    """
    Species-agnostic adapter for MotifDiff output (https://github.com/rezwanhosseini/MotifDiff).
    Reads the .diff matrix (rows = variants, columns = motifs) and yields snp -> TF gene edges
    (label tf_snp) with effect (gain/loss), score, motif_id and taxon_id.

    Two ways of resolving a column header to TF genes; everything else is shared:

    * Human / HOCOMOCO (default, motif_map_path=None): headers look like 'CTCF_HUMAN.H11MO.0.A'.
      The TF symbol is the text before '_HUMAN', mapped to an Ensembl id with HGNCProcessor.
      Among several motifs of one TF the one with the best HOCOMOCO confidence letter (A > B > C > D)
      wins, ties broken by the larger |score|. This is the behavior of the former
      adapters/hsa/motif_diff_adapter.py, unchanged.
    * Any species with a motif map (motif_map_path=FILE.tsv): headers are the whole MOTIF line of the
      MEME file (e.g. 'MA0010.1 br'). The TSV (columns motif_name_in_diff_header, matrix_id, symbol,
      uniprot_ids, fbgns, n_partners, status; written by build_dmel_motif_set.py) gives the gene ids
      of each motif. There is no confidence letter, so the larger |score| per gene wins.

    Defaults reproduce the human KG ids (variant 'DBSNP:<rs>', gene 'ENSEMBL:<id>'). For dmel use
    variant_prefix='', uppercase_variant=True, gene_prefix='' (variant 'RS123', gene 'FBgn...').

    heterodimer ('all' or 'skip', map mode only): a motif with several partners (n_partners > 1,
    e.g. EcR::usp) yields one edge per gene, or is dropped.
    """

    def __init__(self, filepath, hgnc_to_ensembl=None, label=None, write_properties=None,
                 add_provenance=None, threshold=1e-3, hgnc_processor=None, taxon_id=9606,
                 motif_map_path=None, heterodimer='all', variant_prefix='DBSNP:',
                 uppercase_variant=False, gene_prefix='ENSEMBL:'):
        if heterodimer not in ('all', 'skip'):
            raise ValueError("heterodimer must be 'all' or 'skip'")
        self.filepath = filepath
        self.taxon_id = taxon_id
        self.label = label
        self.threshold = threshold
        self.heterodimer = heterodimer
        self.variant_prefix = variant_prefix
        self.uppercase_variant = uppercase_variant
        self.gene_prefix = gene_prefix

        self.motifs = None
        if motif_map_path is None:
            # human mode: HGNC symbol -> Ensembl id
            if hgnc_processor is not None:
                self.hgnc_processor = hgnc_processor
            else:
                self.hgnc_processor = HGNCProcessor()
                self.hgnc_processor.load_or_update()
            self._ensembl_cache = {}
        else:
            if not os.path.exists(motif_map_path):
                raise FileNotFoundError(f"Motif map file not found: {motif_map_path}")

            self.hgnc_processor = None
            self.motifs = {}  # exact header -> (matrix_id, [gene ids], n_partners)
            with open(motif_map_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                required_cols = {'motif_name_in_diff_header', 'matrix_id', 'fbgns', 'n_partners'}
                if not reader.fieldnames or not required_cols.issubset(reader.fieldnames):
                    missing = required_cols - set(reader.fieldnames or [])
                    raise ValueError(f"Malformed motif map file '{motif_map_path}': missing required column(s): {sorted(missing)}")
                for row in reader:
                    genes = [g for g in row['fbgns'].split(';') if g]
                    try:
                        n_partners = int(row['n_partners'])
                    except (ValueError, TypeError):
                        n_partners = len(genes)
                    self.motifs[row['motif_name_in_diff_header']] = (
                        row['matrix_id'], genes, n_partners)

        self.source = 'MotifDiff'
        self.source_url = 'https://github.com/rezwanhosseini/MotifDiff'
        super(MotifDiffAdapter, self).__init__(write_properties, add_provenance)

    def parse_tf_id(self, tf_id):
        if '_HUMAN' in tf_id:
            return tf_id.split('_HUMAN')[0]
        return tf_id

    def _resolve_columns(self, header):
        """column index -> (target keys, rank, motif id). A lower rank is better; keys are HGNC symbols
        (human mode) or gene ids (map mode)."""
        cols, unknown, no_gene, skipped = {}, set(), set(), set()
        for i, name in enumerate(header[1:], 1):
            if self.motifs is None:
                tf_id = name
                parts = tf_id.split('.')
                confidence = parts[-1] if parts else 'Z'
                cols[i] = ([self.parse_tf_id(tf_id)], confidence, tf_id)
                continue
            key = name.strip()
            if key not in self.motifs:
                unknown.add(key)
                continue
            matrix_id, genes, n_partners = self.motifs[key]
            if not genes:
                no_gene.add(key)
            elif n_partners > 1 and self.heterodimer == 'skip':
                skipped.add(key)
            else:
                cols[i] = (genes, '', matrix_id)
        if unknown:
            logger.warning(f"{len(unknown)} .diff column(s) not in the motif map, ignored: {sorted(unknown)[:5]}")
        if no_gene:
            logger.warning(f"{len(no_gene)} motif column(s) without a gene id, ignored: {sorted(no_gene)[:5]}")
        if skipped:
            logger.warning(f"{len(skipped)} multi-partner motif column(s) dropped by heterodimer='skip': {sorted(skipped)[:5]}")
        return cols

    def _target_id(self, key):
        if self.motifs is not None:
            return f"{self.gene_prefix}{key}"
        if key not in self._ensembl_cache:
            self._ensembl_cache[key] = self.hgnc_processor.get_ensembl_id(key)
        ensembl = self._ensembl_cache[key]
        if ensembl is None:
            return None  # TF symbol with no Ensembl id: edge skipped, as before
        return f"{self.gene_prefix}{ensembl}"

    def get_edges(self):
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"MotifDiff input file not found: {self.filepath}")

        with open(self.filepath, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader, None)
            if not header:
                return
            cols = self._resolve_columns(header)

            for row in reader:
                if not row:
                    continue
                rsid = row[0].upper() if self.uppercase_variant else row[0]
                variant = f"{self.variant_prefix}{rsid}"

                best = {}  # target key -> (rank, score, motif id)
                for i, score_str in enumerate(row[1:], 1):
                    col = cols.get(i)
                    if col is None:
                        continue
                    try:
                        score = float(score_str)
                    except ValueError:
                        continue
                    if abs(score) < self.threshold:
                        continue
                    keys, rank, motif_id = col
                    for k in keys:
                        cur = best.get(k)
                        # better (lower) rank wins; same rank: larger |score| wins; first one wins a tie
                        if cur is None or rank < cur[0] or (rank == cur[0] and abs(score) > abs(cur[1])):
                            best[k] = (rank, score, motif_id)

                for k, (rank, score, motif_id) in best.items():
                    target = self._target_id(k)
                    if target is None:
                        continue

                    props = {}
                    if self.write_properties:
                        props['effect'] = "gain" if score > 0 else "loss"
                        props['score'] = to_float(score)
                        props['motif_id'] = motif_id
                        props['taxon_id'] = self.taxon_id
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield variant, target, self.label, props
