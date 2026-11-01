import csv
import gzip
import pickle
from biocypher_metta.adapters import Adapter
from biocypher_metta.processors.hgnc_processor import HGNCProcessor
from biocypher_metta.processors.gene_symbol_ensembl_processor import (
    GeneSymbolEnsemblProcessor,
    GENE_SYMBOL_CONFIGS,
)
try:
    from biocypher._logger import logger
except ImportError:
    import logging
    logger = logging.getLogger("biocypher")

# Human data:
# https://epd.expasy.org/ftp/epdnew/H_sapiens/

# Example EPD bed input file:
##CHRM Start  End   Id  Score Strand -  -
# chr1 959245 959305 NOC2L_1 900 - 959245 959256
# chr1 960583 960643 KLHL17_1 900 + 960632 960643
# chr1 966432 966492 PLEKHN1_1 900 + 966481 966492
# chr1 976670 976730 PERM1_1 900 - 976670 976681

# Fly data:   https://epd.expasy.org/ftp/epdnew/D_melanogaster/
# CEL data:   https://epd.expasy.org/ftp/epdnew/C_elegans/current/
# Mouse data: https://epd.expasy.org/ftp/epdnew/M_musculus/
# Rat data:   https://epd.expasy.org/ftp/epdnew/R_norvegicus/

_SOURCE_URLS = {
    9606:  'https://epd.expasy.org/ftp/epdnew/H_sapiens/',
    7227:  'https://epd.expasy.org/ftp/epdnew/D_melanogaster/',
    6239:  'https://epd.expasy.org/ftp/epdnew/C_elegans/',
    10090: 'https://epd.expasy.org/ftp/epdnew/M_musculus/',
    10116: 'https://epd.expasy.org/ftp/epdnew/R_norvegicus/',
}

class EPDAdapter(Adapter):
    INDEX = {'chr': 0, 'coord_start': 1, 'coord_end': 2, 'gene_id': 3}

    def __init__(
        self,
        filepath,
        label,
        hgnc_to_ensembl_map=None,
        write_properties=None,
        add_provenance=None,
        taxon_id=9606,
        type='promoter',
        delimiter=' ',
        chr=None,
        start=None,
        end=None,
        hgnc_processor=None,
        gene_symbol_processor=None,
        gene_symbol_cache_base_dir='aux_files',
        assembly=None,
    ):
        """
        :param hgnc_to_ensembl_map:          DEPRECATED — legacy pickle path for non-human
                                              species (still supported for dmel).
        :param gene_symbol_processor:         Pre-built GeneSymbolEnsemblProcessor instance.
                                              If None, one is created automatically from
                                              GENE_SYMBOL_CONFIGS for supported species.
        :param gene_symbol_cache_base_dir:    Base dir for auto-created processor caches
                                              (e.g. 'aux_files' → 'aux_files/mmu/gene_symbol_ensembl').
        :param assembly:                      Genome assembly (e.g. 'GRCh38', 'dm6', 'ce6', 'mm10', 'rn6').
        """
        self.filepath = filepath
        self.taxon_id = taxon_id
        self.type = type
        self.label = label
        self.delimiter = delimiter
        self.chr = chr
        self.start = start
        self.end = end
        self.assembly = assembly
        self.source = 'EPD'
        self.version = '006'
        self.source_url = _SOURCE_URLS.get(taxon_id, 'https://epd.expasy.org/')

        # --- ID resolver setup ---
        if taxon_id == 9606:
            # Human: always use HGNCProcessor
            if hgnc_processor is not None:
                self.hgnc_processor = hgnc_processor
            else:
                self.hgnc_processor = HGNCProcessor()
                self.hgnc_processor.load_or_update()
            self._symbol_map = None
            self._gene_symbol_processor = None

        elif gene_symbol_processor is not None:
            # Caller provided a ready-made processor
            self.hgnc_processor = None
            self._gene_symbol_processor = gene_symbol_processor
            self._gene_symbol_processor.load_or_update()
            self._symbol_map = None

        elif hgnc_to_ensembl_map is not None:
            # Legacy pickle (used by existing dmel config)
            self.hgnc_processor = None
            self._gene_symbol_processor = None
            self._symbol_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))

        elif taxon_id in GENE_SYMBOL_CONFIGS:
            # Auto-create processor from registry
            self.hgnc_processor = None
            self._gene_symbol_processor = GeneSymbolEnsemblProcessor.for_taxon(
                taxon_id, cache_base_dir=gene_symbol_cache_base_dir
            )
            self._gene_symbol_processor.load_or_update()
            self._symbol_map = None

        else:
            logger.warning(
                f"EPDAdapter: no symbol resolver configured for taxon_id={taxon_id}. "
                f"All edges will be skipped."
            )
            self.hgnc_processor = None
            self._gene_symbol_processor = None
            self._symbol_map = None

        super(EPDAdapter, self).__init__(write_properties, add_provenance)

    def _resolve_symbol(self, gene_symbol: str):
        """Return (ensembl_gene_id, prefix) or (None, None) if not found."""
        prefix = Adapter.CURIE_PREFIX.get(self.taxon_id, 'ENSEMBL')

        if self.hgnc_processor is not None:
            ensembl_id = self.hgnc_processor.get_ensembl_id(gene_symbol)
            return ensembl_id, prefix

        if self._gene_symbol_processor is not None:
            ensembl_id = self._gene_symbol_processor.get_ensembl_gene(gene_symbol)
            return ensembl_id, prefix

        if self._symbol_map is not None:
            ensembl_id = self._symbol_map.get(gene_symbol)
            return ensembl_id, prefix

        return None, None

    def get_nodes(self):
        from biocypher_metta.adapters.helpers import build_regulatory_region_id, check_genomic_location

        opener = gzip.open if self.filepath.endswith('.gz') else open
        with opener(self.filepath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                chr = line[EPDAdapter.INDEX['chr']]
                coord_start = int(line[EPDAdapter.INDEX['coord_start']]) + 1
                coord_end = int(line[EPDAdapter.INDEX['coord_end']])
                if self.assembly is not None:
                    promoter_id = f"EPD:{build_regulatory_region_id(chr, coord_start, coord_end, self.assembly)}"
                else:
                    promoter_id = f"EPD:{build_regulatory_region_id(chr, coord_start, coord_end)}"

                if check_genomic_location(self.chr, self.start, self.end, chr, coord_start, coord_end):
                    props = {}
                    if self.write_properties:
                        props['chr'] = chr
                        props['start'] = coord_start
                        props['end'] = coord_end
                        props['taxon_id'] = self.taxon_id
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield promoter_id, self.label, props

    def get_edges(self):
        from biocypher_metta.adapters.helpers import build_regulatory_region_id, check_genomic_location

        opener = gzip.open if self.filepath.endswith('.gz') else open
        with opener(self.filepath, 'rt') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            for line in reader:
                chr = line[EPDAdapter.INDEX['chr']]
                coord_start = int(line[EPDAdapter.INDEX['coord_start']]) + 1
                coord_end = int(line[EPDAdapter.INDEX['coord_end']])
                gene_symbol = line[EPDAdapter.INDEX['gene_id']].split('_')[0]

                ensembl_id, prefix = self._resolve_symbol(gene_symbol)
                if ensembl_id is None:
                    continue

                ensembl_gene_id = f"{prefix}:{ensembl_id}"

                if check_genomic_location(self.chr, self.start, self.end, chr, coord_start, coord_end):
                    if self.assembly is not None:
                        promoter_id = f"EPD:{build_regulatory_region_id(chr, coord_start, coord_end, self.assembly)}"
                    else:
                        promoter_id = f"EPD:{build_regulatory_region_id(chr, coord_start, coord_end)}"
                    props = {}
                    if self.write_properties and self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                    yield promoter_id, ensembl_gene_id, self.label, props
