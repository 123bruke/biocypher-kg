'''
# Human:  to be defined…
#
#
# Fly:
# FB https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Alleles_.3C.3D.3E_Genes_.28fbal_to_fbgn_fb_.2A.tsv.29

# FB table columns:
#AlleleID	AlleleSymbol	GeneID	GeneSymbol
FBal0137236	gukh[142]	FBgn0026239	gukh
FBal0137618	Xrp1[142]	FBgn0261113	Xrp1
FBal0092786	Ecol\lacZ[T125]	FBgn0014447	Ecol\lacZ
FBal0100372	Myc[P0]	FBgn0262656	Myc
FBal0009407	kst[01318]	FBgn0004167	kst
FBal0091321	Ecol\lacZ[kst-01318]	FBgn0014447	Ecol\lacZ
FBal0091320	Ecol\lacZ[mam-04615]	FBgn0014447	Ecol\lacZ

# Coordinate convention for 'snp' nodes: chado's featureloc.fmin is 0-based/
# half-open, but EVAAdapter's 'snp' nodes use the VCF's 1-based POS for both
# 'start' and 'end' (start == end, a single point). get_nodes() below converts
# fmin -> fmin + 1 so snp nodes from this adapter and from EVAAdapter share the
# same coordinate system for the same physical base.
'''
from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable
#from flybase_tsv_reader import FlybasePrecomputedTable
from biocypher_metta.adapters import Adapter
import psycopg2


class AlleleAdapter(Adapter):

    def __init__(self, write_properties, add_provenance, dmel_filepath=None, label='allele', taxon_id=7227):
        self.dmel_filepath = dmel_filepath
        self.label = label
        self.source = 'FLYBASE'
        self.source_url = 'https://flybase.org/'
        self.taxon_id = taxon_id
        super(AlleleAdapter, self).__init__(write_properties, add_provenance)
        self.snp_cache = self._load_snp_cache()

    def _load_snp_cache(self):
        """Fetch all SNPs and their locations (including chromosome) from FlyBase in one fast query.

        The chromosome (fs.uniquename, e.g. '2L', '2R', '3L', '3R', '4', 'X') is the
        srcfeature a SNP's featureloc is anchored to — already in the same raw,
        unprefixed arm-name format EVAAdapter uses for its 'chr' property, so no
        further normalization is needed here.
        """
        snp_cache = {}
        conn = None
        try:
            conn = psycopg2.connect(
                host="chado.flybase.org",
                database="flybase",
                user="flybase",
                password="flybase",
                connect_timeout=10
            )
            with conn.cursor() as cursor:
                cursor.execute("""
                    SELECT f.uniquename, fl.fmin, fs.uniquename AS chr
                    FROM feature f
                    LEFT JOIN featureloc fl ON f.feature_id = fl.feature_id
                    LEFT JOIN feature fs ON fl.srcfeature_id = fs.feature_id
                    WHERE f.type_id=733 AND f.is_obsolete=FALSE AND f.is_analysis=FALSE AND f.organism_id=1
                """)
                for uniquename, fmin, chrom in cursor.fetchall():
                    snp_cache[uniquename] = (fmin, chrom)
        except Exception as e:
            print(f"Error connecting to or querying FlyBase: {e}")
        finally:
            if conn is not None:
                conn.close()
        return snp_cache

    def get_nodes(self):
        fbal_table = FlybasePrecomputedTable(self.dmel_filepath)
        self.version = fbal_table.extract_date_string(self.dmel_filepath)
        #header:
        #AlleleID	AlleleSymbol	GeneID	GeneSymbol
        rows = fbal_table.get_rows()

        for row in rows:
            props = {}
            fbal_id = row[0]       # AlleleID e.g. FBal0137236
            allele_symbol = row[1] # AlleleSymbol e.g. gukh[142] — used as uniquename in FlyBase feature table
            allele_id = f'FlyBase:{fbal_id}'
            props['allele_symbol'] = allele_symbol
            props['taxon_id'] = self.taxon_id

            if allele_symbol in self.snp_cache:
                snp_props = props.copy()
                fmin, chrom = self.snp_cache[allele_symbol]
                if fmin is not None:
                    # fmin is 0-based/half-open (chado); convert to the 1-based,
                    # single-point convention EVAAdapter uses (start == end == pos)
                    # so snp nodes from both adapters share one coordinate system.
                    snp_props['start'] = fmin + 1
                    snp_props['end'] = fmin + 1
                if chrom is not None:
                    snp_props['chr'] = chrom
                yield allele_id, 'snp', snp_props
            else:
                yield allele_id, self.label, props      # here label is 'allele'

    def get_edges(self):
        fbal_table = FlybasePrecomputedTable(self.dmel_filepath)
        self.version = fbal_table.extract_date_string(self.dmel_filepath)
        #header:
        #AlleleID	AlleleSymbol	GeneID	GeneSymbol
        rows = fbal_table.get_rows()

        for row in rows:
            props = {}
            fbal_id = row[0]
            allele_symbol = row[1] # AlleleSymbol — used as uniquename in FlyBase feature table
            source = f'FlyBase:{fbal_id}'
            target = f'FlyBase:{row[2]}'
            props['taxon_id'] = self.taxon_id

            is_snp = allele_symbol in self.snp_cache

            if self.label == 'snp_in_gene':
                # Only yield located_in edges for SNPs
                if is_snp:
                    yield source, target, self.label, props
            else:
                # variant_of: mutually exclusive — SNPs get snp_variant_of, others get variant_of
                if is_snp:
                    yield source, target, 'snp_variant_of', props
                else:
                    yield source, target, self.label, props
