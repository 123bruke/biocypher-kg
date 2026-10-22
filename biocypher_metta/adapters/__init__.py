# Author Abdulrahman S. Omar <xabush@singularitynet.io>


class Adapter:
    CURIE_PREFIX = {
        6239: 'WormBase',
        7227: 'FlyBase',
        9606: 'ENSEMBL',
        10090: 'ENSEMBL',
        10116: 'ENSEMBL',
    }
    
    # Mapping of organism taxon_ids to their NCBI category and gene_info filename
    SPECIES_INFO = {
        6239: {
            'category': 'Invertebrates',
            'filename': 'Caenorhabditis_elegans.gene_info.gz',
            'full_name': 'Caenorhabditis elegans',
            'dbxref_prefix': 'WormBase',
            'reactome_prefix': 'R-CEL',
            'ncbi_gene_info_url': "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz",
            'features_data_url': "https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-63/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.63.gtf.gz",
            'entrez_ensembl_cache_directory': 'aux_files/cel/entrez_ensembl',
            'ensembl_uniprot_cache_directory': 'aux_files/cel/ensembl_uniprot',
            'ensembl_uniprot_organism': 'CAEEL_6239',
            'update_interval_hours': 168,
            'eva_source_url': (
                "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_9"
                "/by_species/caenorhabditis_elegans/WBcel235/"
            ),
        },
        7227: {
            'category': 'Invertebrates',
            'filename': 'Drosophila_melanogaster.gene_info.gz',
            'full_name': 'Drosophila melanogaster',
            'dbxref_prefix': 'FLYBASE',
            'reactome_prefix': 'R-DME',
            'ncbi_gene_info_url': "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Drosophila_melanogaster.gene_info.gz",
            'features_data_url': "https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/gtf/dmel-all-r6.67.gtf.gz",
            'entrez_ensembl_cache_directory': 'aux_files/dmel/entrez_ensembl',
            'ensembl_uniprot_cache_directory': 'aux_files/dmel/ensembl_uniprot',
            'ensembl_uniprot_organism': 'DROME_7227',
            'update_interval_hours': 4320,
            'eva_source_url': (
                "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_9"
                "/by_species/drosophila_melanogaster/Release_6_plus_ISO1_MT/"
            ),
        },
        9606: {
            'category': 'Mammalia',
            'filename': 'Homo_sapiens.gene_info.gz',
            'full_name': 'Homo sapiens',
            'dbxref_prefix': 'Ensembl',
            'reactome_prefix': 'R-HSA',
            'ncbi_gene_info_url': "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
            'features_data_url': "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
            'entrez_ensembl_cache_directory': 'aux_files/hsa/entrez_ensembl',
            'ensembl_uniprot_cache_directory': 'aux_files/hsa/ensembl_uniprot',
            'ensembl_uniprot_organism': 'HUMAN_9606',
            'update_interval_hours': 168
        },
        10090: {
            'category': 'Mammalia',
            'filename': 'Mus_musculus.gene_info.gz',
            'full_name': 'Mus musculus',
            'dbxref_prefix': 'Ensembl',
            'reactome_prefix': 'R-MMU',
            'ncbi_gene_info_url': "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
            'features_data_url': "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz",
            'entrez_ensembl_cache_directory': 'aux_files/mmu/entrez_ensembl',
            'ensembl_uniprot_cache_directory': 'aux_files/mmu/ensembl_uniprot',
            'ensembl_uniprot_organism': 'MOUSE_10090',
            'update_interval_hours': 168,
            'eva_source_url': (
                "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_9"
                "/by_species/mus_musculus/GRCm39/"
            ),
        },
        10116: {
            'category': 'Mammalia',
            'filename': 'Rattus_norvegicus.gene_info.gz',
            'full_name': 'Rattus norvegicus',
            'dbxref_prefix': 'RGD',
            'reactome_prefix': 'R-RNO',
            'ncbi_gene_info_url': "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz",
            'features_data_url': "https://ftp.ensembl.org/pub/release-116/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.116.gtf.gz",
            'entrez_ensembl_cache_directory': 'aux_files/rno/entrez_ensembl',
            'ensembl_uniprot_cache_directory': 'aux_files/rno/ensembl_uniprot',
            'ensembl_uniprot_organism': 'RAT_10116',
            'update_interval_hours': 168,
            'eva_source_url': (
                "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_9"
                "/by_species/rattus_norvegicus/GRCr8/"
            ),
        },
        # Add more organisms as needed
    }

    def __init__(self, write_properties, add_provenance):
        self.write_properties = write_properties
        self.add_provenance = add_provenance

    def get_nodes(self):
        pass

    def get_edges(self):
        pass