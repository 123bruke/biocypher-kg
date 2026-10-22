# Adapter and Data Source Catalogue

This catalogue lists every data source integrated into the BioCypher-KG pipeline. Use it to check coverage before adding a new source and to understand what biological domains are represented.

**Source files:**
- [`config/hsa/hsa_data_source_config.yaml`](https://github.com/rejuve-bio/biocypher-kg/blob/main/config/hsa/hsa_data_source_config.yaml) — Human source URLs, version strategies, licenses
- [`config/dmel/dmel_data_source_config.yaml`](https://github.com/rejuve-bio/biocypher-kg/blob/main/config/dmel/dmel_data_source_config.yaml) — Drosophila source URLs
- [`biocypher_metta/adapters/`](https://github.com/rejuve-bio/biocypher-kg/tree/main/biocypher_metta/adapters/) — Adapter implementations

---

## How to add a new data source

See [CONTRIBUTING.md](https://github.com/rejuve-bio/biocypher-kg/blob/main/CONTRIBUTING.md) for the full adapter development guide. In brief:

1. Create a new adapter class in `biocypher_metta/adapters/` (or `hsa/`/`dmel/` for species-specific sources)
2. Implement `get_nodes()` and/or `get_edges()` following the `Adapter` ABC contract
3. Add an entry to the appropriate `*_adapters_config.yaml`
4. Add a URL + version strategy entry to `*_data_source_config.yaml`
5. Run `make run-sample` to validate against sample data

---

> **Species scope note:** Adapters living directly in `biocypher_metta/adapters/` (not in an `hsa/` or `dmel/` subdirectory) are **multi-species** — they are configured per species via the `*_adapters_config.yaml` files and work for any species that includes them. Adapters in `hsa/` or `dmel/` are species-specific implementations. The sections below reflect which species each source is currently configured for, not which species the adapter code supports.

---

## Human (*Homo sapiens*, `hsa`)

### Genomic annotation

| Source | Adapter | Nodes yielded | Edges yielded | URL |
|---|---|---|---|---|
| **GENCODE** | `gencode_gene_adapter.GencodeGeneAdapter` | `gene` | — | [gencodegenes.org](https://www.gencodegenes.org/) |
| **GENCODE** | `gencode_transcript_adapter.GencodeTranscriptAdapter` | `transcript` | `transcribes_to` | same |
| **GENCODE** | `gencode_exon_adapter.GencodeExonAdapter` | `exon` | `part_of_transcript`, `part_of_gene` | same |
| **UniProtKB/Swiss-Prot** | `uniprot_adapter.UniProtAdapter` | `protein` | GO xrefs | [uniprot.org](https://www.uniprot.org/) |
| **UniProtKB/Swiss-Prot** | `uniprot_protein_adapter.UniProtProteinAdapter` | `protein` | protein-xref edges (9 types) | same |

### Protein-protein interactions and coexpression

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **STRING** | `string_ppi_adapter.StringPPIAdapter` | — | `interacts_with` (PPI) | [string-db.org](https://string-db.org/) |
| **STRING** | `string_coexpression_adapter.StringCoexpressionAdapter` | — | `protein_coexpressed_with` | same |
| **CoXpressDB** | `coxpresdb_adapter.CoxpressdbAdapter` | — | `coexpressed_with` | [coxpresdb.jp](https://coxpresdb.jp/) |

### Pathways and reactions

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **Reactome** | `reactome_adapter.ReactomeAdapter` | `pathway`, `reaction` | `reaction_to_pathway` | [reactome.org](https://reactome.org/) |
| **Reactome** | `reactome_edges_adapter.ReactomeEdgesAdapter` | — | gene/protein/small molecule → pathway/reaction | same |
| **Reactome** | `reactome_inference_edges_adapter.ReactomeInferenceEdgesAdapter` | — | inferred pathway relationships | same |
| **Reactome** | `reactome_ppi_adapter.ReactomePPIAdapter` | — | `interacts_with` (physical) | same |
| **Reactome** | `reactome_pathway_go_adapter.ReactomePathwayGoAdapter` | — | `pathway_to_biological_process` etc. | same |

### Gene expression

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **Bgee** | `bgee_adapter.BgeeAdapter` | — | `gene_expressed_in_anatomy`, `gene_expressed_in_developmental_stage` | [bgee.org](https://www.bgee.org/) |
| **GTEx** | `hsa/gtex_expression_adapter.GTExExpressionAdapter` | — | `gene_expressed_in_anatomy` | [gtexportal.org](https://gtexportal.org/) |
| **GTEx eQTL** | `hsa/gtex_eqtl_adapter.GTExEQTLAdapter` | — | variant → gene eQTL edges | same |

### GO annotations

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **GO Annotation File (GAF)** | `gaf_adapter.GAFAdapter` | — | `biological_process_gene`, `molecular_function_gene`, `cellular_component_gene` (and gene_product variants) | [geneontology.org](https://geneontology.org/) |

### Regulatory elements

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **CATLAS** | `hsa/catlas_ccre_adapter.CatlasCCREAdapter` | `regulatory_region` | — | [decoder-genetics.wustl.edu](https://decoder-genetics.wustl.edu/catlasv1/) |
| **CATLAS** | `hsa/catlas_ccre_cell_type_adapter.CatlasCCRECellTypeAdapter` | — | cCRE → cell type | same |
| **CATLAS** | `hsa/catlas_abc_score_adapter.CatlasABCScoreAdapter` | — | `enhancer_gene` (ABC scores by cell type) | same |
| **cCRE (SCREEN)** | `candidate_cis_regulatory_enhancer_adapter.CandidateCisRegulatoryEnhancerAdapter` | `regulatory_region` | `enhancer_gene` | [wenglab.org](https://screen.encodeproject.org/) |
| **cCRE (SCREEN)** | `candidate_cis_regulatory_promoter_adapter.CandidateCisRegulatoryPromoterAdapter` | `regulatory_region` | `promoter_gene` | same |
| **EnhancerAtlas** | `enhancer_atlas_adapter.EnhancerAtlasAdapter` | `regulatory_region` | `enhancer_gene` (tissue-specific) | [singlecelldb.com](http://singlecelldb.com/) |
| **ABC Model** | `hsa/abc_adapter.ABCAdapter` | — | `enhancer_gene` (Activity-By-Contact scores) | [forgedb.cancer.gov](https://forgedb.cancer.gov/) |
| **Roadmap Epigenomics** | `hsa/roadmap_dhs_adapter.RoadmapDHSAdapter` | — | DHS → gene | Broad Institute |
| **Roadmap Epigenomics** | `hsa/roadmap_h3_marks_adapter.RoadmapH3MarksAdapter` | `epigenomic_feature` | — | same |
| **Roadmap Epigenomics** | `hsa/roadmap_state_adapter.RoadmapStateAdapter` | `epigenomic_feature` | — | same |
| **ENCODE RE2G** | `encode_re2g_adapter.EncodeRe2gAdapter` | — | `regulatory_region_gene` | [encodeproject.org](https://www.encodeproject.org/) |
| **TADMap** | `tadmap_adapter.TADMapAdapter` | `tad` | `in_tad_region` | [cb.csail.mit.edu](https://cb.csail.mit.edu/cb/tadmap/) |
| **EPD** | `epd_adapter.EPDAdapter` | `promoter` | `promoter_gene` | [epd.expasy.org](https://epd.expasy.org/) |
| **dbSUPER** | `hsa/dbsuper_adapter.DBSuperAdapter` | `regulatory_region` | `super_enhancer_gene` | [asntech.org/dbsuper](http://asntech.org/dbsuper/) |
| **TFLink** | `tflink_adapter.TFLinkAdapter` | — | `tf_gene` (TF → gene regulatory) | [tflink.net](http://www.tflink.net/) |
| **TFBS/HOCOMOCO** | `tfbs_adapter.TfbsAdapter` | `tfbs` | `gene_tfbs` | [hgdownload.soe.ucsc.edu](https://hgdownload.soe.ucsc.edu/) |
| **HOCOMOCO motifs** | `hocomoco_motif_adapter.HOCOMOCOMotifAdapter` | `motif` | — | [hocomoco11.autosome.org](https://hocomoco11.autosome.org/) |
| **RefSeq closest gene** | `hsa/refseq_closest_gene_adapter.RefSeqClosestGeneAdapter` | — | regulatory element → gene | [forgedb.cancer.gov](https://forgedb.cancer.gov/) |

### Genetic variation

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **dbSNP** | `hsa/dbsnp_adapter.DBSNPAdapter` | `sequence_variant` | variant → gene overlap | [ncbi.nlm.nih.gov/snp](https://www.ncbi.nlm.nih.gov/snp/) |
| **dbVAR** | `hsa/dbvar_adapter.DBVARAdapter` | `structural_variant` | `structural_variant_overlaps_*` | [ncbi.nlm.nih.gov/dbvar](https://www.ncbi.nlm.nih.gov/dbvar/) |
| **DGV** | `hsa/dgv_variant_adapter.DGVVariantAdapter` | `structural_variant` | `structural_variant_overlaps_*` | [dgv.tcag.ca](http://dgv.tcag.ca/) |
| **CADD** | `hsa/cadd_adapter.CADDAdapter` | — | variant pathogenicity scores | [forgedb.cancer.gov](https://forgedb.cancer.gov/) |
| **FAVOR** | `hsa/favor_adapter.FAVORAdapter` | `sequence_variant` | functional annotations | [favor.genohub.org](https://favor.genohub.org/) |
| **PolyPhen-2** | `hsa/polyphen2_adapter.PolyPhen2Adapter` | — | protein mutation effects | [genetics.bwh.harvard.edu](https://genetics.bwh.harvard.edu/pph2/) |
| **GWAS Catalog** | `hsa/gwas_adapter.GWASAdapter` | — | `gene_disease` (via GWAS associations) | [ebi.ac.uk/gwas](https://www.ebi.ac.uk/gwas/) |
| **TopLD** | `hsa/topld_adapter.TopLDAdapter` | — | LD blocks (4 populations: AFR, EAS, EUR, SAS) | [topld.genetics.unc.edu](http://topld.genetics.unc.edu/) |
| **MotifDiff** | `motif_diff_adapter.MotifDiffAdapter` | — | `tf_snp` (motif difference scores) | [github.com/rezwanhosseini/MotifDiff](https://github.com/rezwanhosseini/MotifDiff) |
| **PEREGRINE** | `hsa/peregrine_adapter.PeregrineAdapter` | — | long-range interactions | — |

### Disease and phenotype

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **HPO** | `hsa/human_phenotype_ontology_adapter.HumanPhenotypeOntologyAdapter` | `phenotype` | `hpo_subclass_of` | [hpo.jax.org](https://hpo.jax.org/) |
| **HPO** | `hsa/hpo_gene_disease_adapter.HPOGeneDiseaseAdapter` | — | `gene_disease`, `is_implicated_in` | same |
| **HPO** | `hsa/hpo_gene_phenotype_adapter.HPOGenePhenotypeAdapter` | — | `gene_phenotype` | same |
| **Alliance** | `alliance_gene_disease_adapter.AllianceGeneDiseaseAdapter` | — | `is_implicated_in`, `is_model_of`, `is_marker_for` | [alliancegenome.org](https://www.alliancegenome.org/) |
| **Alliance** | `alliance_gene_orthology_adapter.AllianceGeneOrthologyAdapter` | — | `ortholog_of`, `biomarker_via_orthology`, `implicated_via_orthology` | same |
| **RNA Central** | `rna_central_adapter.RNACentralAdapter` | `non_coding_rna` | — | [rnacentral.org](https://rnacentral.org/) |

### Ontologies (human)

| Ontology | Adapter | Nodes | Edges | CURIE prefix |
|---|---|---|---|---|
| **Gene Ontology** | `gene_ontology_adapter.GeneOntologyAdapter` | `biological_process`, `molecular_function`, `cellular_component` | `*_subclass_of` | `GO:` |
| **Disease Ontology** | `disease_ontology_adapter.DiseaseOntologyAdapter` | `disease` | `do_subclass_of` | `DOID:` |
| **ChEBI** | `chebi_ontology_adapter.ChEBIOntologyAdapter` | `small_molecule` | `chebi_subclass_of` | `CHEBI:` |
| **Cell Ontology** | `cell_ontology_adapter.CellOntologyAdapter` | `cell_type` | `cl_subclass_of`, `cl_capable_of`, `cl_part_of` | `CL:` |
| **Cell Line Ontology** | `cell_line_ontology_adapter.CellLineOntologyAdapter` | `cell_line` | `clo_subclass_of`, `cell_line_is_a_cell_type` | `CLO:` |
| **UBERON** | `uberon_adapter.UberonAdapter` | `anatomy` | `uberon_subclass_of` | `UBERON:` |
| **HsapDv** | `hsapdv_ontology_adapter.HsapDvOntologyAdapter` | `developmental_stage` | hierarchy edges | `HsapDv:` |
| **OMIM** | `omim_ontology_adapter.OmimOntologyAdapter` | `disease` | disease hierarchy | `OMIM:` |
| **BRENDA BTO** | `brenda_tissue_ontology_adapter.BrendaTissueOntologyAdapter` | `tissue` | `bto_subclass_of` | `BTO:` |
| **EFO** | `experimental_factor_ontology_adapter.ExperimentalFactorOntologyAdapter` | `experimental_factor` | `efo_subclass_of` | `EFO:` |
| **MI/PSI-MI** | `molecular_interactions_ontology_adapter.MolecularInteractionsOntologyAdapter` | `molecular_interaction` | `mi_subclass_of` | `MI:` |
| **Sequence Ontology** | `sequence_ontology_adapter.SequenceOntologyAdapter` | `sequence_type` | `so_subclass_of` | `SO:` |
| **Generic OBO** | `ontologies_adapter.OntologiesAdapter` | `ontology_term` | `subclass_of` | varies |

> **Note:** `DiseaseOntologyAdapter` requires `BIOPORTAL_API_KEY` environment variable to download the Disease Ontology from BioPortal.

---

## Drosophila (*Drosophila melanogaster*, `dmel`)

### Genomic annotation

| Source | Adapter | Nodes | Edges | URL |
|---|---|---|---|---|
| **GENCODE/Ensembl** | `gencode_gene_adapter.GencodeGeneAdapter` | `gene` | — | [ensembl.org](https://www.ensembl.org/) |
| **GENCODE/Ensembl** | `gencode_transcript_adapter.GencodeTranscriptAdapter` | `transcript` | `transcribes_to` | same |
| **GENCODE/Ensembl** | `gencode_exon_adapter.GencodeExonAdapter` | `exon` | `part_of_transcript`, `part_of_gene` | same |
| **UniProtKB/Swiss-Prot** | `uniprot_adapter.UniProtAdapter` | `protein` | GO xrefs | [uniprot.org](https://www.uniprot.org/) |
| **RNA Central** | `rna_central_adapter.RNACentralAdapter` | `non_coding_rna` | — | [rnacentral.org](https://rnacentral.org/) |

### FlyBase (comprehensive Drosophila database)

All FlyBase adapters source data from [flybase.org](https://flybase.org/).

| FlyBase data | Adapter | Nodes | Edges |
|---|---|---|---|
| Alleles | `dmel/allele_adapter.AlleleAdapter` | `allele` | — |
| Allele genetic interactions | `dmel/allele_genetic_interaction_adapter.AlleleGeneticInteractionAdapter` | — | genetic interaction edges |
| Gene groups | `dmel/gene_group_adapter.GeneGroupAdapter` | `gene_group` | gene group membership |
| Gene sequence ontology | `dmel/gene_so_adapter.GeneSoAdapter` | — | gene → sequence_type |
| Genetic associations | `dmel/gene_genetic_association_adapter.GeneGeneticAssociationAdapter` | — | genetic association edges |
| Genotype-phenotype sets | `dmel/genotype_phenotype_set_adapter.GenotypePhenotypeSetAdapter` | `genotype`, `phenotype_set` | genotype → phenotype |
| Disease models | `dmel/disease_model_adapter.DiseaseModelAdapter` | `dmel_disease_model` | gene → disease |
| Orthology | `dmel/orthology_adapter.OrthologyAdapter` | — | `ortholog_of` (dmel↔hsa) |
| Paralogy | `dmel/paralogy_adapter.ParalogyAdapter` | — | `paralog_of` |
| Physical interactions (PSI-MI) | `dmel/dmel_physical_interaction_psimi_adapter.DmelPhysicalInteractionPSIMIAdapter` | — | `interacts_with` |
| RNASeq libraries | `dmel/RNASeq_library_adapter.RNASeqLibraryAdapter` | `rnaseq_library` | expression library metadata |
| Expression (expressed_in) | `dmel/expressed_in_adapter.ExpressedInAdapter` | — | `expressed_in` (anatomy/developmental stage) |
| Expression values | `dmel/expression_value_adapter.ExpressionValueAdapter` | — | quantitative expression edges |
| GO annotations | `gaf_adapter.GAFAdapter` (FlyBase GAF) | — | GO annotation edges |

### Shared sources used for Drosophila

| Source | Adapter | Notes |
|---|---|---|
| **Bgee** | `bgee_adapter.BgeeAdapter` | Drosophila-specific expression file |
| **Reactome** | `reactome_*.py` adapters | Dmel-specific pathway entries |
| **STRING** | `string_ppi_adapter.StringPPIAdapter` | taxon 7227 |
| **STRING** | `string_coexpression_adapter.StringCoexpressionAdapter` | taxon 7227 |
| **TFLink** | `tflink_adapter.TFLinkAdapter` | Drosophila interactions |
| **EPD** | `epd_adapter.EPDAdapter` | Drosophila promoters |
| **Alliance** | `alliance_gene_disease_adapter`, `alliance_gene_orthology_adapter` | Multi-species |
| **MotifDiff** | `motif_diff_adapter.MotifDiffAdapter` | SNP to TF gene binding differences (JASPAR motifs, taxon 7227) |

### Drosophila-specific ontologies

| Ontology | Adapter | Nodes | CURIE prefix |
|---|---|---|---|
| **FBcv** (FlyBase Controlled Vocabulary) | `dmel/FBcontrolled_vocabulary_ontology_adapter.FBcontrolledVocabularyOntologyAdapter` | `ontology_term` | `FBcv:` |
| **FBdv** (FlyBase Developmental Stages) | `dmel/FBdevelopment_ontology_adapter.FBdevelopmentOntologyAdapter` | `developmental_stage` | `FBdv:` |
| **FBbt** (FlyBase Anatomy) | `dmel/FBgross_anatomy_ontology_adapter.FBgrossAnatomyOntologyAdapter` | `anatomy` | `FBbt:` |

---

## Cross-species sources

These sources are used for both `hsa` and `dmel` (and potentially other species):

| Source | Description |
|---|---|
| **Alliance of Genome Resources** | Gene-disease, orthology (multi-species) |
| **Reactome** | Pathways and reactions (human + dmel releases) |
| **Gene Ontology** | GO terms and annotations via GAF files |
| **STRING** | PPI and coexpression (species-specific protein links files) |
| **UniProt** | Protein sequences and cross-references |
| **MotifDiff** | SNP-to-transcription factor binding affinity differences (HOCOMOCO for human, JASPAR for Drosophila) |

---

## Coverage matrix

| Domain | hsa | dmel | mmu | cel | rno |
|---|---|---|---|---|---|
| Gene / transcript / exon | Yes | Yes | - | - | - |
| Protein | Yes | Yes | - | - | - |
| Protein-protein interactions | Yes | Yes | - | - | - |
| Pathways / reactions | Yes | Yes | - | - | - |
| GO annotations | Yes | Yes | - | - | - |
| Gene expression | Yes | Yes | - | - | - |
| Regulatory elements | Yes | Yes (MotifDiff) | - | - | - |
| Genetic variants | Yes | Yes (alleles) | - | - | - |
| Disease / phenotype | Yes | Yes | - | - | - |
| Orthology / paralogy | Yes | Yes | - | - | - |
| Ontologies (DO, HPO, ChEBI...) | Yes | Partial | - | - | - |
| Non-coding RNA | Yes | Yes | - | - | - |

`-` = in development; `—` = not applicable or not yet implemented.

---

## Adapter count summary

| Namespace | Count |
|---|---|
| Shared adapters (`biocypher_metta/adapters/`) | ~40 (including ontology adapters) |
| Human-specific (`biocypher_metta/adapters/hsa/`) | 24 |
| Drosophila-specific (`biocypher_metta/adapters/dmel/`) | 16 adapters + 2 utilities |
| **Total adapter classes** | **~80** |

Verification: `find biocypher_metta/adapters -name "*.py" | grep -v __pycache__ | wc -l` (subtract `__init__.py`, `helpers.py`, `reactome_constants.py`, `flybase_tsv_reader.py`).
