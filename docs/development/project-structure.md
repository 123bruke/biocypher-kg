# Project Structure

This document describes the repository layout and the purpose of every significant directory and file. For a higher-level view of how components interact, see [system-overview.md](../architecture/system-overview.md).

---

## Root level

```
biocypher-kg/
├── create_knowledge_graph.py   # Main pipeline entry point (Typer CLI, ~1,670 lines)
├── checkpoint_manager.py       # CheckpointManager — resume interrupted runs
├── pyproject.toml              # Package metadata, dependencies (106 direct), pytest config
├── requirements.txt            # Compiled transitive dependency list (uv pip compile output)
├── uv.lock                     # Exact lock file for reproducible installs (uv)
├── Makefile                    # Primary user interface: make run, make test, make neo4j-up ...
├── docker-compose.yml          # 4-service compose: build, import, deploy, mork
├── docker-variables.env        # BC_TABLE_NAME, NEO4J_AUTH, FILL_DB_ON_STARTUP
├── .env.example                # Template for secrets — copy to .env and fill in
├── .python-version             # Pinned Python version for pyenv / uv
├── README.md                   # Quickstart, make commands, Neo4j deployment
├── CONTRIBUTING.md             # Adapter development guide, schema workflow, PR checklist
├── AGENTS.md                   # Documentation architecture instructions for AI agents
└── MORK_README.md              # MORK/Hyperon integration overview
```

---

## `biocypher_metta/` — Core KG generation package

The heart of the pipeline. Contains adapters, writers, and ID-mapping processors.

```
biocypher_metta/
├── __init__.py                 # BaseWriter ABC: write_nodes(), write_edges(), finalize()
├── metta_writer.py             # MeTTaWriter — .metta files for OpenCog Hyperon
├── neo4j_writer.py             # Neo4jWriter — direct Bolt driver writes
├── neo4j_csv_writer.py         # Neo4jCSVWriter — TSV + Cypher for neo4j-admin import
├── prolog_writer.py            # PrologWriter — .pl logic files
├── parquet_writer.py           # ParquetWriter — Apache Parquet columnar files
├── kgx_writer.py               # KGXWriter — KGX JSON format
├── networkx_writer.py          # NetworkXWriter — in-memory NetworkX graph
├── provenance.py               # ProvenanceLookup — reads download_manifest.json
│
├── adapters/                   # 80 adapter classes across three namespaces
│   ├── __init__.py             # Adapter ABC: get_nodes() → (id, label, props); get_edges() → (src, tgt, label, props)
│   ├── helpers.py              # Shared utilities (type conversions, data cleaning)
│   ├── reactome_constants.py   # Reactome-specific constants
│   │
│   ├── [40 shared adapters]    # Cross-species adapters (human + dmel both use these)
│   │   ├── gene_ontology_adapter.py
│   │   ├── disease_ontology_adapter.py
│   │   ├── uniprot_adapter.py / uniprot_protein_adapter.py
│   │   ├── gencode_gene/transcript/exon adapters (3)
│   │   ├── string_ppi_adapter.py / string_coexpression_adapter.py
│   │   ├── reactome_*.py (5 adapters)
│   │   ├── bgee_adapter.py / coxpresdb_adapter.py
│   │   ├── [ontology adapters: chebi, cell, cell_line, uberon, hsapdv, omim, brenda, efo, mi, seq, generic obo, gaf]
│   │   ├── [regulatory: cCRE enhancer/promoter, EnhancerAtlas, TFBS, HOCOMOCO, TFLink, ENCODE RE2G, TADMap, EPD, MotifDiff]
│   │   └── [cross-species: Alliance gene-disease, Alliance orthology, RNAcentral]
│   │
│   ├── hsa/                    # Human-specific adapters (24 files)
│   │   ├── dbsnp_adapter.py    # dbSNP rsID variants
│   │   ├── dbvar_adapter.py    # Structural variants
│   │   ├── dgv_variant_adapter.py
│   │   ├── cadd_adapter.py     # CADD pathogenicity scores
│   │   ├── favor_adapter.py    # FAVOR variant annotations
│   │   ├── polyphen2_adapter.py
│   │   ├── gtex_eqtl_adapter.py / gtex_expression_adapter.py
│   │   ├── gwas_adapter.py
│   │   ├── topld_adapter.py    # Linkage disequilibrium blocks
│   │   ├── abc_adapter.py      # ABC enhancer predictions
│   │   ├── catlas_ccre_adapter.py / catlas_ccre_cell_type_adapter.py / catlas_abc_score_adapter.py
│   │   ├── roadmap_dhs_adapter.py / roadmap_h3_marks_adapter.py / roadmap_state_adapter.py
│   │   ├── hpo_gene_disease_adapter.py / hpo_gene_phenotype_adapter.py / human_phenotype_ontology_adapter.py
│   │   ├── peregrine_adapter.py
│   │   ├── dbsuper_adapter.py  # Super-enhancer database
│   │   └── refseq_closest_gene_adapter.py
│   │
│   └── dmel/                   # Drosophila-specific adapters (16 adapters + 2 utilities)
│       ├── __init__.py
│       ├── flybase_tsv_reader.py         # FlyBase TSV parser utility (not an adapter)
│       ├── allele_adapter.py / allele_genetic_interaction_adapter.py
│       ├── disease_model_adapter.py
│       ├── dmel_physical_interaction_psimi_adapter.py
│       ├── expressed_in_adapter.py / expression_value_adapter.py
│       ├── FBcontrolled_vocabulary_ontology_adapter.py
│       ├── FBdevelopment_ontology_adapter.py / FBgross_anatomy_ontology_adapter.py
│       ├── gene_genetic_association_adapter.py / gene_group_adapter.py / gene_so_adapter.py
│       ├── genotype_phenotype_set_adapter.py
│       ├── orthology_adapter.py / paralogy_adapter.py
│       └── RNASeq_library_adapter.py
│
└── processors/                 # 5 ID-mapping processors
    ├── __init__.py             # Exports DBSNPProcessor
    ├── base_mapping_processor.py  # BaseMappingProcessor ABC
    ├── dbsnp_processor.py      # rsID → genomic position (SQLite cache)
    ├── ensembl_uniprot_processor.py  # ENSP → UniProt AC
    ├── entrez_ensembl_processor.py   # Entrez Gene ID → ENSEMBL
    ├── go_subontology_processor.py   # GO term → subontology (BP/MF/CC)
    ├── hgnc_processor.py       # HGNC ID → gene symbol + xrefs
    └── README.md               # Processor usage, update strategies, cache files
```

---

## `config/` — All configuration files

```
config/
├── biocypher_config.yaml           # BioCypher core: offline mode, biolink model path, Neo4j delimiters
├── biocypher_docker_config.yaml    # Docker-specific variant
├── species_config.yaml             # Species registry: maps hsa/dmel/mmu/cel/rno to config paths
├── primer_schema_config.yaml       # Shared base schema: 36 node types, 108 edge types
├── biolink-model.yaml              # Biolink semantic model specification
├── biolink-model.owl.ttl           # OWL/TTL representation of Biolink
├── yaml_loader.py                  # YAML loader with !include directive support
│
├── hsa/                            # Human configuration
│   ├── hsa_schema_config.yaml      # Human-specific node/edge types (extends primer schema)
│   ├── hsa_adapters_config.yaml    # Full production adapter registry
│   ├── hsa_adapters_config_sample.yaml  # Sample adapter registry (used in CI)
│   └── hsa_data_source_config.yaml # 40+ data sources: URL, version strategy, license
│
├── dmel/                           # Drosophila configuration
│   ├── dmel_schema_config.yaml     # Drosophila-specific types (allele, genotype, rnaseq_library...)
│   ├── dmel_adapters_config.yaml
│   ├── dmel_adapters_config_sample.yaml
│   ├── dmel_data_source_config.yaml
│   └── net_act_adapters_config.yaml  # Network activity adapter config
│
├── cel/                            # C. elegans (data source config only — no adapters yet)
│   └── cel_data_source_config.yaml
├── mmu/                            # Mouse (data source config only — no adapters yet)
│   └── mmu_data_source_config.yaml
└── rno/                            # Rat (data source config only — no adapters yet)
    └── rno_data_source_config.yaml
```

For a complete reference of every configuration field and environment variable, see [configuration.md](../operations/configuration.md).

---

## `biocypher_dataset_downloader/` — Dataset acquisition

```
biocypher_dataset_downloader/
├── __init__.py
├── download_manager.py         # DownloadManager: HTTP/FTP acquisition, version resolution
├── download_data.py            # Download + caching logic
├── protocols/
│   ├── __init__.py
│   ├── base.py                 # Protocol ABC
│   └── http.py                 # HTTP/FTP download implementation
└── versioning/
    ├── __init__.py
    ├── base.py                 # VersionInfo dataclass
    ├── strategies.py           # StaticGetter, UrlRegexGetter, HttpHeadGetter
    ├── registry.py             # Strategy registry
    ├── cli.py                  # check-versions CLI (--species flag)
    └── README.md               # VersionInfo/VersionGetter contract
```

See [ingestion-pipeline.md](../knowledge-graph/ingestion-pipeline.md) for the download → adapter flow.

---

## `biocypher_cli/` — Interactive CLI

```
biocypher_cli/
├── __init__.py
├── cli.py                      # questionary TUI: organism selection, config discovery
├── modules/
│   ├── adapters.py
│   ├── config.py
│   ├── init.py
│   └── utils.py
└── README.md
```

Run with `uv run python biocypher_cli/cli.py`. Supports Human and Drosophila melanogaster, default and custom configuration workflows.

---

## `kg-service/` — REST API and database management

```
kg-service/
├── backend/
│   ├── api/
│   │   ├── main.py             # FastAPI app, CORS, 72-hour cache refresh scheduler
│   │   └── routes/
│   │       ├── __init__.py
│   │       ├── summary.py      # GET /api/summary
│   │       ├── updates.py      # GET /api/updates
│   │       ├── graph_info.py   # GET /api/graph-info
│   │       ├── versions.py     # Version history routes
│   │       └── databases.py    # Database management routes
│   └── core/
│       ├── config.py           # Settings (Pydantic): NEO4J_URI, MORK_URL, ARCHIVE_BASE, etc.
│       ├── neo4j_client.py     # Neo4j driver wrapper
│       └── graph_info_cache.py # 72-hour cached graph_info.json reader
├── neo4j_loader.py             # Neo4jLoader: bulk CSV import + surgical dataset delete
├── version_manager.py          # MD5 hashing, change detection, incremental reload
├── mork_loader.py              # MORK/MeTTa data loader
└── version_diff.py             # Delta computation between KG versions
```

> **Known issue:** `kg-service/backend/api/main.py` line 8 imports `meta` and `entities` from `backend.api.routes`, but neither file exists in the `routes/` directory. The service cannot start until this is resolved. See [endpoints.md](../api/endpoints.md).

---

## `biocypher-mork/` — Rust MORK query server

```
biocypher-mork/
├── Dockerfile                  # Multi-stage Rust build (rustlang/rust:nightly base)
├── docker-compose.yml          # Exposes port 8027 (HOST_PORT env var)
├── entrypoint.sh
├── client.py                   # Python MeTTa query client
├── wal_client.py               # Write-ahead log client
├── data_loader.py              # Loads .metta files into MORK
├── mork_persist/
│   ├── snapshot.paths          # Snapshot file registry
│   └── wal.metta               # Write-ahead log
└── README.md                   # -> MORK_README.md
```

---

## `docker/` — Neo4j deployment

```
docker/
├── docker-compose.neo4j.yml    # Parameterised Neo4j 5.21 deployment
├── neo4j.env.example           # Template — copy to neo4j.env
├── biocypher_entrypoint_patch.sh
├── create_table.sh
└── import.sh
```

---

## `data_source_schemas/` — Auto-generated per-source YAML schemas

Generated by `schema_generator/generate_data_source_schemas.py` on every full pipeline run. **Do not edit manually.** 34 files for `hsa/`, 11 for `dmel/`.

```
data_source_schemas/
├── hsa/      # 34 YAML files: ABC.yaml, STRING.yaml, UniProt.yaml, ...
└── dmel/     # 11 YAML files
```

---

## `aux_files/` — Pre-built ID mapping tables

Pickle files and SQLite databases built by the processor scripts (see `scripts/update_dbsnp.py`, `scripts/dmel_create_entrez_to_ensembl_map.py`). These are read by the five ID-mapping processors at adapter instantiation time.

```
aux_files/
├── hsa/
│   ├── sample_dbsnp/           # ~150 KB bundled dbSNP sample (ships with repo)
│   ├── ensembl_uniprot/        # ENSP → UniProt AC mapping
│   ├── entrez_ensembl/         # Entrez Gene ID → ENSEMBL
│   ├── go_subontology/         # GO term → subontology
│   ├── hgnc/                   # HGNC symbol mappings
│   └── [tissue/regulatory mapping pkl files]
└── dmel/
    ├── dmel_entrez_to_ensembl.pkl
    ├── dmel_string_ensembl_to_uniprot_map.pkl
    └── [FlyBase synonym mappings]
```

---

## `samples/` — Test data

Small representative files used by CI (`test-adapters.yml`) and `make run-sample`. Contains real-format subsets of production data files.

```
samples/
├── hsa/        # 30+ gzipped sample files (GENCODE, UniProt, STRING, GWAS, ...)
├── dmel/       # Drosophila sample files
└── reactome/   # Reactome pathway files
```

---

## `scripts/` — Utility scripts (22 files)

| Script | Purpose |
|---|---|
| `build.sh` | Main Docker build entry point |
| `build_act.sh` / `build_paths.sh` | Build variants |
| `import.sh` | Neo4j import wrapper |
| `neo4j_loader.py` | Direct Neo4j load (no versioning) |
| `mork_loader.py` / `mork_repl.py` / `mork_client.py` | MORK interaction |
| `update_dbsnp.py` | Rebuild dbSNP SQLite cache |
| `dmel_create_entrez_to_ensembl_map.py` | Build Drosophila ID mapping |
| `create_catlas_*.py` (3 scripts) | Build CATLAS mapping pickles |
| `build_mapping.py` | Build mapping configs |
| `convert_topaths.py` | Convert to paths format |
| `load_paths.py` | Path loading utilities |
| `reorganize_schemas.py` / `sort_schemas.py` | Schema maintenance |
| `pkl_cat.py` | Pickle file inspection |
| `count_reactome_protein_roles.py` | Reactome analysis utility |
| `metta_space_import.py` / `das_metta_loader.py` | MeTTa loading utilities |

---

## `test/` — Test suite

```
test/
├── conftest.py         # pytest fixtures: --adapters-config, --primer-schema-config, --adapter-test-mode
├── test.py             # Main adapter/writer/schema validation tests
├── test_provenance.py  # ProvenanceLookup tests
└── test_versioning.py  # VersionInfo / strategy tests
```

Run with `make test` or `uv run pytest`. Two modes: `smoke` (skips heavy ontology adapters) and `full`. See [testing.md](testing.md).

---

## `.github/` — CI/CD

```
.github/
├── workflows/
│   ├── test-adapters.yml           # Smart matrix: test changed adapters + all writers
│   ├── test-schema.yml             # Schema validation: smoke (PR) vs full (main)
│   ├── check-dataset-versions.yml  # Weekly staleness check, creates issue on change
│   ├── ci-quality-analysis.yml     # SonarQube code quality scan
│   └── dependabot-auto-merge.yml   # Auto-merge Dependabot PRs
├── scripts/
│   ├── detect_config_changes.py    # Parse YAML diff → changed adapter blocks
│   ├── detect_heavy_adapter_changes.py  # Identify ontology adapter changes
│   └── prepare_config.py           # Generate filtered test adapter config
└── PULL_REQUEST_TEMPLATE.md
```

---

## `schema_generator/` — Schema auto-generation

```
schema_generator/
├── generate_data_source_schemas.py  # SchemaGenerator class
└── README.md
```

Generates `data_source_schemas/<species>/*.yaml` from raw data source files. Called automatically at the end of every full pipeline run. See `schema_generator/README.md` for usage.

---

## `docs/knowledge-graph/` — Versioning documentation

```
docs/knowledge-graph/
└── dataset-versioning.md    # Versioning strategies, manifest format, provenance in Neo4j
```

The versioning guide was moved from `doc/` into the main `docs/` documentation hierarchy. See [ingestion-pipeline.md](../knowledge-graph/ingestion-pipeline.md) for the end-to-end flow.

---

## `notebooks/` — Jupyter notebooks

```
notebooks/
├── cellxgene_exp_corr.ipynb   # CellxGene expression correlation analysis
├── finemapping_susie.ipynb    # SuSiE fine-mapping analysis
└── llm_connect.ipynb          # LLM integration exploration
```

Exploratory notebooks — not part of the pipeline.

---

## `hgnc_gene_data/` — HGNC gene reference cache

```
hgnc_gene_data/
├── hgnc_data.pkl      # Cached HGNC gene symbol/xref data (built by HGNCProcessor)
└── hgnc_version.txt   # Version record for the cached data
```
