# Sample generation algorithm — connected_full_sample_data branch

## Motivation

Generating and loading the **full** knowledge graph for a species takes
hours to days — too slow to answer a simple question: *is the KG actually
one connected graph, or does it fragment into isolated islands?* This branch
exists to answer that in **minutes** instead: a small, per-species sample
that still exercises every configured adapter, built so the resulting graph
is guaranteed to be a single connected component. Load the sample into
Neo4j and get the same connectivity signal the full build would give,
without paying the full build's runtime cost.

Context: `connected_full_sample_data` branch of biocypher-kg. Goal — build a
"sample" file set (per species, under `./samples/<species>/`) small enough to
run every configured adapter, but that produces a **connected** knowledge
graph when loaded into Neo4j (via `neo4j_csv_writer`). First target species:
dmel. The algorithm must be generic enough to reapply to hsa, cel, mmu, rno.

This document is the reference for `scripts/generate_connected_sample.py` and
for adding *new* sample sources (new adapters, new edges on existing
adapters, new node types) without breaking the guarantee this branch exists
to provide: **every species' sample is a single connected component
containing every anchor/closure gene.**

It is referenced (as a relative path comment) from every
`config/<species>/<species>_anchor_genes.yaml` and every
`config/<species>/<species>_adapters_config_sample.yaml` — keep it in sync
with the actual code if the algorithm changes.

## Pipeline findings that motivated the design

- `create_knowledge_graph.py` (mode `--species <sp> --dataset sample
  --writer-type neo4j`) loads `config/species_config.yaml` → adapters/schema
  from `config/<species>/`, merges the schema with
  `config/primer_schema_config.yaml`, runs every adapter active in
  `<species>_adapters_config_sample.yaml`, writes CSV+Cypher per label via
  `Neo4jCSVWriter` (`biocypher_metta/neo4j_csv_writer.py`).
- **Neo4jCSVWriter does not connect to a real Neo4j.** It only generates
  `nodes_*.csv`/`edges_*.csv` + `.cypher` (LOAD CSV + `apoc.periodic.iterate`).
  Edges use `MATCH` (not `MERGE`) to find source/target — if the node doesn't
  exist, the edge is simply not created, even with a correct CSV row.
  **Therefore: every node type referenced as source/target of some edge
  needs an active node adapter**, or the graph ends up with silently
  dangling edges. This single fact is why the whole algorithm exists: a
  sample that's merely "small" isn't enough — every kept edge needs both its
  endpoints to also survive filtering, everywhere, or connectivity breaks
  silently instead of loudly.

## Algorithm: ego-network (k-hop neighborhood) closure

In graph-theory terms, Phases 1-3 below construct the **ego network** of the
anchor gene set — the subgraph reachable within `hops` steps of the anchors —
then rank *which* one-hop neighbors to keep by how many independent sources
connect them back to the existing set (coverage-weighted expansion, not a
uniform BFS frontier). Phases 4-5 materialize that node/edge set into
filtered per-adapter files, tagging any adapter left with zero surviving
real rows for a synthetic fallback instead.

### Inputs
- `species` (dmel, hsa, cel, mmu, rno, ...)
- `anchor_genes`: a manually curated list — housekeeping genes
  (transcription/translation, broadly expressed) **plus a few well-studied
  classic genes** (to cover phenotype/disease/genetic-interaction adapters,
  which tend to be sparse for purely housekeeping genes)
- `size_budget`: approximate node target (soft cap — can go a bit over)
- `hops`: expansion depth (default 1)
- a registry of active adapters + paths to the **complete, real** source
  files (e.g. `/mnt/hdd_1/biocypher-kg/input/dmel/flybase/`)

### Phase 0 — Setup (generic part, doesn't change per species)
Isolate what's species-specific behind a small interface:
- `resolve_anchor_ids(species, symbols) -> gene_ids`
- `load_id_backbone(species) -> {gene→transcripts, gene/transcript→proteins, gene/protein→uniprot}`
- `per_format_filter(adapter) -> function that extracts the key column(s) for
  that specific format` (GTF, FlyBase TSV, GAF, PSI-MI tab, STRING links,
  pkl)

The core of the algorithm (Phases 1-5) doesn't know whether it's FBgn or
ENSG — it only manipulates `id_set` and "source file → filtered rows".

### Phase 1 — Resolve anchors
```
gene_ids = resolve_anchor_ids(species, anchor_genes)
assert each gene_id appears in the species' GTF/GFF   # fail early if an anchor doesn't exist
```

### Phase 2 — Close the ID backbone (gene → transcript → protein → uniprot)
```
id_set = {gene_ids}
id_set |= transcripts_of(gene_ids)      # via fbgn_fbtr_fbpp_expanded
id_set |= proteins_of(id_set)           # same
id_set |= uniprot_of(id_set)            # via fbgn_uniprot / *_ensembl_to_uniprot_map.pkl
```

### Phase 3 — One-hop expansion by coverage (the key connectivity point)
For every real relational source (STRING PPI, coexpression, coxpresdb,
TFLink, GAF, physical interactions, orthology, paralogy, genetic
interaction, genotype-phenotype, disease model...):
```
partner_hits = Counter()   # partner_id -> how many distinct sources connect it to id_set
for source in relational_sources:
    for row in stream(source.real_file):
        a, b = extract_ids(row, source.format)
        if a in id_set and b not in id_set: partner_hits[b] += 1
        if b in id_set and a not in id_set: partner_hits[a] += 1
```
This gives, for every candidate to enter the sample, **how many different
adapters would connect it** — it's not random: it picks whoever maximizes
cross-connections.
```
ranked = sorted(partner_hits, key=hits, reverse=True)
id_set |= top_k(ranked, until size_budget reached)
```

### Phase 4 — Materialize the filtered file per adapter
Second pass, now with the expanded `id_set` — keeps only **closed** edges
(both ends inside the set, otherwise Neo4jCSVWriter's `MATCH` won't match):
```
for source in all_sources:            # includes node sources: GTF, uniprot .dat, etc.
    rows = filter(source.real_file, keep_if=lambda row: all(id in id_set for id in extract_ids(row)))
    write(./samples/<species>/<same relative path>/<same name>, rows)
    manifest[source].real_rows = len(rows)
```

### Phase 5 — Detect gaps and (only if necessary) force a synthetic link
```
for source, count in manifest.items():
    if count == 0:
        log_gap(source)
        synthetic_rows = build_minimal_rows(source.schema, ids=sample(id_set, k=2..3))
        write(./samples/<species>/.../<name>.synthetic.tsv.gz, synthetic_rows)
        manifest[source].synthetic = True
```
Tag it in the filename (`.synthetic.`) + log it in the manifest — never mix
it with the filtered "real" file.

### Phase 6 — Re-enable needed adapters
Re-enable, in `*_adapters_config_sample.yaml`, the node adapters that are
currently commented out but that Phase 4 can now feed (gene/transcript/exon
in dmel's case), pointing at the new files under `./samples/<species>/`.

### Phase 7 — Validate connectivity before going to Neo4j
Run `create_knowledge_graph.py --species <sp> --dataset sample --writer-type
neo4j`, then:
- check `graph_info.json`: every node_type referenced as source/target of
  some edge must have count > 0
- optional/cheap: build a graph with `networkx` from the generated CSVs
  (without needing Neo4j running) and check `nx.is_weakly_connected(G)` —
  fails fast if something is still orphaned

  **Superseded (see §6 below):** for the actual closure-membership check on
  a full-size build, don't use networkx — it caused two machine freezes on
  a 6M-node/45M-edge graph. Use the streaming Union-Find approach instead.

### Generic vs. per-species

| Generic (core) | Per-species (plugin) |
|---|---|
| Phases 1-5 (closure, coverage ranking, cutoff, synthetic tag) | `resolve_anchor_ids`, `load_id_backbone`, `extract_ids` per format |
| Phase 7 (connectivity validation) | list of "relational sources" participating in Phase 3 |
| manifest schema | curated list of anchor genes |

## Candidate anchor_genes (verified against real data)

Source: `/mnt/hdd_1/biocypher-kg/input/multi/alliance/ORTHOLOGY-ALLIANCE_COMBINED_7.tsv.gz`
(Alliance of Genome Resources, multi-algorithm pairwise orthology).
Methodology: for each human gene, required a **reciprocal best-score**
ortholog (best hit in both directions) simultaneously in mmu, rno, dmel and
cel — 5452 human genes passed this filter. Of those, filtered by classic
housekeeping families (translation, ribosome, proteasome, cytoskeleton,
chaperones, DNA replication, splicing, transcription), leaving 187. The
final list below is a curation of ~22 of those, covering distinct
functional categories (avoiding redundancy, e.g. not picking 20 ribosomal
proteins).

| Category | hsa | dmel (FBgn) | cel (WBGene) | mmu | rno |
|---|---|---|---|---|---|
| Cytoskeleton | ACTB | Act5C (FBgn0000042) | act-2 (WBGene00000064) | Actb | Actb |
| Cytoskeleton | TUBA1B | αTub85E (FBgn0003886) | mec-12 (WBGene00003175) | Tuba1b | Tuba1b |
| Glycolysis | GAPDH | Gapdh1 (FBgn0001091) | gpd-1 (WBGene00001683) | Gapdh | Gapdh |
| Glycolysis | ENO1 | Eno (FBgn0000579) | enol-1 (WBGene00011884) | Eno1 | Eno1 |
| Translation | EEF1A1 | eEF1α2 (FBgn0000557) | eef-1A.2 (WBGene00001169) | Eef1a1 | Eef1a1 |
| Translation | EEF2 | eEF2 (FBgn0000559) | eef-2 (WBGene00001167) | Eef2 | Eef2 |
| Ribosome (large) | RPL7 | RpL7 (FBgn0005593) | rpl-7 (WBGene00004418) | Rpl7 | Rpl7 |
| Ribosome (large) | RPL32 | RpL32 (FBgn0002626) | rpl-32 (WBGene00004446) | Rpl32 | Rpl32 |
| Ribosome (small) | RPS6 | RpS6 (FBgn0261592) | rps-6 (WBGene00004475) | Rps6 | Rps6 |
| Ribosome (small) | RPS2 | RpS2 (FBgn0004867) | rps-2 (WBGene00004471) | Rps2 | Rps2 |
| Proteasome α | PSMA1 | Prosα6T (FBgn0032492) | pas-6 (WBGene00003927) | Psma1 | Psma1 |
| Proteasome β | PSMB5 | Prosβ5 (FBgn0029134) | pbs-5 (WBGene00003951) | Psmb5 | Psmb5 |
| Proteasome ATPase | PSMC1 | Rpt2 (FBgn0015282) | rpt-2 (WBGene00004502) | Psmc1 | Psmc1 |
| Ubiquitin | UBB | CG11700 (FBgn0029856) | F52C6.3 (WBGene00018660) | Ubb | Ubb |
| Ubiquitin | UBA52 | RpL40 (FBgn0003941) | ubq-2 (WBGene00006728) | Uba52 | Uba52 |
| Chaperone | HSPA8 | Hsc70-4 (FBgn0266599) | hsp-1 (WBGene00002005) | Hspa8 | Hspa8 |
| Chaperone | CCT2 | CCT2 (FBgn0030086) | cct-2 (WBGene00000378) | Cct2 | Cct2 |
| Histone | H2BC1 | His2B:CG33868 (FBgn0053868) | his-34 (WBGene00001908) | H2bc1 | H2bc1 |
| DNA replication | PCNA | PCNA (FBgn0005655) | pcn-1 (WBGene00003955) | Pcna | Pcna |
| DNA replication | MCM2 | Mcm2 (FBgn0014861) | mcm-2 (WBGene00003154) | Mcm2 | Mcm2 |
| Transcription (Pol II) | POLR2A | Polr2A (FBgn0003277) | ama-1 (WBGene00000123) | Polr2a | Polr2a |
| Splicing | SF3B1 | Sf3b1 (FBgn0031266) | sftb-1 (WBGene00011605) | Sf3b1 | Sf3b1 |
| Antioxidant | SOD1 | Sod1 (FBgn0003462) | sod-1 (WBGene00004930) | Sod1 | Sod1 |

Caveats:
- Several human paralogs collapse onto the same gene in dmel/cel (e.g.
  `H4C1`...`H4C16` all point to a single `His4r` in the fly; `RPS6KA1/2/3/6`
  all point to `S6kII`) — real species annotation (fewer gene copies there),
  not an error. Already avoided in the list above.
- The confidence score (number of agreeing algorithms) ties at ~38 for
  practically all 187 — not a fine-grained ranking, more of a binary filter
  "reliable ortholog across all 5" vs. not.
- Phenotype/disease/genetic-interaction coverage wasn't tested specifically
  for these genes yet — may need to supplement with 2-3 well-studied
  classic genes (e.g. `white`, `hh`, `per` in dmel) in Phase 3, since pure
  housekeeping genes tend to have little disease-model / genotype-phenotype
  annotation.

## Decisions (all resolved)

- **Cross-validation**: 100% coverage in STRING/TFLink/GAF for the 23
  anchors; 22/23 in coxpresdb (only His2B:CG33868, a histone, has no file —
  expected).
- **size_budget**: 180 nodes (soft cap).
- **`fbgn_uniprot` missing from `dmel_data_source_config.yaml`**: fixed, URL
  added to the `flybase:` block.
- **Synthetic tag**: `.SYNTHETIC.` (uppercase, stands out well) inserted
  before the real extension chain, e.g.
  `disease_model_annotations.SYNTHETIC.tsv.gz`. Manifest:
  `samples/<species>/sample_manifest.json` with `real_rows`/`synthetic`/
  `synthetic_rows`/`reason` per source.

## dmel v1 implemented and validated (2026-07-18/19)

`scripts/generate_connected_sample.py` + `config/dmel/dmel_anchor_genes.yaml`
+ `config/dmel/dmel_adapters_config_sample.yaml` (rewritten, core + FlyBase
extension). `scripts/check_kg_connectivity.py` created to validate.

Result: 180 genes (23 anchors + 157 by coverage expansion), 0 sources needed
a synthetic fallback. Final KG: **679,241 nodes, 883,582 edges**, 817,317
resolved (94%), **180 genes 100% in the same connected component** (giant
component of 110,761 nodes).

Two pre-existing bugs found and fixed (also affect production):
- `biocypher_metta/adapters/dmel/expression_value_adapter.py:98` — missing
  `:` in `f'FlyBase{row[11]...}'`, zeroed out the 513,662
  `expression_value` edges.
- `biocypher_metta/adapters/dmel/expressed_in_adapter.py:100` — spurious
  `FlyBase:` prefix on the target (already normalized), zeroed out the
  61,151 `expressed_in` edges.
Also: the `aux_files/dmel/entrez_ensembl/` cache was empty (same WAF
bot-challenge on s3ftp, hitting `EntrezEnsemblProcessor` via
`Adapter.SPECIES_INFO[7227]['features_data_url']`) — fixed with
`scripts/build_dmel_entrez_ensembl_cache.py` (local bootstrap, no network).

`scripts/generate_connected_sample.py` also gained: file resolution by
**glob** (`resolve_source_file`) instead of exact name, to survive release
bumps (`_fb_2026_01` → `_fb_2026_02` already tested on another machine),
output always under a **stable, version-free** name, and a directory
fallback (`AUX_FALLBACK_ROOTS`) for files that `download_data.py` moves
outside `input_dir` (`fbal_to_fbgn` → `aux_files/dmel/`). `download_data.py`
now calls `generate_connected_sample.py` automatically after a complete
download (not on a partial `--source`), saving to `samples/<species>/`.

## Generalization to the other 4 species (2026-07-20, in progress)

Structural finding: the core (gencode gene/transcript/exon, uniprot, STRING
PPI+coexpression, coxpresdb, TFLink, GAF, GO) uses the **same generic
modules** across all 5 species. dmel's FlyBase-specific backbone tables
(`fbgn_fbtr_fbpp_expanded`, `fbgn_uniprot`) and its extra adapters
(gene_group, disease_model, genotype_phenotype, allele,
physical_interaction, gene_genetic_interaction) don't exist for the other 4
— refactored into `build_generic_backbone_maps()`: gene→transcript via the
GTF's own "transcript" rows, gene→protein→uniprot via `DR   Ensembl;
ENST...; ENSP...; ENSG...;` lines in `uniprot.dat`, gene→entrez via
`EntrezEnsemblProcessor` (the same cache the real adapters use).
`SOURCE_FILES`/`SPECIES` became per-species dicts. `build_global_maps()`
(FlyBase) stayed reserved for dmel only.

GAF per species does **not** always key on the gene ID — hsa uses UniProt
accession, mmu uses MGI (dbXrefs "MGI:MGI:x" → "MGI:x"), rno uses RGD
**without a prefix** (dbXrefs "RGD:x" → bare "x"), cel/dmel use the native
gene ID directly. Parameterized via `gaf_id_space` per species.

dmel post-refactor regression: identical numbers (679,241/883,582, 94%, 1
component) — confirmed, zero regression.

Stable IDs for the 23 anchor genes resolved across the 4 species via
`ORTHOLOGY-ALLIANCE_COMBINED_7.tsv.gz` (human→mmu/rno/cel, reciprocal-best) +
`hgnc_mapping.pkl` (hsa: HGNC→ENSG) + local mmu/rno `gene_info.gz`
(symbol→Ensembl dbXrefs) — `config/{hsa,mmu,rno,cel}_anchor_genes.yaml`
written, same pattern as dmel.

Real data confirmed available in all 4 (checked with `find`/`du`, not just
`ls`): gencode/ensembl GTF, uniprot.dat, STRING×2, coxpresdb, TFLink, GAF —
all present with a name matching the production config, except rno
(uniprot loose at the root of `input_dir` instead of `uniprot/` — handled
via `AUX_FALLBACK_ROOTS[("rno","uniprot_dat")] = "INPUT_ROOT"`).

**hsa implemented and validated**:
`config/hsa/hsa_adapters_config_connected_sample.yaml` (new file, doesn't
collide with the `hsa_adapters_config_sample.yaml` used by tests) + new
`connected_sample` dataset tier in `species_config.yaml` +
`samples/hsa_connected/` (not `samples/hsa/`, which are the test fixtures).
Result: 180 genes, all ~7 core sources with real data (no synthetic
fallback — gencode_gtf=101,597, uniprot=179, gaf=18,432,
string_ppi/coexpr=12,804, coxpresdb=31,152, tflink=7,849). KG: 85,716 nodes,
136,276 edge rows. **179/180 genes in the giant component (43,843 nodes)** —
only `POLR2A` (ENSG00000181222) ended up isolated (gene+transcript+exon
only, no link to the rest); root cause: UniProt itself (`P24928`) has a
stale `DR Ensembl` cross-reference pointing at `ENSG00000047315` (an
old/different ID), not at GENCODE v49's current `ENSG00000181222` (which
even carries the official `reference_genome_error` tag) — a genuine
inconsistency between databases, not a generator bug. Accepted as an
isolated case (99.4% success), documented, not a blocker.

Measurement note: `check_kg_connectivity.py` uses `nx.Graph()` (not
`MultiGraph`), so the reported "resolved edges / total edges" is a **lower
bound** — node pairs with multiple real edge rows (e.g. GO's `is_a`)
collapse into 1 edge in the graph, artificially inflating the "unresolved"
count when naively compared to the row count. The metric that actually
matters is `dangling_instances` (an exact, pair-by-pair count) — for hsa
that was only 461 out of 136,276 (99.66%).

**mmu implemented and validated**:
`config/mmu/mmu_adapters_config_sample.yaml` (dataset tier `sample`, same
pattern as dmel — no name collision, mmu had no pre-existing test-sample
config). Result: 180 genes, core 100% real (gencode_gtf=14,289, uniprot=179,
gaf=10,349, string_ppi/coexpr=15,468, coxpresdb=31,862, tflink=2,518).
**180/180 genes in the giant component** — no isolated case, a perfect
result.

**rno implemented and validated**:
`config/rno/rno_adapters_config_sample.yaml` (dataset tier `sample`, rno had
no `sample` tier before — only `full`). Confirmed in practice that the
`AUX_FALLBACK_ROOTS[("rno","uniprot_dat")] = "INPUT_ROOT"` fallback works
(uniprot_sprot_rodents.dat.gz loose at the root of input_dir, found
correctly). Result: 180 genes, core 100% real (gencode_gtf=9,608,
uniprot=170, gaf=9,655, string_ppi/coexpr=7,508, coxpresdb=31,152,
tflink=352). KG: 42,777 nodes, 119,848 edges. **179/180 genes in the giant
component** (40,937 nodes) — only `Uba52` (ENSRNOG00000090524) ended up
isolated (gene+transcript+exon only). Root cause confirmed by grepping the
`.dat.gz` directly: **no** UniProt entry references that ENSRNOG via `DR
Ensembl` — a genuine data gap (a ubiquitin-ribosomal fusion gene, missing/
stale mapping in UniProt), same pattern as hsa's POLR2A case. Accepted as an
isolated case (99.4%), not a blocker.

**cel implemented and validated**:
`config/cel/cel_adapters_config_sample.yaml` (dataset tier `sample`
re-enabled in `species_config.yaml`, it had been commented out; kept the
pre-existing stub's `wbbt_anatomy`/`wbbt_subclass_of` adapters since they
cost nothing). Result: 180 genes, core 100% real (gencode_gtf=5,137,
uniprot=167, gaf=3,189, string_ppi/coexpr=10,236, coxpresdb=31,862,
tflink=3,017). KG: 40,877 nodes, 114,886 edges. **180/180 genes in the giant
component** (38,910 nodes) — no isolated case, a second perfect result
(tied with mmu).

Two generalization bugs found and fixed while implementing cel (the worm
breaks assumptions that had gone unnoticed in the other 4 species, all
Ensembl-style):
- **UniProt DR line**: cel uses `DR   EnsemblMetazoa;` (not `DR
  Ensembl;`), and the file `uniprot_sprot_invertebrates.dat.gz` is
  **shared** across several invertebrate species (Ciona, insects,
  molluscs...) — without an extra marker, `gene_to_protein`/
  `gene_to_uniprot` came out empty (0) for cel. Fixed with `uniprot_dr_db`
  (the DR database to match, per species) + `uniprot_gene_id_marker`
  (requires `"WBGene"` in the gene field, rejecting DR lines from other
  species in the same file) in `SPECIES["cel"]`, used by
  `_parse_uniprot_dat_gene_protein_uniprot()`.
- **Version suffix assumed on every ID**: `_parse_gtf_gene_to_transcript`,
  `_parse_uniprot_dat_gene_protein_uniprot` and
  `find_expansion_candidates` did `.split(".")[0]` or
  `.rsplit(".", 1)[-1]` on transcript/protein IDs assuming the
  Ensembl/GENCODE "ID.version" pattern (only 1 dot). Native WormBase IDs
  (`"T11F9.4a.1"`, from the GTF/STRING/DR line itself) have dots as part of
  the identifier — that assumption truncated it to `"T11F9"`, zeroing out
  `string_ppi`/`string_coexpression` (protein=0) even with the DR-line fix
  above already in place. Fixed with `strip_dot_version` (a per-species
  flag, `False` only for cel) in `_strip_dot_version()`, and replaced
  `find_expansion_candidates`'s `rsplit(".", 1)` with a taxon-prefix strip
  (`"6239."`) matching what `filter_string_links` already did correctly —
  this change is strictly more correct for the other 4 species too (it was
  only harmless there because their IDs have no internal dot), not a
  regression.
- Post-fix regression confirmed clean: dmel (same numbers as always), hsa
  and mmu (same `string_ppi`/`uniprot`/`gaf`/`coxpresdb` counts as before
  the fix) — both bugs only affected cel.

## State as of 2026-07-20: 5/5 species implemented and validated

| Species | Genes in giant component | Isolated case | Dataset tier / config |
|---|---|---|---|
| dmel | 180/180 (100%) | — | `sample` / `dmel_adapters_config_sample.yaml` |
| hsa | 179/180 (99.4%) | POLR2A (stale DR Ensembl in UniProt) | `connected_sample` / `hsa_adapters_config_connected_sample.yaml` |
| mmu | 180/180 (100%) | — | `sample` / `mmu_adapters_config_sample.yaml` |
| rno | 179/180 (99.4%) | Uba52 (no DR Ensembl in UniProt) | `sample` / `rno_adapters_config_sample.yaml` |
| cel | 180/180 (100%) | — | `sample` / `cel_adapters_config_sample.yaml` |

All 5 use the same generic engine (`scripts/generate_connected_sample.py` +
`check_kg_connectivity.py`), with dmel as the only one using
`build_global_maps` (FlyBase-specific backbone) — the other 4 use
`build_generic_backbone_maps`.

## Goal clarified (2026-07-24): full coverage, not a curated subset

The target is **full coverage** — every adapter block present in each
species' full `config/<species>/<species>_adapters_config.yaml` must
eventually be represented in that species' sample config with a
corresponding filtered sample file, not a curated subset. Don't treat any
adapter as "out of scope" or "deferred indefinitely" without checking back
— the working assumption is that everything in the full config eventually
needs a sample counterpart.

## Phases 4-8 (2026-07-26 to 2026-08-01): full coverage per species

Summary — full detail is in the assistant's project memory; here is only
what stays permanently useful for the algorithm/methodology.

- **Phase 4** (hsa): large batch of new hsa-only sources — GTEx, HOCOMOCO,
  dbSuper, Peregrine, dbsnp_snps, CADD, RefSeq closest-gene, TopLD (×4
  populations), EnhancerAtlas, cCRE, Catlas (directory-based, multi-file
  sampling — cCRE master TSV + per-cell-type ABC score files + cell
  ontology), motif_diff. Two **scope** bugs (not adapter bugs) found only
  after the first build showed low resolved-edge coverage: (1) the `snp`
  node universe (`dbsnp_snps`) was restricted to rsids referenced only by
  GWAS+GTEx, while ABC/CADD/RefSeq/TopLD reference a much larger rsid set —
  fixed by widening to the full union (`all_snp_rsids`); (2) HOCOMOCO (the
  source of `motif` nodes) was filtered to only the ~12 TFs whose own gene
  falls in the 180-gene closure, but motif_diff's edges reference HOCOMOCO's
  full ~400-TF universe regardless of closure membership — fixed by copying
  HOCOMOCO wholesale (it's only ~1.8MB).
- **Phase 5** (hsa): Roadmap Epigenomics (chromatin_state/h3_marks/dhs) —
  **corrected a wrong assumption from earlier phases** that no real data was
  available; the correct URLs were already documented as a comment inside
  the adapter files themselves
  (`biocypher_metta/adapters/hsa/roadmap_*_adapter.py`, `self.source_url`),
  only `hsa_data_source_config.yaml`'s were wrong/stale. **Lesson: before
  writing off a source as "unavailable," re-check the URL documented in the
  adapter's own code, not just the data-source-config comments** — those
  can go stale while the adapter itself stays correct. Also: the 4 missing
  GAF `_gene` blocks (the production config always had 8 GAF blocks —
  `_gene_product` and `_gene` variants — the sample only had the first 4),
  and Reactome's `reaction_to_pathway` + the 5 protein-role adapters
  (input/output/catalyst/negative/positive), which need
  `reactome_reaction_exporter_All_species.txt`.
- **Phase 6** (dmel/mmu/rno/cel): closed the same Reactome+GAF backlog for
  the other 4 species. Important finding: a good chunk of the block-by-block
  diff between production config and sample was a **false positive** —
  different block names for the same data/adapter (e.g.
  `alliance_gene_disease_biomarker_orthology` ≈ the already-existing
  `alliance_gene_disease_biomarker_via_orthology`; dmel's `coexpression` ≈
  the already-existing `coxpresdb_coexpression`). Real and closed:
  `uniprot_dbxref_bgee_gene` (real data already exists in all 4 species'
  filtered `.dat.gz` — a case-sensitivity bug in my own verification `grep`
  had made it look like 0: `DR   BGEE` vs. the real `DR   Bgee`),
  `emapa_anatomy`/`emapa_subclass_of` (mmu, full ontology, same pattern as
  BTO/CL/UBERON). The public reactome.org
  `reactome_reaction_exporter_All_species.txt` is **human-only** (confirmed:
  133,223 rows, all R-HSA) — worked for hsa in Phase 5 by scope coincidence,
  but blocked the other 4 species. Closed in Phase 7 when a genuinely
  multi-species re-export (from an as-yet-unmerged Reactome PR, 718,196 rows
  covering R-HSA/R-MMU/R-RNO/R-DME/R-CEL/others) became available —
  replicated (hard-linked) to the remaining 4 species. Two real bugs found
  in this phase: (1) the synthetic-record generator
  `write_synthetic_uniprot_chebi_part_of()` had `taxon_id` hardcoded to
  9606 — worked for hsa by coincidence, silently zeroed out for the other 4
  species (the adapter skips records whose taxonomy doesn't match); fixed
  by parameterizing. **Lesson: every "write a fake record for adapter X"
  helper must parameterize every field the adapter branches on, not just
  the ones the first call site needed.** (2) hsa-specific code (Reactome GO
  cross-links) had been placed by mistake in the function shared across
  species in Phase 5 — a latent bug that would have broken dmel/mmu/rno/cel
  with a `KeyError` the moment they were re-run; fixed by moving it to the
  hsa-only function.
- **Phase 7** (dmel): `fca2`/`afca` (expression) were missing from
  `samples/dmel/` with no raw source anywhere — restored first as a small
  static copy, then with the **real, full raw source** supplied by the user
  (fca2 154MB, afca 6.2GB) and a proper filter added: `filter_afca()` is the
  reference example of the "symbol/alt-ID resolution" pattern (see §methodology
  below) — the source keys on gene symbol, not FBgn, so the filter
  replicates exactly the same symbol→FBgn resolution the real adapter
  (`ExpressionValueAdapter`) already does internally via
  `fbgn_fbtr_fbpp_expanded.tsv.gz`.
- **Phase 8** (hsa): `samples/hsa_connected/` had grown to 9.9GB, almost
  entirely (9.6GB) from `motif_diff`. See "Shrinking a large source" below —
  it's the reference case for this technique.

## How to add a new sample source (general checklist)

The graph doesn't need every new node/edge type to invent its own
connectivity mechanism. It needs every row you keep to land on a node
that's **already reachable from the closure** — either directly, or via
something *derived* from the closure by an earlier, already-filtered step.
Connectivity is a property of the whole pipeline's filtering, not something
each new source has to re-establish on its own.

**Step 1 — read the adapter.** What does `get_nodes()`/`get_edges()`
actually yield (`source, target, label`), and which column does that come
from? That tells you which `id_sets[...]` to filter against (`gene`,
`protein`, `transcript`, `uniprot`, `entrez`, ...) — computed once per
species from the closure, before any adapter-specific filtering runs.

**Step 2 — pick the filtering pattern**, in the order worth trying:

1. **Direct membership** (the common case): the row's gene/protein/
   transcript column is already in `id_sets`. Reuse an existing generic
   helper before writing a bespoke filter: `filter_by_column_in_set`,
   `filter_by_two_columns_in_set`, `filter_by_any_column_in_set`,
   `copy_header_comments_and_filter` (the primitive underneath the other
   three — write a custom `keep_fn` if none of the wrappers fit, but keep
   using this primitive for the `#`-comment/blank-line passthrough
   behavior).

2. **Derived/propagated ID sets** (two-hop): the new source references
   something that isn't itself a gene but is *reachable* from one via an
   already-filtered file. Filter the first hop by `id_sets["gene"]`, capture
   the resulting IDs it touched, then filter every *subsequent* file
   referencing that same entity type against the derived set — not against
   genes directly. Example: Reactome as a whole (see
   `filter_reactome_gene_pathway_or_reaction`,
   `filter_reactome_reaction_exporter` — the latter requires both the
   pathway id AND the reaction id already in their respective retained sets,
   so neither side of the edge can ever dangle).

3. **Borrowed relevance via a shared cache** (rsid-anchored sources with no
   gene column at all): ABC/CADD/RefSeq/GTEx/Roadmap/motif_diff key on SNP
   rsid, with nothing to filter on directly. Build one shared
   `all_snp_rsids` set from the sources that DO resolve rsid → gene
   (ABC/CADD/RefSeq/GTEx), then filter every other rsid-keyed source by
   membership in that set. The same cache (`build_dbsnp_cache`) also backs a
   small SQLite `dbsnp_rsid_map`/`dbsnp_pos_map` built directly from those
   rows (not a genome-wide VCF) for adapters that need rsid → chr/pos at
   KG-build time.

4. **Wholesale copy** for compact standalone ontologies (GO, BTO, CL,
   UBERON, HOCOMOCO, CHEBI, EMAPA) — these load in full regardless of gene
   closure, matching the full production config. No filtering, no scoping.

5. **Symbol/alt-ID resolution**: if the new source keys on something other
   than your closure's primary ID space (e.g. dmel's `afca` file keys on
   gene *symbol*, not FBgn), don't invent a new mapping — read the real
   adapter's own resolution code and replicate it exactly (same aux file,
   same column indices), so the filtered subset matches precisely what the
   adapter would keep at build time. See `filter_afca()`.

**Step 3 — keep node and edge files in sync.** If an edge references a
derived ID (pathway, reaction, chebi term), the corresponding *node* file
must be filtered by that same retained-ID set (or wholesale-loaded), so
nothing ends up pointing at a node absent from the sample. This is silent
dangling — it doesn't throw, `Neo4jCSVWriter`'s `MATCH`-based Cypher just
quietly fails to create the edge. `scripts/check_kg_connectivity.py`'s
"dangling types/instances" report catches this after the fact — don't rely
on it catching everything preemptively; check the node/edge pairing logic
yourself when wiring a new derived-ID source.

**Step 4 — if real data comes back empty, verify before assuming it's a
genuine gap.

Only after confirming genuine absence: leave it empty with
`record(..., empty_ok_reason="...")`, or — if the edge type is
schema-important enough to want at least one real example in the sample —
write a minimal synthetic record using **real, already-present node IDs**
(never fabricate new nodes), clearly tagged and logged as synthetic. See
`write_synthetic_tsv()` (simple TSV/line-based formats) and
`write_synthetic_uniprot_chebi_part_of()` (a full worked example: a
hand-built, `Bio.SwissProt`-parseable minimal SwissProt record with a
`BINDING` feature carrying `/ligand_id` + `/ligand_part_id` qualifiers
pointing at two real CHEBI ids already in the sample's own Reactome
small-molecule output). **Any such "write a fake record for adapter X"
helper must parameterize every field the adapter branches on** (taxon_id
included) — don't hardcode anything that only happens to match for one
species.

**Step 5 — shrinking a large source.** See the dedicated section below —
don't assume it's connectivity-critical without testing first.

**Step 6 — always validate at the end, never assume.** Regenerate the
sample for the affected species → rebuild the KG → run the connectivity
checker (§6) → confirm every anchor/closure gene is still in the single
largest connected component. Every phase in this branch's history has been
gated on this check passing, including after changes that *could*
plausibly affect it but don't look related at first glance (e.g. moving
code between per-species and shared generation functions).

## Wiring checklist (mechanical steps once you've picked a filter pattern)

1. Add a `SOURCE_FILES[species][name] = (glob_pattern, stable_output_name)`
   entry. Use a glob (`*`) wherever the real filename embeds a release
   version; use the exact name otherwise (e.g. non-versioned files like
   `afca_afca_annotation_group_by_mean.tsv.gz`).
2. Write or reuse a filter function (§ above) and call it + `record(name, n,
   empty_ok_reason=..., synthetic_info=...)` from the right generation
   function:
   - `_generate_core` — universal, every species (GENCODE/UniProt/GAF/
     STRING/coxpresdb/TFLink).
   - `_generate_phase2_shared` — universal, applies the same way per species
     (bgee, alliance, epd, reactome core, rna_central).
   - `_generate_dmel_extension` / `_generate_hsa_phase2_extension` —
     species-specific adapters that only exist for that one species. (mmu/
     rno/cel don't currently have their own extension functions — if a
     species-specific adapter shows up for one of them, add a
     `_generate_<species>_extension` following the same shape.)
3. Add the config block to
   `config/<species>/<species>_adapters_config_sample.yaml`. All 5 species
   (including hsa, since the hsa-directory-unification migration) use this
   same single file, which is also the default `test/conftest.py` uses —
   there is no separate connected-sample file anymore. `filepath:` must
   point at the **stable** name from step 1, not the versioned real
   filename.
4. Run `python scripts/generate_connected_sample.py --species <species>`,
   then rebuild (`create_knowledge_graph.py --species <species> --dataset
   sample --output-dir output_<species>_sample --writer-type neo4j
   --no-checkpoint`), then validate (§6). hsa's `dbsnp_cache_root`
   (`aux_files/hsa/sample_dbsnp`) is already wired into
   `config/species_config.yaml`'s `sample` tier, so no extra dbsnp flags are
   needed on the command line.

## Shrinking a large source without breaking connectivity

If a source is large, **don't assume it sustains connectivity — test it.**
Before pruning/capping anything: with the source already included in a
build, run the connectivity checker *excluding* that edge type's output
file from the union-find pass, and see if the closure still holds. If it
does, you have full freedom to cap it aggressively rather than keeping
every matching row.

Reference case: `motif_diff` (hsa). The filtered sample had 915K rows but
was 10GB — each row carries 770 TF-model score columns (~10.8KB/row), so
**row count**, not filtering scope, drove the size. Excluding
`edges_snp_alters_binding_motif.csv` entirely from the union-find pass and
re-checking confirmed the 180 genes stay in the largest component —
SNP↔gene connectivity flows through ABC/CADD/RefSeq/TopLD/GWAS/eQTL, not
motif_diff. `filter_motif_diff()` therefore gained a `max_rows` parameter
applied as an **early-exit during the scan** (stop reading the 10GB source
the moment the cap is hit — much faster than reading all of it and
truncating after), not post-hoc truncation. Capping hsa to 2000 rows took
`samples/hsa_connected/` from 9.9GB → 338MB, KG edges from 45.3M → 13.3M,
the `motif_diff_tf_snp` build step from 1h5m → 8.5s, and Neo4j load time
from ~25min → ~7min — with 180/180 genes still connected. This is the
template for any future "shrink X while staying connected" request: never
assume a large file is connectivity-critical without testing; most
breadth/coverage sources (especially ones with no `dbsnp_rsid_map`/gene-
anchored filter of their own) are prunable near-arbitrarily once the
closure's core connectivity comes from elsewhere.

## Validating connectivity without freezing the machine

`scripts/check_kg_connectivity.py` (`load_graph()`) is the general-purpose
node/edge loader + dangling-reference reporter, for any species'
`output_<species>_.../` build directory. It loads the *entire* graph,
though — **do not** build a second, ad-hoc validation script around
`networkx` (`nx.Graph()`) just to check closure membership. A prior version
of this specific check did exactly that and caused two machine freezes:
networkx's per-edge attribute dict + per-node adjacency dict-of-dicts
overhead grew past 10GB RSS on a 6.1M-node/45M-edge graph while
system-wide available memory collapsed toward zero (this machine has only
2GB swap against 31GB RAM, so a bad spike hard-freezes the OS via swap
thrashing rather than degrading gracefully).

Use a streaming **Union-Find (disjoint-set)** instead: map every node key
to a small int in one pass over the node CSVs (the only thing that scales
with node count — no adjacency structure kept), then stream every edge file
once, `union()`-ing endpoint pairs and discarding the row immediately.
Always set a hard `resource.setrlimit(resource.RLIMIT_AS, ...)` cap at the
top of any such script as a safety net, so a wrong estimate kills the
*process* with a clean `MemoryError` instead of taking the whole machine
down. This checks the same 180/180-gene-closure question on a
6.1M-node/45.3M-edge graph while staying under ~1GB RSS (vs. 10GB+ and a
frozen machine with the networkx version).

To test whether a specific edge type is connectivity-critical before
pruning it (see above), add an `exclude_substr` parameter that skips any
edge file whose path contains it during the union-find pass, then compare
## Shared MotifDiff Adapter integration (Human + Drosophila)

### Overview

`biocypher_metta/adapters/motif_diff_adapter.py` consolidates the species-specific MotifDiff adapters into a unified, species-agnostic adapter. It reads a MotifDiff `.diff` matrix (where rows represent genomic variants and columns represent transcription factor binding motifs) and yields `tf_snp` edges (`snp to motif association` in schema, written as `alters_binding` in graph outputs).

### Operating Modes

1. **Human mode (`motif_map_path=None`)**:
   - Parses human HOCOMOCO model identifiers (e.g. `ANDR_HUMAN.H11MO.0.A`).
   - Resolves HGNC symbols to Ensembl Gene IDs (`ENSEMBL:<id>`) via `HGNCProcessor`.
   - Formats variants as `DBSNP:<rs>`.
   - Selects the best motif per gene based on HOCOMOCO confidence tier rank (`A > B > C > D`) and absolute score.
   - Retains exact backward compatibility across all combinations of `write_properties` and `add_provenance`.

2. **Drosophila / mapping-based mode (`motif_map_path` specified)**:
   - Resolves motif column headers (e.g. `MA0010.1 br`, `MA0530.1 EcR::usp`) using a TSV mapping file (e.g. `jaspar/dmel_motif_map_v4.tsv`).
   - Maps motifs to FlyBase gene IDs (e.g. `FBgn0000210`).
   - Grouping: selects the motif column with the highest absolute MotifDiff score for each TF gene.
   - Heterodimers: controlled by `heterodimer: "all"` (expands to both partner genes) or `"skip"` (omits heterodimer motifs).
   - Formats variants as bare `RS<digits>` (`variant_prefix: ""`, `uppercase_variant: True`) and genes as bare `FBgn...` (`gene_prefix: ""`).
   - Taxon ID: configured as `7227` (`taxon_id: 7227`, threshold: `0.001`).

### Shared Schema Specification

The `snp to motif association` edge schema is defined in `config/primer_schema_config.yaml`:
- `input_label`: `tf_snp`
- `output_label`: `alters_binding`
- `source`: `snp`
- `target`: `motif`
- `properties`:
  - `effect`: `gain` or `loss`
  - `score`: float (MotifDiff score)
  - `motif_id`: str (motif accession, e.g. `MA0010.2`)
  - `taxon_id`: int (`9606` for human, `7227` for dmel)
  - `source`: str (`primary_knowledge_source`)
  - `source_url`: str (`source_web_page`)

### Production Dataset Instructions (Bizon Server)

The full Drosophila merged dataset is located on the Bizon storage server:
`/mnt/hdd_1/biocypher-kg/input/dmel/motif_diff/_mono_probNorm_average.diff`

To merge the 22 variant chunks:
```bash
cd /mnt/hdd_1/biocypher-kg/input/dmel/jaspar
# 1. Verify all chunk headers match
for f in chunks/v_*average.diff; do head -1 "$f" | md5sum; done | sort | uniq -c
# Expected: exactly 1 line with count 22

# 2. Concatenate into merged output
mkdir -p ../motif_diff
out=../motif_diff/_mono_probNorm_average.diff
head -1 chunks/v_000_*average.diff > "$out"
for f in chunks/v_*average.diff; do tail -n +2 "$f" >> "$out"; done

# 3. Validate line count
wc -l "$out"
# Expected: 5,327,113 lines (1 header + 5,327,112 variants across 357 motif columns)
```

### Running and Validating

- Run unit & integration tests:
  ```bash
  pytest test/test_motif_diff_adapter.py -v
  ```
- Build Drosophila sample graph:
  ```bash
  python create_knowledge_graph.py --species dmel --dataset sample --writer-type neo4j --no-checkpoint
  ```

## Open / to revisit

- `gene-tf_snp-snp` (hsa): confirmed to be a legacy/different-version schema
  shape on a remote comparison database, not a gap in the current schema's
  design — do not implement without checking with the user first (would
  mean reintroducing an edge shape the current schema deliberately replaced
  with `snp-alters_binding-motif`).
- `dmel_snp_located_in_gene`: needs a real schema decision (a new
  `snp_in_gene`-labeled association added to `dmel_schema_config.yaml`, or
  confirm it should just reuse `variant_of` and the production config's
  separate block is redundant/a mistake) — check with the user before
  touching, it's a schema-design question, not just wiring a block. The
  adapter's own code confirms (a comment in
  `biocypher_metta/adapters/dmel/allele_adapter.py`) that it was only ever
  built/tested for the `variant_of` case.
- Bgee's expression edges (`gene_expressed_in_anatomical_entity` and
  similar) still dangle for mmu/rno/cel in some cases (no matching
  anatomy/developmental-stage ontology loaded for that specific edge's
  expected node type) — not chased down in Phases 6-8, pre-existing.
- Original Phase 2 backlog (bgee, rna_central, reactome, chebi_ontology,
  epd, alliance, hsa-only GWAS/dbSNP/ENCODE/ABC/Catlas/Roadmap/HPO) —
  **completed** in Phases 4-8 above; kept here only as a historical record
  that it was once "open."
- dmel phenotype/disease/genetic-interaction coverage hasn't been tested
  specifically yet — may need 2-3 well-studied classic genes (`white`,
  `hh`, `per`) if this is ever needed beyond what's already covered.
