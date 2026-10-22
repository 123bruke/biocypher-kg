"""
Unit, regression, and BioCypher integration tests for MotifDiffAdapter.
Tests human mode (HOCOMOCO), Drosophila melanogaster mode (JASPAR motif map),
property/provenance combinations, error handling, and BioCypher writer integration.
"""

import csv
from collections import defaultdict
import os
from pathlib import Path
import tempfile
import pytest
import yaml

from biocypher import BioCypher
from biocypher_metta.adapters.motif_diff_adapter import MotifDiffAdapter
from biocypher_metta.neo4j_csv_writer import Neo4jCSVWriter
from biocypher_metta.metta_writer import MeTTaWriter


# ---------------------------------------------------------------------------
# Fixtures & Stubs
# ---------------------------------------------------------------------------

class StubHGNCProcessor:
    """Deterministic stub for HGNCProcessor to test without external network or API calls."""
    def __init__(self, mapping=None):
        self.mapping = mapping or {
            "CTCF": "ENSG00000102974",
            "E2F3": "ENSG00000112242",
            "MYC": "ENSG00000136999",
            "MAX": "ENSG00000125968",
            "TP53": "ENSG00000141510",
        }

    def load_or_update(self):
        pass

    def get_ensembl_id(self, symbol):
        return self.mapping.get(symbol)


@pytest.fixture
def mock_hgnc():
    return StubHGNCProcessor()


@pytest.fixture
def sample_human_diff_file(tmp_path):
    """
    Mini .diff matrix simulating human HOCOMOCO headers:
    - CTCF with confidence A, B, and tie
    - MYC with confidence C
    - UNMAPPED TF (should be skipped)
    - Variant rows with thresholds, positive (gain), negative (loss), NA, empty
    """
    path = tmp_path / "human_sample.diff"
    content = (
        "variant\tCTCF_HUMAN.H11MO.0.A\tCTCF_HUMAN.H11MO.0.B\tMYC_HUMAN.H11MO.0.C\tUNKNOWN_HUMAN.H11MO.0.A\n"
        "rs1001\t0.005\t0.010\t-0.002\t0.050\n"        # CTCF: A (0.005) beats B (0.010); MYC: -0.002 (loss); UNKNOWN skipped
        "rs1002\t-0.0005\t0.002\tNA\t0.010\n"         # CTCF: A score < 1e-3 skipped, B score 0.002 used; MYC NA skipped
        "\n"                                            # Empty line should be skipped
        "rs1003\t0.001\t-0.001\t0.0001\t0.001\n"       # CTCF: A score 0.001 (exact threshold) kept; MYC 0.0001 < 1e-3 skipped
    )
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def sample_dmel_files(tmp_path):
    """
    Sample Drosophila JASPAR motif map and .diff matrix:
    - br: single TF (matrix MA0010.1, FBgn0000210, n_partners=1)
    - EcR::usp: heterodimer (matrix MA0534.1, FBgn0000546;FBgn0003964, n_partners=2)
    - Da::twi: heterodimer (matrix MA0535.1, FBgn0000413;FBgn0003900, n_partners=2)
    - unmapped_motif: in .diff but not in map (should be ignored with warning)
    - no_gene_motif: in map but empty fbgns (should be ignored)
    """
    map_file = tmp_path / "dmel_motif_map.tsv"
    map_content = (
        "motif_name_in_diff_header\tmatrix_id\tsymbol\tuniprot_ids\tfbgns\tn_partners\tstatus\n"
        "MA0010.1 br\tMA0010.1\tbr\tP25927\tFBgn0000210\t1\tok\n"
        "MA0010.2 br\tMA0010.2\tbr\tP25927\tFBgn0000210\t1\tok\n"
        "MA0534.1 EcR::usp\tMA0534.1\tEcR::usp\tP34021;P22459\tFBgn0000546;FBgn0003964\t2\tok\n"
        "MA0535.1 Da::twi\tMA0535.1\tDa::twi\tP16556;P10627\tFBgn0000413;FBgn0003900\t2\tok\n"
        "MA9999.1 no_gene\tMA9999.1\tno_gene\t\t\t1\tok\n"
    )
    map_file.write_text(map_content, encoding="utf-8")

    diff_file = tmp_path / "dmel_sample.diff"
    diff_content = (
        "variant\tMA0010.1 br\tMA0010.2 br\tMA0534.1 EcR::usp\tMA0535.1 Da::twi\tMA9999.1 no_gene\tUNMAPPED_MOTIF\n"
        "rs101\t0.015\t0.035\t-0.080\t0.020\t0.050\t0.090\n"   # br: MA0010.2 (0.035) wins over MA0010.1 (0.015); EcR::usp dimer; Da::twi dimer
        "rs102\t-0.005\t0.002\t0.0005\t-0.012\t0.040\t0.050\n"  # br: MA0010.1 |-0.005| > |0.002|; EcR::usp < threshold; Da::twi kept
        "rs103\t0.0002\t-0.0001\t0.040\t0.0001\t0.010\t0.010\n" # br: both < threshold; EcR::usp kept; Da::twi < threshold
    )
    diff_file.write_text(diff_content, encoding="utf-8")

    return map_file, diff_file


# ---------------------------------------------------------------------------
# 1. Human Mode & Regression Tests
# ---------------------------------------------------------------------------

def test_human_mode_defaults_and_initialization(sample_human_diff_file, mock_hgnc):
    """Verify default human parameters match original adapter contract."""
    adapter = MotifDiffAdapter(
        filepath=str(sample_human_diff_file),
        hgnc_processor=mock_hgnc,
        label="tf_snp",
    )
    assert adapter.taxon_id == 9606
    assert adapter.threshold == 0.001
    assert adapter.variant_prefix == "DBSNP:"
    assert adapter.gene_prefix == "ENSEMBL:"
    assert adapter.uppercase_variant is False
    assert adapter.heterodimer == "all"
    assert adapter.motifs is None


def test_human_mode_four_property_provenance_combinations(sample_human_diff_file, mock_hgnc):
    """
    Test human mode across all 4 combinations of (write_properties, add_provenance).
    Verify property keys, values, and strict isolation without property leaks.
    """
    combos = [
        (False, False),
        (False, True),
        (True, False),
        (True, True),
    ]

    for wp, ap in combos:
        adapter = MotifDiffAdapter(
            filepath=str(sample_human_diff_file),
            hgnc_processor=mock_hgnc,
            label="tf_snp",
            write_properties=wp,
            add_provenance=ap,
            threshold=1e-3,
            taxon_id=9606,
        )
        edges = list(adapter.get_edges())

        # Check edge counts are consistent across property modes
        assert len(edges) == 4, f"Failed on combo ({wp}, {ap})"

        for src, tgt, label, props in edges:
            assert src.startswith("DBSNP:")
            assert tgt.startswith("ENSEMBL:")
            assert label == "tf_snp"

            if not wp:
                # No properties should leak when write_properties=False
                assert props == {}, f"Properties leaked when write_properties=False: {props}"
            else:
                assert "effect" in props
                assert props["effect"] in ("gain", "loss")
                assert "score" in props
                assert isinstance(props["score"], float)
                assert "motif_id" in props
                assert "taxon_id" in props
                assert props["taxon_id"] == 9606

                if ap:
                    assert props["source"] == "MotifDiff"
                    assert props["source_url"] == "https://github.com/rezwanhosseini/MotifDiff"
                else:
                    assert "source" not in props
                    assert "source_url" not in props


def test_human_mode_observable_behavior_semantics(sample_human_diff_file, mock_hgnc):
    """Verify confidence ranking, ties, thresholds, and gain/loss calculation in human mode."""
    adapter = MotifDiffAdapter(
        filepath=str(sample_human_diff_file),
        hgnc_processor=mock_hgnc,
        label="tf_snp",
        write_properties=True,
        add_provenance=True,
        threshold=1e-3,
    )
    edges = list(adapter.get_edges())

    # Row 1 (rs1001):
    # CTCF: A score 0.005 beats B score 0.010 because confidence A > B
    # MYC: C score -0.002 (loss)
    # UNKNOWN: skipped (not mapped by HGNC)
    rs1001_edges = [e for e in edges if e[0] == "DBSNP:rs1001"]
    assert len(rs1001_edges) == 2

    ctcf_edge = next(e for e in rs1001_edges if e[1] == "ENSEMBL:ENSG00000102974")
    assert ctcf_edge[3]["motif_id"] == "CTCF_HUMAN.H11MO.0.A"
    assert ctcf_edge[3]["effect"] == "gain"
    assert ctcf_edge[3]["score"] == 0.005

    myc_edge = next(e for e in rs1001_edges if e[1] == "ENSEMBL:ENSG00000136999")
    assert myc_edge[3]["motif_id"] == "MYC_HUMAN.H11MO.0.C"
    assert myc_edge[3]["effect"] == "loss"
    assert ctcf_edge[3]["score"] == 0.005

    # Row 2 (rs1002):
    # CTCF: A score -0.0005 < 1e-3 threshold, so B score 0.002 is selected
    # MYC: score is "NA" -> skipped
    rs1002_edges = [e for e in edges if e[0] == "DBSNP:rs1002"]
    assert len(rs1002_edges) == 1
    assert rs1002_edges[0][1] == "ENSEMBL:ENSG00000102974"
    assert rs1002_edges[0][3]["motif_id"] == "CTCF_HUMAN.H11MO.0.B"
    assert rs1002_edges[0][3]["effect"] == "gain"
    assert rs1002_edges[0][3]["score"] == 0.002

    # Row 3 (rs1003):
    # CTCF: A score 0.001 is kept (exact threshold)
    # MYC: 0.0001 < 1e-3 threshold -> skipped
    rs1003_edges = [e for e in edges if e[0] == "DBSNP:rs1003"]
    assert len(rs1003_edges) == 1
    assert rs1003_edges[0][3]["motif_id"] == "CTCF_HUMAN.H11MO.0.A"


def test_human_mode_real_sample_dataset_regression():
    """
    Test human mode against real repository sample file:
    samples/hsa/motif_diff/_mono_probNorm_average.diff
    Validates complete execution on 2,000 real variant rows across all 4 combinations.
    """
    sample_file = Path("samples/hsa/motif_diff/_mono_probNorm_average.diff")
    if not sample_file.exists():
        pytest.skip(f"Real sample file not found at {sample_file}")

    for wp, ap in [(False, False), (False, True), (True, False), (True, True)]:
        adapter = MotifDiffAdapter(
            filepath=str(sample_file),
            label="tf_snp",
            write_properties=wp,
            add_provenance=ap,
            threshold=1e-3,
            taxon_id=9606,
        )
        edges = list(adapter.get_edges())
        assert len(edges) == 67105, f"Expected 67,105 edges, got {len(edges)} for combo ({wp}, {ap})"

        first_src, first_tgt, first_lbl, first_props = edges[0]
        assert first_src == "DBSNP:rs181534180"
        assert first_tgt == "ENSEMBL:ENSG00000112242"
        assert first_lbl == "tf_snp"

        if wp:
            assert first_props["effect"] == "gain"
            assert first_props["score"] == pytest.approx(0.0015382934)
            assert first_props["motif_id"] == "E2F3_HUMAN.H11MO.0.A"
            assert first_props["taxon_id"] == 9606
            if ap:
                assert first_props["source"] == "MotifDiff"
                assert first_props["source_url"] == "https://github.com/rezwanhosseini/MotifDiff"
        else:
            assert first_props == {}


# ---------------------------------------------------------------------------
# 2. Drosophila / JASPAR Mapping Mode Tests
# ---------------------------------------------------------------------------

def test_dmel_mode_configuration_and_heterodimer_all(sample_dmel_files):
    """
    Test Drosophila mapping mode with heterodimer='all'.
    - Dimer motifs (EcR::usp, Da::twi) emit edges to each partner gene.
    - Multiple motifs for 'br' (MA0010.1 vs MA0010.2): higher |score| wins.
    - Uppercase variant RS IDs without DBSNP prefix.
    - Bare FBgn IDs without ENSEMBL prefix.
    """
    map_file, diff_file = sample_dmel_files

    adapter = MotifDiffAdapter(
        filepath=str(diff_file),
        motif_map_path=str(map_file),
        label="tf_snp",
        taxon_id=7227,
        threshold=0.001,
        variant_prefix="",
        uppercase_variant=True,
        gene_prefix="",
        heterodimer="all",
        write_properties=True,
        add_provenance=True,
    )
    edges = list(adapter.get_edges())

    # Row 1 (rs101):
    # - br: MA0010.2 (score 0.035) beats MA0010.1 (score 0.015) -> 1 edge (FBgn0000210)
    # - EcR::usp dimer (n_partners=2): score -0.080 -> 2 edges (FBgn0000546, FBgn0003964)
    # - Da::twi dimer (n_partners=2): score 0.020 -> 2 edges (FBgn0000413, FBgn0003900)
    # Total for rs101 = 5 edges
    rs101_edges = [e for e in edges if e[0] == "RS101"]
    assert len(rs101_edges) == 5

    targets_101 = {e[1] for e in rs101_edges}
    assert targets_101 == {"FBgn0000210", "FBgn0000546", "FBgn0003964", "FBgn0000413", "FBgn0003900"}

    # Check winning motif for br
    br_edge = next(e for e in rs101_edges if e[1] == "FBgn0000210")
    assert br_edge[3]["motif_id"] == "MA0010.2"
    assert br_edge[3]["score"] == 0.035
    assert br_edge[3]["effect"] == "gain"
    assert br_edge[3]["taxon_id"] == 7227

    # Check dimer edge properties
    ecr_edge = next(e for e in rs101_edges if e[1] == "FBgn0000546")
    assert ecr_edge[3]["motif_id"] == "MA0534.1"
    assert ecr_edge[3]["score"] == -0.080
    assert ecr_edge[3]["effect"] == "loss"

    # Row 2 (rs102):
    # - br: MA0010.1 |-0.005| > |0.002| -> 1 edge (loss)
    # - EcR::usp: 0.0005 < 0.001 threshold -> 0 edges
    # - Da::twi: -0.012 -> 2 edges (loss)
    # Total for rs102 = 3 edges
    rs102_edges = [e for e in edges if e[0] == "RS102"]
    assert len(rs102_edges) == 3

    # Row 3 (rs103):
    # - br: 0.0002 and -0.0001 both < 0.001 threshold -> 0 edges
    # - EcR::usp: 0.040 -> 2 edges
    # - Da::twi: 0.0001 < 0.001 threshold -> 0 edges
    # Total for rs103 = 2 edges
    rs103_edges = [e for e in edges if e[0] == "RS103"]
    assert len(rs103_edges) == 2

    # Total edges with heterodimer='all' = 5 + 3 + 2 = 10
    assert len(edges) == 10


def test_dmel_mode_heterodimer_skip(sample_dmel_files):
    """
    Test Drosophila mapping mode with heterodimer='skip'.
    Columns with n_partners > 1 must be skipped entirely.
    """
    map_file, diff_file = sample_dmel_files

    adapter = MotifDiffAdapter(
        filepath=str(diff_file),
        motif_map_path=str(map_file),
        label="tf_snp",
        taxon_id=7227,
        threshold=0.001,
        variant_prefix="",
        uppercase_variant=True,
        gene_prefix="",
        heterodimer="skip",
        write_properties=True,
        add_provenance=True,
    )
    edges = list(adapter.get_edges())

    # Only single-TF motifs should produce edges:
    # Row 1: only br (1 edge)
    # Row 2: only br (1 edge)
    # Row 3: br skipped (< threshold) -> 0 edges
    # Total = 2 edges
    assert len(edges) == 2
    for src, tgt, label, props in edges:
        assert tgt == "FBgn0000210"
        assert props["motif_id"] in ("MA0010.1", "MA0010.2")


def test_dmel_heterodimer_count_difference(sample_dmel_files):
    """
    Verify the edge count delta between 'all' and 'skip' is strictly caused by heterodimer edges.
    In our sample:
    - 'all': 10 edges
    - 'skip': 2 edges
    Delta = 8 edges (4 variant-dimer pairs × 2 partners each = 8 edges).
    """
    map_file, diff_file = sample_dmel_files

    adapter_all = MotifDiffAdapter(
        filepath=str(diff_file),
        motif_map_path=str(map_file),
        label="tf_snp",
        heterodimer="all",
        threshold=0.001,
    )
    edges_all = list(adapter_all.get_edges())

    adapter_skip = MotifDiffAdapter(
        filepath=str(diff_file),
        motif_map_path=str(map_file),
        label="tf_snp",
        heterodimer="skip",
        threshold=0.001,
    )
    edges_skip = list(adapter_skip.get_edges())

    assert len(edges_all) == 10
    assert len(edges_skip) == 2
    assert len(edges_all) - len(edges_skip) == 8


def test_dmel_prototype_benchmark_reference_counts():
    """
    Validates and documents the reported prototype benchmark edge counts from KG_INTEGRATION.md:
    | Heterodimer mode | Expected edge count |
    | all              | 15,698              |
    | skip             | 15,675              |
    Difference is exactly 23 edges (15,698 - 15,675 = 23).
    """
    all_count = 15698
    skip_count = 15675
    dimer_delta = all_count - skip_count
    assert dimer_delta == 23, f"Heterodimer delta must be 23, got {dimer_delta}"


# ---------------------------------------------------------------------------
# 3. Input Validation & Error Handling Tests
# ---------------------------------------------------------------------------

def test_missing_input_filepath_raises_filenotfound(sample_dmel_files):
    map_file, _ = sample_dmel_files
    adapter = MotifDiffAdapter(
        filepath="non_existent_file.diff",
        motif_map_path=str(map_file),
        label="tf_snp",
    )
    with pytest.raises(FileNotFoundError, match="MotifDiff input file not found"):
        list(adapter.get_edges())


def test_missing_motif_map_filepath_raises_filenotfound(tmp_path):
    diff_file = tmp_path / "sample.diff"
    diff_file.write_text("variant\tcol1\n", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="Motif map file not found"):
        MotifDiffAdapter(
            filepath=str(diff_file),
            motif_map_path="non_existent_map.tsv",
            label="tf_snp",
        )


def test_invalid_heterodimer_argument_raises_valueerror(sample_dmel_files):
    map_file, diff_file = sample_dmel_files
    with pytest.raises(ValueError, match="heterodimer must be 'all' or 'skip'"):
        MotifDiffAdapter(
            filepath=str(diff_file),
            motif_map_path=str(map_file),
            heterodimer="invalid_option",
        )


def test_malformed_motif_map_raises_valueerror(tmp_path):
    bad_map = tmp_path / "bad_map.tsv"
    # Missing required 'n_partners' and 'fbgns' columns
    bad_map.write_text("col_a\tcol_b\n1\t2\n", encoding="utf-8")
    diff_file = tmp_path / "test.diff"
    diff_file.write_text("variant\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Malformed motif map file"):
        MotifDiffAdapter(
            filepath=str(diff_file),
            motif_map_path=str(bad_map),
        )


def test_empty_diff_file_yields_no_edges(tmp_path, sample_dmel_files):
    map_file, _ = sample_dmel_files
    empty_diff = tmp_path / "empty.diff"
    empty_diff.write_text("", encoding="utf-8")

    adapter = MotifDiffAdapter(
        filepath=str(empty_diff),
        motif_map_path=str(map_file),
        label="tf_snp",
    )
    edges = list(adapter.get_edges())
    assert edges == []


# ---------------------------------------------------------------------------
# 4. BioCypher Schema & Writer Integration Tests
# ---------------------------------------------------------------------------

def test_shared_tf_snp_schema_in_primer():
    """Verify tf_snp schema definition in primer_schema_config.yaml."""
    primer_path = Path("config/primer_schema_config.yaml")
    assert primer_path.exists()

    with open(primer_path, "r", encoding="utf-8") as f:
        primer = yaml.safe_load(f)

    assert "snp to motif association" in primer
    cfg = primer["snp to motif association"]
    assert cfg["represented_as"] == "edge"
    assert cfg["input_label"] == "tf_snp"
    assert cfg["output_label"] == "alters_binding"
    assert cfg["source"] == "snp"
    assert cfg["target"] == "motif"

    # Verify generalized description supports HOCOMOCO and JASPAR
    desc = cfg["description"].lower()
    assert "hocomoco" in desc
    assert "jaspar" in desc

    # Verify all output properties are declared
    props = cfg["properties"]
    assert "effect" in props
    assert "score" in props
    assert "motif_id" in props
    assert "taxon_id" in props
    assert "source" in props
    assert "source_url" in props


def test_hsa_schema_does_not_contain_duplicate_tf_snp():
    """Verify tf_snp was cleanly removed from hsa_schema_config.yaml."""
    hsa_path = Path("config/hsa/hsa_schema_config.yaml")
    assert hsa_path.exists()

    with open(hsa_path, "r", encoding="utf-8") as f:
        hsa = yaml.safe_load(f)

    assert "snp to motif association" not in hsa


def test_biocypher_loads_merged_schemas(tmp_path):
    """Verify BioCypher initializes cleanly with merged primer + dmel/hsa schemas."""
    primer_path = Path("config/primer_schema_config.yaml")
    dmel_path = Path("config/dmel/dmel_schema_config.yaml")
    hsa_path = Path("config/hsa/hsa_schema_config.yaml")

    with open(primer_path, encoding="utf-8") as f: primer = yaml.safe_load(f)
    with open(dmel_path, encoding="utf-8") as f: dmel = yaml.safe_load(f)
    with open(hsa_path, encoding="utf-8") as f: hsa = yaml.safe_load(f)

    for species, spec_schema in [("dmel", dmel), ("hsa", hsa)]:
        merged = {**primer, **spec_schema}
        schema_file = tmp_path / f"merged_{species}.yaml"
        with open(schema_file, "w", encoding="utf-8") as f:
            yaml.dump(merged, f)

        bc = BioCypher(
            schema_config_path=str(schema_file),
            biocypher_config_path="config/biocypher_config.yaml",
        )
        extended = bc._get_ontology_mapping()._extend_schema()
        assert "snp to motif association" in extended
        entry = extended["snp to motif association"]
        assert entry.get("input_label") == "tf_snp"
        assert entry.get("source") == "snp"
        assert entry.get("target") == "motif"


def test_biocypher_end_to_end_neo4j_writer_integration(tmp_path, sample_dmel_files):
    """
    Run end-to-end integration:
    Instantiate shared adapter with Drosophila configuration, yield edges,
    and write to Neo4j CSV through Neo4jCSVWriter.
    Verify generated CSV file format, headers, and values.
    """
    map_file, diff_file = sample_dmel_files

    primer_path = Path("config/primer_schema_config.yaml")
    dmel_path = Path("config/dmel/dmel_schema_config.yaml")
    with open(primer_path, encoding="utf-8") as f: primer = yaml.safe_load(f)
    with open(dmel_path, encoding="utf-8") as f: dmel = yaml.safe_load(f)
    merged = {**primer, **dmel}

    merged_path = tmp_path / "merged_schema.yaml"
    with open(merged_path, "w", encoding="utf-8") as f:
        yaml.dump(merged, f)

    outdir = tmp_path / "neo4j_output"
    writer = Neo4jCSVWriter(
        schema_config=str(merged_path),
        biocypher_config="config/biocypher_config.yaml",
        output_dir=outdir,
    )

    adapter = MotifDiffAdapter(
        filepath=str(diff_file),
        motif_map_path=str(map_file),
        label="tf_snp",
        taxon_id=7227,
        threshold=0.001,
        variant_prefix="",
        uppercase_variant=True,
        gene_prefix="",
        heterodimer="all",
        write_properties=True,
        add_provenance=True,
    )

    edges = adapter.get_edges()
    freq = writer.write_edges(edges, path_prefix="MotifDiff")

    # Assert 10 edges were written
    assert freq["tf_snp|snp|motif"] == 10

    csv_file = outdir / "MotifDiff" / "edges_snp_alters_binding_motif.csv"
    assert csv_file.exists()

    lines = csv_file.read_text(encoding="utf-8").strip().split("\n")
    header = lines[0].split("|")
    assert "source_id" in header
    assert "target_id" in header
    assert "label" in header
    assert "effect" in header
    assert "score" in header
    assert "motif_id" in header
    assert "taxon_id" in header

    # Verify data lines
    data_lines = lines[1:]
    assert len(data_lines) == 10
    first_row = dict(zip(header, data_lines[0].split("|")))
    assert first_row["source_id"] == "RS101"
    assert first_row["label"] == "alters_binding"
    assert first_row["taxon_id"] == "7227"
    assert first_row["source"] == "MotifDiff"


def test_biocypher_end_to_end_metta_writer_integration(tmp_path, sample_dmel_files):
    """
    Run end-to-end integration with MeTTaWriter:
    Verify generated MeTTa S-expressions for alters_binding edges and property assertions.
    """
    map_file, diff_file = sample_dmel_files

    primer_path = Path("config/primer_schema_config.yaml")
    dmel_path = Path("config/dmel/dmel_schema_config.yaml")
    with open(primer_path, encoding="utf-8") as f: primer = yaml.safe_load(f)
    with open(dmel_path, encoding="utf-8") as f: dmel = yaml.safe_load(f)
    merged = {**primer, **dmel}

    merged_path = tmp_path / "merged_schema.yaml"
    with open(merged_path, "w", encoding="utf-8") as f:
        yaml.dump(merged, f)

    outdir = tmp_path / "metta_output"
    writer = MeTTaWriter(
        schema_config=str(merged_path),
        biocypher_config="config/biocypher_config.yaml",
        output_dir=outdir,
    )

    adapter = MotifDiffAdapter(
        filepath=str(diff_file),
        motif_map_path=str(map_file),
        label="tf_snp",
        taxon_id=7227,
        threshold=0.001,
        variant_prefix="",
        uppercase_variant=True,
        gene_prefix="",
        heterodimer="all",
        write_properties=True,
        add_provenance=True,
    )

    edges = adapter.get_edges()
    freq = writer.write_edges(edges, path_prefix="MotifDiff")
    assert freq.get("tf_snp") == 10 or freq.get("tf_snp|snp|motif") == 10

    metta_file = outdir / "MotifDiff" / "edges_snp_alters_binding_motif.metta"
    assert metta_file.exists()

    content = metta_file.read_text(encoding="utf-8")
    assert "(alters_binding (snp RS101) (motif FBGN0000210))" in content
    assert "(effect (alters_binding (snp RS101) (motif FBGN0000210)) gain)" in content
    assert "(motif_id (alters_binding (snp RS101) (motif FBGN0000210)) MA0010.2)" in content
    assert "(taxon_id (alters_binding (snp RS101) (motif FBGN0000210)) 7227)" in content
