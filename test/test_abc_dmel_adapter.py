import csv
import gzip
from pathlib import Path
import pickle
import pytest
import tempfile

from biocypher_metta.adapters.dmel.abc_dmel_adapter import (
    ABCDmelAdapter,
    DEFAULT_TISSUE_TO_ONTOLOGY_ID_MAP,
)
from biocypher_metta.adapters.hsa.abc_adapter import ABCAdapter


@pytest.fixture
def sample_csv_data():
    headers = [
        'rsid', 'chromosome', 'position', 'element_chr', 'element_start',
        'element_end', 'element_name', 'target_gene', 'target_gene_tss',
        'abc_score', 'cell_type'
    ]
    rows = [
        ['rs1001', '2L', '5000', '2L', '4800', '5200', 'genic|2L:4800-5200', 'FBgn0000042', '5100', '0.045', 'Kc167'],
        ['RS1002', '2L', '6000', '2L', '5800', '6200', 'intergenic|2L:5800-6200', 'FBgn0001091', '6050', '0.080', 'S2'],
        ['rs1003', '3R', '12000', '3R', '11800', '12200', 'genic|3R:11800-12200', 'FBgn0002626', '12050', '0.025', 'Kc167, S2'],
        ['rs1004', '2L', '7000', '2L', '6800', '7200', 'genic|2L:6800-7200', 'FBgn0003462', '7100', '0.010', 'Kc167'],
    ]
    return headers, rows


@pytest.fixture
def sample_csv_file(tmp_path, sample_csv_data):
    headers, rows = sample_csv_data
    file_path = tmp_path / 'abc_test.csv'
    with open(file_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    return file_path


@pytest.fixture
def sample_gz_file(tmp_path, sample_csv_data):
    headers, rows = sample_csv_data
    file_path = tmp_path / 'abc_test.csv.gz'
    with gzip.open(file_path, 'wt', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    return file_path


def test_import_paths():
    """Verify both canonical and re-export import paths work."""
    from biocypher_metta.adapters.dmel.abc_dmel_adapter import ABCDmelAdapter as Cls1
    from biocypher_metta.adapters.abc_dmel_adapter import ABCDmelAdapter as Cls2
    from biocypher_metta.adapters.dmel import ABCDmelAdapter as Cls3
    assert Cls1 is Cls2
    assert Cls1 is Cls3


def test_valid_edge_generation(sample_gz_file):
    """Test standard edge generation and tuple structure."""
    adapter = ABCDmelAdapter(filepath=sample_gz_file)
    edges = list(adapter.get_edges())
    assert len(edges) == 4

    src, tgt, label, props = edges[0]
    assert src == 'RS1001'
    assert tgt == 'FBgn0000042'
    assert label == 'activity_by_contact'
    assert props['score'] == pytest.approx(0.045)
    assert props['biological_context'] == 'FBbt:00005349'
    assert props['taxon_id'] == 7227


def test_uppercase_variant_ids(sample_csv_file):
    """Verify rsIDs are uppercased by default."""
    adapter = ABCDmelAdapter(filepath=sample_csv_file, variant_case='upper')
    edges = list(adapter.get_edges())
    for src, _, _, _ in edges:
        assert src.startswith('RS')
        assert src.isupper()


def test_target_gene_bare_id(sample_csv_file):
    """Verify gene IDs remain in bare FBgn format."""
    adapter = ABCDmelAdapter(filepath=sample_csv_file, gene_prefix='')
    edges = list(adapter.get_edges())
    for _, tgt, _, _ in edges:
        assert tgt.startswith('FBgn')
        assert not tgt.startswith('FlyBase:FBgn')


def test_kc167_and_s2_mappings(sample_csv_file):
    """Verify Kc167 and S2 map to expected FBbt identifiers."""
    adapter = ABCDmelAdapter(filepath=sample_csv_file)
    edges = list(adapter.get_edges())

    # row 0 is Kc167 -> FBbt:00005349
    assert edges[0][3]['biological_context'] == 'FBbt:00005349'
    # row 1 is S2 -> FBbt:00005737
    assert edges[1][3]['biological_context'] == 'FBbt:00005737'
    # row 2 is Kc167, S2 -> FBbt:00005349,FBbt:00005737
    assert 'FBbt:00005349' in edges[2][3]['biological_context']
    assert 'FBbt:00005737' in edges[2][3]['biological_context']


def test_custom_ontology_map_dict(sample_csv_file):
    """Test custom tissue map provided as a dictionary."""
    custom_map = {'Kc167': 'CUSTOM:001', 'S2': 'CUSTOM:002'}
    adapter = ABCDmelAdapter(filepath=sample_csv_file, tissue_to_ontology_id_map=custom_map)
    edges = list(adapter.get_edges())
    assert edges[0][3]['biological_context'] == 'CUSTOM:001'
    assert edges[1][3]['biological_context'] == 'CUSTOM:002'


def test_custom_ontology_map_pickle(tmp_path, sample_csv_file):
    """Test custom tissue map loaded from a pickle file."""
    custom_map = {'Kc167': 'PICKLE:001', 'S2': 'PICKLE:002'}
    pickle_path = tmp_path / 'custom_map.pkl'
    with open(pickle_path, 'wb') as f:
        pickle.dump(custom_map, f)

    adapter = ABCDmelAdapter(filepath=sample_csv_file, tissue_to_ontology_id_map=pickle_path)
    edges = list(adapter.get_edges())
    assert edges[0][3]['biological_context'] == 'PICKLE:001'
    assert edges[1][3]['biological_context'] == 'PICKLE:002'


def test_single_and_multiple_filepaths(tmp_path, sample_csv_data):
    """Test adapter handles both a single path and a list of paths."""
    headers, rows = sample_csv_data
    file1 = tmp_path / 'f1.csv'
    file2 = tmp_path / 'f2.csv'
    for fpath in (file1, file2):
        with open(fpath, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(headers)
            w.writerows(rows[:2])

    adapter_single = ABCDmelAdapter(filepath=file1)
    assert len(list(adapter_single.get_edges())) == 2

    adapter_multi = ABCDmelAdapter(filepath=[file1, file2])
    assert len(list(adapter_multi.get_edges())) == 4


def test_threshold_none_keeps_all(sample_csv_file):
    """Test that threshold=None does not filter any rows."""
    adapter = ABCDmelAdapter(filepath=sample_csv_file, threshold=None)
    assert len(list(adapter.get_edges())) == 4


def test_numeric_threshold_filtering(sample_csv_file):
    """Test score threshold filtering and boundary conditions."""
    # Rows have scores: 0.045, 0.080, 0.025, 0.010
    adapter = ABCDmelAdapter(filepath=sample_csv_file, threshold=0.03)
    edges = list(adapter.get_edges())
    assert len(edges) == 2  # 0.045 and 0.080 survive

    # Boundary test: exact match should pass (>= threshold)
    adapter_boundary = ABCDmelAdapter(filepath=sample_csv_file, threshold=0.025)
    edges_boundary = list(adapter_boundary.get_edges())
    assert len(edges_boundary) == 3  # 0.045, 0.080, and 0.025 survive


def test_genomic_location_filter(sample_csv_file):
    """Test chromosome, start, and end filtering."""
    # Row 0: 2L:5000, Row 1: 2L:6000, Row 2: 3R:12000, Row 3: 2L:7000
    adapter_chr = ABCDmelAdapter(filepath=sample_csv_file, chr='2L')
    assert len(list(adapter_chr.get_edges())) == 3

    # Also test 'chr2L' format normalization
    adapter_chr_prefix = ABCDmelAdapter(filepath=sample_csv_file, chr='chr2L')
    assert len(list(adapter_chr_prefix.get_edges())) == 3

    # Position range
    adapter_range = ABCDmelAdapter(filepath=sample_csv_file, chr='2L', start=5500, end=6500)
    edges_range = list(adapter_range.get_edges())
    assert len(edges_range) == 1
    assert edges_range[0][0] == 'RS1002'


def test_invalid_input_handling(tmp_path):
    """Test missing or malformed rows are safely skipped."""
    headers = [
        'rsid', 'chromosome', 'position', 'element_chr', 'element_start',
        'element_end', 'element_name', 'target_gene', 'target_gene_tss',
        'abc_score', 'cell_type'
    ]
    bad_rows = [
        ['', '2L', '100', '2L', '50', '150', 'e1', 'FBgn001', '100', '0.05', 'Kc167'],  # missing rsid
        ['rs1', '2L', '100', '2L', '50', '150', 'e1', '', '100', '0.05', 'Kc167'],       # missing gene
        ['rs2', '2L', '100', '2L', '50', '150', 'e1', 'FBgn002', '100', 'not_a_num', 'Kc167'],  # bad score
        ['rs3', '2L', '100', '2L', '50', '150', 'e1', 'FBgn003', '100', '0.05', 'UnknownCell'],  # unmapped cell
        ['rs4', '2L', 'not_pos', '2L', '50', '150', 'e1', 'FBgn004', '100', '0.05', 'Kc167'],   # bad pos with chr filter
        ['rs5', '2L', '100', '2L', '50', '150', 'genic|2L:50-150', 'FBgn005', '100', '0.05', 'Kc167'],  # valid row
    ]
    test_file = tmp_path / 'bad_data.csv'
    with open(test_file, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(bad_rows)

    adapter = ABCDmelAdapter(filepath=test_file, chr='2L', start=10, end=200)
    edges = list(adapter.get_edges())
    assert len(edges) == 1
    assert edges[0][0] == 'RS5'
    assert edges[0][1] == 'FBgn005'


def test_element_name_with_pipe(sample_csv_file):
    """Ensure element_name containing pipes does not disrupt parsing."""
    adapter = ABCDmelAdapter(filepath=sample_csv_file)
    edges = list(adapter.get_edges())
    # All 4 rows in sample fixture have element_name containing '|'
    assert len(edges) == 4


def test_generated_sample_file_integration():
    """Verify the real generated sample file can be read by ABCDmelAdapter."""
    sample_path = Path('samples/dmel/abc/abc_dmel_variants_sample.csv.gz')
    assert sample_path.exists(), f"Sample file not found at {sample_path}"

    adapter = ABCDmelAdapter(filepath=sample_path)
    edges = list(adapter.get_edges())
    assert len(edges) == 16
    for src, tgt, label, props in edges:
        assert src.startswith('RS')
        assert tgt.startswith('FBgn')
        assert label == 'activity_by_contact'
        assert isinstance(props['score'], float)
        assert 'FBbt:' in props['biological_context']
        assert props['taxon_id'] == 7227


def test_human_abc_adapter_regression():
    """Regression test ensuring human ABCAdapter is preserved and functional."""
    # Verify signature and attributes
    import inspect
    sig = inspect.signature(ABCAdapter.__init__)
    assert 'filepath' in sig.parameters
    assert 'hgnc_to_ensembl_map' in sig.parameters
    assert 'tissue_to_ontology_id_map' in sig.parameters
    assert 'dbsnp_rsid_map' in sig.parameters
    assert 'taxon_id' in sig.parameters
    assert sig.parameters['taxon_id'].default == 9606


def test_schema_compatibility():
    """Verify Dmel schema config has activity by contact edge between snp and gene."""
    import yaml
    schema_path = Path("config/dmel/dmel_schema_config.yaml")
    assert schema_path.exists()

    with open(schema_path, "r") as f:
        schema = yaml.safe_load(f)

    assert "activity by contact" in schema
    abc_cfg = schema["activity by contact"]
    assert abc_cfg["represented_as"] == "edge"
    assert abc_cfg["input_label"] == "activity_by_contact"
    assert abc_cfg["source"] == "snp"
    assert abc_cfg["target"] == "gene"
    assert "score" in abc_cfg["properties"]
    assert "biological_context" in abc_cfg["properties"]
    assert "taxon_id" in abc_cfg["properties"]

    score_type = abc_cfg["properties"]["score"]
    assert (score_type["type"] if isinstance(score_type, dict) else score_type) == "float"

    ctx_type = abc_cfg["properties"]["biological_context"]
    assert (ctx_type["type"] if isinstance(ctx_type, dict) else ctx_type) == "str"

    tax_type = abc_cfg["properties"]["taxon_id"]
    assert (tax_type["type"] if isinstance(tax_type, dict) else tax_type) == "int"


def test_biocypher_schema_initialization(tmp_path):
    """Verify BioCypher initializes cleanly with the merged Dmel schema including activity_by_contact."""
    import yaml
    from biocypher import BioCypher

    primer_path = Path("config/primer_schema_config.yaml")
    dmel_path = Path("config/dmel/dmel_schema_config.yaml")
    assert primer_path.exists()
    assert dmel_path.exists()

    with open(primer_path, "r") as f:
        primer = yaml.safe_load(f)
    with open(dmel_path, "r") as f:
        dmel = yaml.safe_load(f)

    merged = {**primer, **dmel}
    temp_schema = tmp_path / "merged_schema.yaml"
    with open(temp_schema, "w") as f:
        yaml.dump(merged, f)

    bc = BioCypher(
        schema_config_path=str(temp_schema),
        biocypher_config_path="config/biocypher_config.yaml",
    )
    schema = bc._get_ontology_mapping()._extend_schema()
    assert "activity by contact" in schema
    entry = schema["activity by contact"]
    assert entry.get("source") == "snp"
    assert entry.get("target") == "gene"
