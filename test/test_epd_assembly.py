import pytest
from pathlib import Path
from biocypher_metta.adapters.helpers import (
    ALLOWED_ASSEMBLIES,
    build_regulatory_region_id,
    build_variant_id,
)
from biocypher_metta.adapters.epd_adapter import EPDAdapter


ROOT_DIR = Path(__file__).resolve().parent.parent
SAMPLES_DIR = ROOT_DIR / "samples"


class MockHGNCProcessor:
    """Mock processor for tests to avoid external network/data calls."""
    def load_or_update(self):
        pass

    def get_ensembl_id(self, symbol):
        return f"ENSG_{symbol}"


class MockGeneSymbolProcessor:
    """Mock processor for non-human species symbol resolution."""
    def load_or_update(self):
        pass

    def get_ensembl_gene(self, symbol):
        return f"ENS_MOCK_{symbol}"


def test_allowed_assemblies_contains_all_species():
    """Verify that ALLOWED_ASSEMBLIES includes GRCh38 and all supported non-human assemblies."""
    expected = {'GRCh38', 'dm6', 'ce6', 'mm10', 'rn6'}
    assert expected.issubset(set(ALLOWED_ASSEMBLIES))
    assert 'GRCh38' in ALLOWED_ASSEMBLIES


def test_build_regulatory_region_id_defaults_to_grch38():
    """Verify default backward-compatible behavior when assembly is omitted."""
    region_id = build_regulatory_region_id("chr1", 1000, 2000)
    assert region_id == "chr1_1000_2000_GRCh38"


def test_build_regulatory_region_id_supported_assemblies():
    """Verify build_regulatory_region_id with each supported assembly."""
    for asm in ['GRCh38', 'dm6', 'ce6', 'mm10', 'rn6']:
        region_id = build_regulatory_region_id("chr1", 1000, 2000, assembly=asm)
        assert region_id == f"chr1_1000_2000_{asm}"


def test_build_regulatory_region_id_unsupported_assembly():
    """Verify that unsupported assemblies raise ValueError('Assembly not supported')."""
    with pytest.raises(ValueError, match="Assembly not supported"):
        build_regulatory_region_id("chr1", 1000, 2000, assembly="hg19")

    with pytest.raises(ValueError, match="Assembly not supported"):
        build_regulatory_region_id("chr1", 1000, 2000, "unknown_assembly")


def test_build_variant_id_unsupported_assembly():
    """Verify that build_variant_id enforces assembly_check decorator validation."""
    with pytest.raises(ValueError, match="Assembly not supported"):
        build_variant_id("chr1", 1000, "A", "T", assembly="invalid_assembly")


def test_human_epd_backward_compatibility():
    """
    Verify that:
    1. Human EPD IDs remain identical before and after the change.
    2. Node and edge outputs are identical between assembly=None and assembly='GRCh38'.
    """
    human_bed = SAMPLES_DIR / "hsa" / "epd" / "Hs_EPDnew.bed.gz"
    assert human_bed.exists(), f"Sample file not found: {human_bed}"

    mock_hgnc = MockHGNCProcessor()

    # Default adapter (assembly=None)
    adapter_default = EPDAdapter(
        filepath=str(human_bed),
        label="promoter",
        taxon_id=9606,
        write_properties=True,
        add_provenance=True,
        hgnc_processor=mock_hgnc,
    )

    # Explicit adapter (assembly='GRCh38')
    adapter_explicit = EPDAdapter(
        filepath=str(human_bed),
        label="promoter",
        taxon_id=9606,
        write_properties=True,
        add_provenance=True,
        assembly="GRCh38",
        hgnc_processor=mock_hgnc,
    )

    nodes_default = list(adapter_default.get_nodes())
    nodes_explicit = list(adapter_explicit.get_nodes())

    assert len(nodes_default) > 0
    assert len(nodes_default) == len(nodes_explicit)

    # Every single node ID, label, and property dict must match 100%
    for (id_def, lbl_def, props_def), (id_exp, lbl_exp, props_exp) in zip(nodes_default, nodes_explicit):
        assert id_def == id_exp
        assert id_def.endswith("_GRCh38")
        assert lbl_def == lbl_exp
        assert props_def == props_exp

    # Now verify edges
    adapter_edges_default = EPDAdapter(
        filepath=str(human_bed),
        label="promoter_gene",
        type="promoter to gene association",
        taxon_id=9606,
        write_properties=True,
        add_provenance=True,
        hgnc_processor=mock_hgnc,
    )

    adapter_edges_explicit = EPDAdapter(
        filepath=str(human_bed),
        label="promoter_gene",
        type="promoter to gene association",
        taxon_id=9606,
        write_properties=True,
        add_provenance=True,
        assembly="GRCh38",
        hgnc_processor=mock_hgnc,
    )

    edges_default = list(adapter_edges_default.get_edges())
    edges_explicit = list(adapter_edges_explicit.get_edges())

    assert len(edges_default) > 0
    assert len(edges_default) == len(edges_explicit)

    for (src_def, dst_def, lbl_def, props_def), (src_exp, dst_exp, lbl_exp, props_exp) in zip(edges_default, edges_explicit):
        assert src_def == src_exp
        assert src_def.endswith("_GRCh38")
        assert dst_def == dst_exp
        assert lbl_def == lbl_exp
        assert props_def == props_exp


@pytest.mark.parametrize(
    "species,sample_file,assembly,taxon_id",
    [
        ("dmel", "Dm_EPDnew.bed.gz", "dm6", 7227),
        ("cel", "Ce_EPDnew.bed.gz", "ce6", 6239),
        ("mmu", "Mm_EPDnew.bed.gz", "mm10", 10090),
        ("rno", "Rn_EPDnew.bed.gz", "rn6", 10116),
    ]
)
def test_non_human_epd_uses_configured_assembly(species, sample_file, assembly, taxon_id):
    """
    Verify that non-human EPD nodes and edges:
    1. Contain the configured species assembly in their IDs.
    2. Never receive a human assembly tag (GRCh38).
    """
    bed_path = SAMPLES_DIR / species / "epd" / sample_file
    assert bed_path.exists(), f"Sample file not found: {bed_path}"

    mock_resolver = MockGeneSymbolProcessor()

    adapter = EPDAdapter(
        filepath=str(bed_path),
        label="promoter",
        taxon_id=taxon_id,
        assembly=assembly,
        write_properties=True,
        gene_symbol_processor=mock_resolver,
    )

    nodes = list(adapter.get_nodes())
    assert len(nodes) > 0

    for node_id, label, props in nodes:
        assert node_id.startswith("EPD:")
        assert node_id.endswith(f"_{assembly}")
        assert "_GRCh38" not in node_id
        assert "hg38" not in node_id
        assert label == "promoter"

    # Edge test
    adapter_edges = EPDAdapter(
        filepath=str(bed_path),
        label="promoter_gene",
        type="promoter to gene association",
        taxon_id=taxon_id,
        assembly=assembly,
        write_properties=True,
        gene_symbol_processor=mock_resolver,
    )

    edges = list(adapter_edges.get_edges())
    assert len(edges) > 0

    for src_id, dst_id, label, props in edges:
        assert src_id.startswith("EPD:")
        assert src_id.endswith(f"_{assembly}")
        assert "_GRCh38" not in src_id
        assert "hg38" not in src_id


def test_epd_adapter_unsupported_assembly_raises():
    """Verify that instantiating EPDAdapter with an unsupported assembly raises ValueError during execution."""
    human_bed = SAMPLES_DIR / "hsa" / "epd" / "Hs_EPDnew.bed.gz"
    adapter = EPDAdapter(
        filepath=str(human_bed),
        label="promoter",
        assembly="unsupported_assembly_123",
    )
    with pytest.raises(ValueError, match="Assembly not supported"):
        next(adapter.get_nodes())

