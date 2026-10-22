#!/usr/bin/env python3
"""
Automated Data Source Schema Generator v4

Generates data source schema YAML files by analyzing:
1. Adapter implementations (to extract properties)
2. Schema configuration (to get descriptions and validate properties)
3. Adapter configuration (to map adapters to data sources)

Key features:
- Groups by actual adapter provenance (self.source)
- Uses specific source_url for each node/relationship
- Handles all property assignment patterns (props={}, props[], props.update())
- One schema file per unique data source
"""

import argparse
import ast
import sys
import warnings
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Any, Optional
import urllib.parse

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from config.yaml_loader import load_yaml_with_includes
except ImportError as exc:  # pragma: no cover
    warnings.warn(
        f"Falling back to yaml.safe_load because config.yaml_loader could not be imported: {exc}",
        RuntimeWarning,
        stacklevel=2,
    )
    load_yaml_with_includes = yaml.safe_load


class FlowStyleListDumper(yaml.SafeDumper):
    """Custom YAML dumper that outputs lists in flow style (rectangular brackets)."""
    pass


def represent_list(dumper, data):
    """Represent lists in flow style for source/target fields."""
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


FlowStyleListDumper.add_representer(list, represent_list)


class AdapterAnalyzer:
    def __init__(self, adapter_path: Path):
        self.adapter_path = adapter_path
        self.source_code = adapter_path.read_text()
        self.tree = ast.parse(self.source_code)
        self.class_attributes = self._extract_class_attributes()
        self.label_attributes = self._extract_label_attributes()
        # metadata attr -> (class dict attribute name, constructor arg name used
        # as the subscript key), for self.attr = ClassName.DICT[ctor_arg] patterns
        # that a specific adapter config block's args can resolve precisely
        # (see get_metadata_from_init's Subscript case and its caller).
        self.dict_lookups: Dict[str, tuple] = {}

    @staticmethod
    def _is_label_name(name: str) -> bool:
        label_names = {'label', 'input_label', 'node_label', 'edge_label', 'relationship_label'}
        return name in label_names or name.endswith('_label')

    @classmethod
    def _extract_string_values(cls, node: ast.AST) -> Set[str]:
        values = set()
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            values.add(node.value)
        elif isinstance(node, ast.BoolOp):
            for value in node.values:
                values.update(cls._extract_string_values(value))
        elif isinstance(node, ast.IfExp):
            values.update(cls._extract_string_values(node.body))
            values.update(cls._extract_string_values(node.orelse))
        return values

    def _extract_class_attributes(self) -> Dict[str, Any]:
        class_attrs = {}
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                for item in node.body:
                    if isinstance(item, ast.Assign):
                        for target in item.targets:
                            if isinstance(target, ast.Name):
                                if isinstance(item.value, ast.Dict):
                                    dict_value = {}
                                    for key, value in zip(item.value.keys, item.value.values):
                                        if not isinstance(value, ast.Constant):
                                            continue
                                        if isinstance(key, ast.Constant):
                                            dict_value[key.value] = value.value
                                        elif isinstance(key, ast.Tuple):
                                            key_parts = []
                                            for elt in key.elts:
                                                if isinstance(elt, ast.Constant):
                                                    key_parts.append(elt.value)
                                            if len(key_parts) == len(key.elts):
                                                dict_value[tuple(key_parts)] = value.value
                                    if dict_value:
                                        class_attrs[target.id] = dict_value
                                elif isinstance(item.value, ast.Constant):
                                    class_attrs[target.id] = item.value.value
        return class_attrs

    def get_adapter_class_name(self) -> str:
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                return node.name
        return self.adapter_path.stem

    def _extract_label_attributes(self) -> Dict[str, Set[str]]:
        label_attrs = defaultdict(set)
        for node in ast.walk(self.tree):
            if not isinstance(node, ast.Assign):
                continue
            values = self._extract_string_values(node.value)
            if not values:
                continue
            for target in node.targets:
                target_name = None
                if isinstance(target, ast.Name):
                    target_name = target.id
                elif isinstance(target, ast.Attribute):
                    target_name = target.attr
                if target_name and self._is_label_name(target_name):
                    label_attrs[target_name].update(values)
        return dict(label_attrs)

    def get_class_dict_string_values(self) -> Set[str]:
        """Return string values from class-level dictionaries used as label maps."""
        values = set()
        for attr_value in self.class_attributes.values():
            if isinstance(attr_value, dict):
                for value in attr_value.values():
                    if isinstance(value, str):
                        values.add(value)
        return values

    def get_yield_string_labels(self, method_name: str, tuple_index: int) -> Set[str]:
        """Extract constant string labels yielded at a tuple position in a method."""
        labels = set()
        for node in ast.walk(self.tree):
            if isinstance(node, ast.FunctionDef) and node.name == method_name:
                for stmt in ast.walk(node):
                    if isinstance(stmt, (ast.Yield, ast.YieldFrom)):
                        value = stmt.value
                        if isinstance(value, ast.Tuple) and len(value.elts) > tuple_index:
                            label_node = value.elts[tuple_index]
                            if isinstance(label_node, ast.Constant) and isinstance(label_node.value, str):
                                labels.add(label_node.value)
                            elif isinstance(label_node, ast.Attribute):
                                labels.update(self.label_attributes.get(label_node.attr, set()))
        return labels

    def get_label_literals(self) -> Set[str]:
        """Return strings assigned to variables or attributes that are explicitly label-like."""
        labels = set()
        for node in ast.walk(self.tree):
            if not isinstance(node, ast.Assign):
                continue
            values = self._extract_string_values(node.value)
            if not values:
                continue
            for target in node.targets:
                target_name = None
                if isinstance(target, ast.Name):
                    target_name = target.id
                elif isinstance(target, ast.Attribute):
                    target_name = target.attr
                if target_name and self._is_label_name(target_name):
                    labels.update(values)
        return labels

    def get_feature_overlap_labels(self, feature_labels: List[str], variant_label: str) -> List[str]:
        """Resolve labels generated from feature overlap templates in the adapter."""
        labels = []
        source = self.source_code
        has_feature_to_variant = "_overlaps_" in source and f"_overlaps_{variant_label}" in source
        has_variant_to_feature = "_overlaps_{" in source or f"{variant_label}_overlaps_" in source

        for feature_label in feature_labels:
            if has_feature_to_variant:
                labels.append(f"{feature_label}_overlaps_{variant_label}")
            if has_variant_to_feature:
                labels.append(f"{variant_label}_overlaps_{feature_label}")
        return labels

    def inherits_from_ontology_adapter(self) -> bool:
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                for base in node.bases:
                    if isinstance(base, ast.Name) and base.id == 'OntologyAdapter':
                        return True
                    if isinstance(base, ast.Attribute) and base.attr == 'OntologyAdapter':
                        return True
        return False

    def get_metadata_from_init(self) -> Dict[str, Any]:
        metadata = {}
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                for item in node.body:
                    if isinstance(item, ast.FunctionDef) and item.name == '__init__':
                        # Extract default values from function signature
                        if item.args.defaults:
                            arg_names = [arg.arg for arg in item.args.args]
                            num_defaults = len(item.args.defaults)
                            default_start_idx = len(arg_names) - num_defaults
                            for i, default in enumerate(item.args.defaults):
                                arg_idx = default_start_idx + i
                                arg_name = arg_names[arg_idx]
                                if arg_name in ['label', 'version'] and isinstance(default, ast.Constant):
                                    if arg_name not in metadata:
                                        metadata[arg_name] = default.value

                        # Extract from assignments (these override defaults)
                        for stmt in ast.walk(item):
                            if isinstance(stmt, ast.Assign):
                                for target in stmt.targets:
                                    if isinstance(target, ast.Attribute):
                                        if target.attr in ['source', 'source_url', 'version', 'label']:
                                            if isinstance(stmt.value, ast.Constant):
                                                metadata[target.attr] = stmt.value.value
                                            elif isinstance(stmt.value, ast.JoinedStr):
                                                url_parts = []
                                                for value in stmt.value.values:
                                                    if isinstance(value, ast.Constant):
                                                        url_parts.append(str(value.value))
                                                if url_parts:
                                                    metadata[target.attr] = ''.join(url_parts)
                                            elif isinstance(stmt.value, ast.Subscript):
                                                if isinstance(stmt.value.value, ast.Attribute):
                                                    dict_name = stmt.value.value.attr
                                                    if dict_name in self.class_attributes:
                                                        dict_value = self.class_attributes[dict_name]
                                                        if isinstance(dict_value, dict) and dict_value:
                                                            # Naive fallback (first value) -- may not match this
                                                            # specific adapter config block's actual constructor
                                                            # arg. Record the (dict, index-var) pair too, so the
                                                            # caller can re-resolve precisely once it knows this
                                                            # block's real args (see self.dict_lookups above).
                                                            metadata[target.attr] = next(iter(dict_value.values()))
                                                            index = stmt.value.slice
                                                            if isinstance(index, ast.Name):
                                                                self.dict_lookups[target.attr] = (dict_name, index.id)
                                            elif isinstance(stmt.value, ast.Call):
                                                # e.g. self.source_url = _PER_SPECIES_URLS.get(taxon_id, 'https://example.org/')
                                                # -- use the literal fallback as the representative value.
                                                func = stmt.value.func
                                                if (isinstance(func, ast.Attribute) and func.attr == 'get'
                                                        and len(stmt.value.args) >= 2
                                                        and isinstance(stmt.value.args[1], ast.Constant)):
                                                    metadata[target.attr] = stmt.value.args[1].value
        return metadata

    def extract_properties_from_dict(self, dict_node: ast.Dict) -> Set[str]:
        properties = set()
        for key in dict_node.keys:
            if isinstance(key, ast.Constant) and isinstance(key.value, str):
                if key.value not in ['source', 'source_url']:
                    properties.add(key.value)
        return properties

    @staticmethod
    def _is_property_var_name(name: str) -> bool:
        return name in ['props', '_props', 'properties'] or name.endswith('_props')

    def get_property_variables_from_method(self, method_node: ast.FunctionDef, tuple_index: int) -> Set[str]:
        property_vars = set()
        for stmt in ast.walk(method_node):
            if isinstance(stmt, ast.Assign):
                for target in stmt.targets:
                    if isinstance(target, ast.Name) and self._is_property_var_name(target.id):
                        property_vars.add(target.id)

            if isinstance(stmt, (ast.Yield, ast.YieldFrom)):
                value = stmt.value
                if isinstance(value, ast.Tuple) and len(value.elts) > tuple_index:
                    props_node = value.elts[tuple_index]
                    if isinstance(props_node, ast.Name):
                        property_vars.add(props_node.id)

        return property_vars

    def extract_properties_for_vars_from_method(
        self,
        method_node: ast.FunctionDef,
        property_vars: Set[str],
    ) -> Set[str]:
        properties = set()
        for stmt in ast.walk(method_node):
            # Pattern 1: props = {...}
            if isinstance(stmt, ast.Assign):
                for target in stmt.targets:
                    if isinstance(target, ast.Name) and target.id in property_vars:
                        if isinstance(stmt.value, ast.Dict):
                            properties.update(self.extract_properties_from_dict(stmt.value))

            # Pattern 2: props['key'] = value
            if isinstance(stmt, ast.Assign):
                for target in stmt.targets:
                    if isinstance(target, ast.Subscript):
                        if isinstance(target.value, ast.Name) and target.value.id in property_vars:
                            if isinstance(target.slice, ast.Constant):
                                prop_name = target.slice.value
                                if prop_name not in ['source', 'source_url']:
                                    properties.add(prop_name)

            # Pattern 3: props.update({...})
            if isinstance(stmt, ast.Expr):
                if isinstance(stmt.value, ast.Call):
                    if isinstance(stmt.value.func, ast.Attribute):
                        if (stmt.value.func.attr == 'update' and
                            isinstance(stmt.value.func.value, ast.Name) and
                            stmt.value.func.value.id in property_vars):
                            if stmt.value.args and isinstance(stmt.value.args[0], ast.Dict):
                                properties.update(self.extract_properties_from_dict(stmt.value.args[0]))
        return properties

    def get_properties_from_method(self, method_name: str, tuple_index: int) -> Set[str]:
        properties = set()
        for node in ast.walk(self.tree):
            if isinstance(node, ast.FunctionDef) and node.name == method_name:
                property_vars = self.get_property_variables_from_method(node, tuple_index)
                properties.update(self.extract_properties_for_vars_from_method(node, property_vars))
        return properties

    def get_properties_for_label(
        self,
        method_name: str,
        label: str,
        label_index: int,
        props_index: int,
    ) -> Set[str]:
        properties = set()
        for node in ast.walk(self.tree):
            if not isinstance(node, ast.FunctionDef) or node.name != method_name:
                continue

            property_vars = set()
            for stmt in ast.walk(node):
                if not isinstance(stmt, (ast.Yield, ast.YieldFrom)):
                    continue
                value = stmt.value
                if not isinstance(value, ast.Tuple) or len(value.elts) <= max(label_index, props_index):
                    continue
                label_node = value.elts[label_index]
                props_node = value.elts[props_index]
                if (
                    isinstance(label_node, ast.Constant)
                    and label_node.value == label
                    and isinstance(props_node, ast.Name)
                ):
                    property_vars.add(props_node.id)

            if property_vars:
                properties.update(self.extract_properties_for_vars_from_method(node, property_vars))

        return properties

    def get_all_class_properties(self) -> Set[str]:
        """Extract properties from all methods in the class (including helper methods)."""
        properties = set()
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                for item in node.body:
                    if isinstance(item, ast.FunctionDef):
                        for stmt in ast.walk(item):
                            # Pattern 1: props = {...}
                            if isinstance(stmt, ast.Assign):
                                for target in stmt.targets:
                                    if isinstance(target, ast.Name) and self._is_property_var_name(target.id):
                                        if isinstance(stmt.value, ast.Dict):
                                            properties.update(self.extract_properties_from_dict(stmt.value))

                            # Pattern 2: props['key'] = value
                            if isinstance(stmt, ast.Assign):
                                for target in stmt.targets:
                                    if isinstance(target, ast.Subscript):
                                        if isinstance(target.value, ast.Name) and self._is_property_var_name(target.value.id):
                                            if isinstance(target.slice, ast.Constant):
                                                prop_name = target.slice.value
                                                if prop_name not in ['source', 'source_url']:
                                                    properties.add(prop_name)

                            # Pattern 3: props.update({...})
                            if isinstance(stmt, ast.Expr):
                                if isinstance(stmt.value, ast.Call):
                                    if isinstance(stmt.value.func, ast.Attribute):
                                        if (stmt.value.func.attr == 'update' and
                                            isinstance(stmt.value.func.value, ast.Name) and
                                            self._is_property_var_name(stmt.value.func.value.id)):
                                            if stmt.value.args and isinstance(stmt.value.args[0], ast.Dict):
                                                properties.update(self.extract_properties_from_dict(stmt.value.args[0]))
        return properties

    def get_node_properties(self) -> Set[str]:
        """Get properties from get_nodes method and all helper methods."""
        main_props = self.get_properties_from_method('get_nodes', 2)
        all_props = self.get_all_class_properties()
        return main_props.union(all_props)

    def get_edge_properties(self) -> Set[str]:
        """Get properties from get_edges method and all helper methods."""
        main_props = self.get_properties_from_method('get_edges', 3)
        all_props = self.get_all_class_properties()
        return main_props.union(all_props)

    def get_parent_class_properties(self, method_name: str, parent_adapter_path: Path) -> Set[str]:
        if not parent_adapter_path.exists():
            return set()
        try:
            parent_analyzer = AdapterAnalyzer(parent_adapter_path)
            tuple_index = 2 if method_name == 'get_nodes' else 3
            return parent_analyzer.get_properties_from_method(method_name, tuple_index)
        except Exception as e:
            print(f"Warning: Could not analyze parent class {parent_adapter_path}: {e}")
            return set()


class SchemaGenerator:
    SPECIES_CONFIGS = {
        'hsa': {
            'schema_config': 'config/hsa/hsa_schema_config.yaml',
            'adapter_config': 'config/hsa/hsa_adapters_config.yaml',
            'output_dir': 'data_source_schemas/hsa',
        },
        'dmel': {
            'schema_config': 'config/dmel/dmel_schema_config.yaml',
            'adapter_config': 'config/dmel/dmel_adapters_config.yaml',
            'output_dir': 'data_source_schemas/dmel',
        },
    }
    SOURCE_FILENAME_ALIASES = {
        'HOCOMOCOv11': 'HOCOMOCO',
        'Reactome': 'REACTOME',
    }
    SOURCE_NAME_ALIASES = {
        'Reactome': 'REACTOME',
    }
    # Top-level keys in an adapters_config YAML that are metadata, not adapter entries.
    ADAPTER_CONFIG_METADATA_KEYS = {'input_dir'}

    @staticmethod
    def load_schema_config(schema_config_path: str, include_primer: bool = False) -> Dict:
        schema_config = {}
        if include_primer:
            # Resolved relative to this module's own file location (not schema_config_path
            # and not the process cwd), so this works regardless of which repo checkout or
            # working directory the caller runs from. schema_config_path's own depth isn't a
            # reliable anchor: create_knowledge_graph.py's --species mode always passes an
            # already-primer-merged temp file sitting directly under config/ (one level deep),
            # not the two-level-deep config/<species>/<species>_schema_config.yaml this used to
            # assume via `.resolve().parent.parent` -- that previously resolved to the repo root
            # instead of config/, raising FileNotFoundError for primer_schema_config.yaml.
            primer_path = Path(__file__).resolve().parent.parent / 'config' / 'primer_schema_config.yaml'
            with open(primer_path) as f:
                primer_config = load_yaml_with_includes(f) or {}
            schema_config.update(primer_config)

        with open(schema_config_path) as f:
            species_config = load_yaml_with_includes(f) or {}
        schema_config.update(species_config)
        return schema_config

    @staticmethod
    def load_commented_adapter_config(adapter_config_path: str) -> Dict:
        """Load adapter entries that are commented out as YAML blocks.

        Species-level schema generation should document all configured human/fly
        sources, including adapters that are temporarily disabled in the KG build
        config. This keeps those entries generated from adapter code/config rather
        than maintaining datasource YAMLs by hand.
        """
        inactive_adapters = {}
        current_block = []

        def flush_block():
            if not current_block:
                return
            text = '\n'.join(current_block)
            current_block.clear()
            try:
                parsed = yaml.safe_load(text)
            except yaml.YAMLError:
                return
            if not isinstance(parsed, dict):
                return
            for name, config in parsed.items():
                if (
                    isinstance(name, str)
                    and isinstance(config, dict)
                    and isinstance(config.get('adapter'), dict)
                    and ('nodes' in config or 'edges' in config)
                ):
                    inactive_adapters[name] = config

        for line in Path(adapter_config_path).read_text().splitlines():
            stripped = line.lstrip()
            if stripped.startswith('#'):
                content = stripped[1:]
                if content.startswith(' '):
                    content = content[1:]
                current_block.append(content)
            else:
                flush_block()
        flush_block()

        return inactive_adapters

    def __init__(
        self,
        schema_config_path: str,
        adapter_config_path: str,
        adapters_dir: str,
        output_dir: str,
        adapter_config_data: Optional[Dict] = None,
        schema_config_data: Optional[Dict] = None,
    ):
        self.schema_config_path = Path(schema_config_path)
        self.adapter_config_path = Path(adapter_config_path)
        self.adapters_dir = Path(adapters_dir)
        self.output_dir = Path(output_dir)

        if schema_config_data is not None:
            self.schema_config = schema_config_data
        else:
            # Always merge config/primer_schema_config.yaml, matching both --species
            # mode below and the main KG build pipeline (create_knowledge_graph.py) --
            # otherwise types defined only in the primer schema (e.g. enhancer,
            # promoter) silently resolve to "no schema config found" here.
            self.schema_config = self.load_schema_config(
                str(self.schema_config_path), include_primer=True,
            )
        if adapter_config_data is not None:
            self.adapter_config = adapter_config_data
        else:
            with open(self.adapter_config_path) as f:
                self.adapter_config = load_yaml_with_includes(f)

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.adapter_cache = {}

    def get_adapter_file_from_module(self, module: str) -> str:
        """Resolve a configured adapter module to a path relative to adapters_dir."""
        adapter_prefix = 'biocypher_metta.adapters.'
        if module.startswith(adapter_prefix):
            return module.removeprefix(adapter_prefix).replace('.', '/') + '.py'

        module_path = module.replace('.', '/') + '.py'
        if (self.adapters_dir / module_path).exists():
            return module_path

        adapter_file = module.split('.')[-1] + '.py'
        matches = sorted(self.adapters_dir.rglob(adapter_file))
        if len(matches) == 1:
            return str(matches[0].relative_to(self.adapters_dir))
        if len(matches) > 1:
            print(
                f"Warning: Adapter module '{module}' matched multiple files; "
                f"using top-level fallback '{adapter_file}'"
            )
        elif '.' in module:
            print(
                f"Warning: Adapter module '{module}' does not use the expected "
                f"'{adapter_prefix}' prefix and no matching subpackage file was found"
            )
        return adapter_file

    def get_labels_for_adapter_config(
        self,
        adapter_name: str,
        adapter_cfg: Dict,
        adapter_data: Dict,
    ) -> List[str]:
        """Resolve schema input labels represented by one adapter config entry."""
        adapter_info = adapter_cfg.get('adapter', {})
        adapter_args = adapter_info.get('args', {})
        writes_nodes = adapter_cfg.get('nodes', False)
        writes_edges = adapter_cfg.get('edges', False)
        analyzer = adapter_data['analyzer']

        candidates = []
        seen_candidates = set()

        def add_candidate(value):
            if isinstance(value, str) and value not in seen_candidates:
                candidates.append(value)
                seen_candidates.add(value)

        add_candidate(adapter_args.get('label'))

        metadata_label = adapter_data['metadata'].get('label')
        if metadata_label:
            add_candidate(metadata_label)

        for value in adapter_args.values():
            add_candidate(value)

        for value in sorted(analyzer.get_class_dict_string_values()):
            add_candidate(value)
        for value in sorted(analyzer.get_label_literals()):
            add_candidate(value)

        feature_labels = [
            feature_config.get('label')
            for feature_config in adapter_args.get('feature_files', [])
            if isinstance(feature_config, dict)
        ]
        variant_label = adapter_args.get('label')
        if writes_edges and feature_labels and variant_label:
            for value in analyzer.get_feature_overlap_labels(feature_labels, variant_label):
                add_candidate(value)

        if writes_nodes:
            for value in sorted(analyzer.get_yield_string_labels('get_nodes', 1)):
                add_candidate(value)
        if writes_edges:
            for value in sorted(analyzer.get_yield_string_labels('get_edges', 2)):
                add_candidate(value)

        add_candidate(adapter_name)

        labels = []
        seen_labels = set()
        for candidate in candidates:
            if candidate in seen_labels:
                continue

            type_info = self.get_schema_type_info(candidate)
            if type_info:
                represented_as = type_info['config'].get('represented_as')
                if writes_nodes and represented_as == 'node':
                    labels.append(candidate)
                    seen_labels.add(candidate)
                elif writes_edges and represented_as == 'edge':
                    labels.append(candidate)
                    seen_labels.add(candidate)
                continue

            # The adapter may compute its actual per-row label as f'{candidate}_<category>'
            # (e.g. one adapter emitting 'chromatin_state_anatomy', 'chromatin_state_tissue', ...
            # depending on data determined at runtime). Match those against the configured
            # base label by prefix, since a static analyzer can't resolve the runtime suffix.
            for represented_as, enabled in (('node', writes_nodes), ('edge', writes_edges)):
                if not enabled:
                    continue
                for label in self.get_schema_labels_with_prefix(candidate, represented_as):
                    if label not in seen_labels:
                        labels.append(label)
                        seen_labels.add(label)

        if not labels:
            labels.extend(self.infer_labels_from_adapter_properties(adapter_cfg, analyzer))

        return labels or [adapter_name]

    def get_schema_labels_with_prefix(self, prefix: str, represented_as: str) -> List[str]:
        """Find literal input_labels of the form '{prefix}_<suffix>' for a given represented_as kind."""
        marker = f'{prefix}_'
        matches = []
        seen = set()
        for type_name, type_config in self.schema_config.items():
            if not isinstance(type_config, dict):
                continue
            if type_config.get('represented_as') != represented_as:
                continue
            labels = type_config.get('input_label')
            label_values = labels if isinstance(labels, list) else [labels] if labels else []
            for label in label_values:
                if isinstance(label, str) and label.startswith(marker) and label not in seen:
                    seen.add(label)
                    matches.append(label)
        return matches

    def infer_labels_from_adapter_properties(self, adapter_cfg: Dict, analyzer: AdapterAnalyzer) -> List[str]:
        """Infer schema labels when config labels are absent/stale but properties are distinctive."""
        candidates = []
        modes = []
        if adapter_cfg.get('nodes', False):
            modes.append(('node', analyzer.get_node_properties()))
        if adapter_cfg.get('edges', False):
            modes.append(('edge', analyzer.get_edge_properties()))

        for represented_as, adapter_props in modes:
            if not adapter_props:
                continue
            scored = []
            for type_name, type_config in self.schema_config.items():
                if not isinstance(type_config, dict):
                    continue
                if type_config.get('represented_as') != represented_as:
                    continue
                labels = type_config.get('input_label')
                if labels is None:
                    continue
                schema_props = set(self.get_all_properties_for_type(type_name))
                overlap = adapter_props.intersection(schema_props)
                if overlap:
                    label_values = labels if isinstance(labels, list) else [labels]
                    scored.append((
                        len(overlap),
                        self.get_type_depth(type_name),
                        len(schema_props),
                        type_name,
                        label_values,
                    ))
            if not scored:
                continue
            best_score = max(score for score, _, _, _, _ in scored)
            best_depth = max(depth for score, depth, _, _, _ in scored if score == best_score)
            best_schema_size = max(
                schema_size
                for score, depth, schema_size, _, _ in scored
                if score == best_score and depth == best_depth
            )
            for score, depth, schema_size, _, label_values in sorted(scored, key=lambda item: item[3]):
                if score == best_score and depth == best_depth and schema_size == best_schema_size:
                    for label in label_values:
                        if label not in candidates:
                            candidates.append(label)

        return candidates

    def get_type_depth(self, type_name: str, visited: Set[str] = None) -> int:
        if visited is None:
            visited = set()
        if type_name in visited:
            return 0
        visited.add(type_name)
        type_info = self.get_type_by_name(type_name)
        if not type_info:
            return 0
        parents = type_info['config'].get('is_a')
        if not parents:
            return 0
        if isinstance(parents, str):
            parents = [parents]
        if not isinstance(parents, list):
            return 0
        return 1 + max((self.get_type_depth(parent, visited) for parent in parents), default=0)

    def get_schema_type_info(self, input_label: str) -> Optional[Dict]:
        matches = self.get_schema_type_infos(input_label)
        return matches[0] if matches else None

    def get_schema_type_infos(self, input_label: str) -> List[Dict]:
        matches = []
        for type_name, type_config in self.schema_config.items():
            if not isinstance(type_config, dict):
                continue
            labels = type_config.get('input_label')
            if labels == input_label or (isinstance(labels, list) and input_label in labels):
                matches.append({
                    'name': type_name,
                    'config': type_config
                })
        return matches

    def get_type_by_name(self, type_name: str) -> Optional[Dict]:
        if type_name in self.schema_config:
            config = self.schema_config[type_name]
            if isinstance(config, dict):
                return {'name': type_name, 'config': config}
        return None

    def get_all_properties_for_type(self, type_name: str, visited: Set[str] = None) -> Dict[str, Any]:
        if visited is None:
            visited = set()
        if type_name in visited:
            return {}
        visited.add(type_name)

        type_info = self.get_type_by_name(type_name)
        if not type_info:
            return {}

        config = type_info['config']
        all_props = {}

        if config.get('inherit_properties', False):
            parent_type = config.get('is_a')
            if parent_type:
                if isinstance(parent_type, str):
                    parent_types = [parent_type]
                elif isinstance(parent_type, list):
                    parent_types = parent_type
                else:
                    parent_types = []
                for parent in parent_types:
                    parent_props = self.get_all_properties_for_type(parent, visited)
                    all_props.update(parent_props)

        direct_props = config.get('properties', {})
        all_props.update(direct_props)
        return all_props

    def get_valid_properties(
        self,
        input_label: str,
        adapter_props: Set[str],
        type_info: Optional[Dict] = None,
    ) -> Dict[str, str]:
        type_info = type_info or self.get_schema_type_info(input_label)
        if not type_info:
            return {}

        type_name = type_info['name']
        schema_props = self.get_all_properties_for_type(type_name)
        valid_props = {}

        for prop in sorted(adapter_props):
            if prop in schema_props:
                prop_config = schema_props[prop]
                if isinstance(prop_config, dict):
                    prop_type = prop_config.get('type', 'str')
                else:
                    prop_type = 'str'
                valid_props[prop] = self.normalize_property_type(prop_type)

        return valid_props

    @staticmethod
    def normalize_property_type(prop_type: Any) -> str:
        """Normalize BioCypher property types to the shorthand used in datasource schemas."""
        if not isinstance(prop_type, str):
            return 'str'
        return {
            'string': 'str',
            'integer': 'int',
            'double': 'float',
            'boolean': 'bool',
        }.get(prop_type, prop_type)

    def extract_base_url(self, url: str) -> str:
        if not url:
            return ""
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme and parsed.netloc:
            return f"{parsed.scheme}://{parsed.netloc}"
        return url

    def analyze_adapter_file(self, adapter_file: str) -> Optional[Dict]:
        if adapter_file in self.adapter_cache:
            return self.adapter_cache[adapter_file]

        adapter_path = self.adapters_dir / adapter_file
        if not adapter_path.exists():
            print(f"Warning: Adapter file not found: {adapter_path}")
            return None

        analyzer = AdapterAnalyzer(adapter_path)
        metadata = analyzer.get_metadata_from_init()

        source = metadata.get('source', adapter_file.replace('_adapter.py', '').replace('_', ' ').title())
        source_url = metadata.get('source_url', '')

        # Override source for ontology adapters to group them under OBO Foundry
        if analyzer.inherits_from_ontology_adapter():
            source = 'OBO Foundry'
            if not source_url or 'obolibrary.org' in source_url or 'ebi.ac.uk/efo' in source_url:
                source_url = 'http://www.obofoundry.org/'

        result = {
            'analyzer': analyzer,
            'metadata': metadata,
            'source_name': source,
            'source_url': source_url,
        }

        self.adapter_cache[adapter_file] = result
        return result

    def discover_species_adapters(
        self,
        species: str,
        existing_adapter_config: Dict,
    ) -> Dict:
        for name, config in existing_adapter_config.items():
            if name in self.ADAPTER_CONFIG_METADATA_KEYS:
                continue
            if not isinstance(config, dict):
                print(f"Warning: skipping non-dict adapter config entry '{name}' ({type(config).__name__})")
        existing_modules = {
            (config.get('adapter') or {}).get('module')
            for config in existing_adapter_config.values()
            if isinstance(config, dict)
        }
        discovered = {}
        species_adapter_dir = self.adapters_dir / species
        if not species_adapter_dir.exists():
            return discovered

        for adapter_path in sorted(species_adapter_dir.rglob('*_adapter.py')):
            adapter_file = str(adapter_path.relative_to(self.adapters_dir))
            module = (
                'biocypher_metta.adapters.'
                + adapter_path.relative_to(self.adapters_dir).with_suffix('').as_posix().replace('/', '.')
            )
            if module in existing_modules:
                continue

            adapter_data = self.analyze_adapter_file(adapter_file)
            if not adapter_data:
                continue

            analyzer = adapter_data['analyzer']
            has_nodes = any(
                isinstance(node, ast.FunctionDef) and node.name == 'get_nodes'
                for node in ast.walk(analyzer.tree)
            )
            has_edges = any(
                isinstance(node, ast.FunctionDef) and node.name == 'get_edges'
                for node in ast.walk(analyzer.tree)
            )
            if not has_nodes and not has_edges:
                continue

            adapter_name = adapter_path.stem.removesuffix('_adapter')
            discovered[adapter_name] = {
                'adapter': {
                    'module': module,
                    'cls': analyzer.get_adapter_class_name(),
                    'args': {},
                },
                'nodes': has_nodes,
                'edges': has_edges,
            }

        return discovered

    def generate_all_schemas(self, filter_adapters: Optional[List[str]] = None, filter_modules: Optional[List[str]] = None, filter_sources: Optional[List[str]] = None):
        if filter_adapters:
            print(f"Analyzing specific adapters: {', '.join(filter_adapters)}")
        elif filter_modules:
            print(f"Analyzing adapters by module: {', '.join(filter_modules)}")
        elif filter_sources:
            print(f"Analyzing specific data sources: {', '.join(filter_sources)}")
        else:
            print(f"Analyzing adapters and generating schemas...")

        by_source = defaultdict(list)

        for adapter_name, config in self.adapter_config.items():
            if adapter_name in self.ADAPTER_CONFIG_METADATA_KEYS:
                continue
            if not isinstance(config, dict):
                print(f"Warning: skipping non-dict adapter config entry '{adapter_name}' ({type(config).__name__})")
                continue

            # Filter by adapter name if specified
            if filter_adapters and adapter_name not in filter_adapters:
                continue

            adapter_info = config.get('adapter', {})
            module = adapter_info.get('module', '')
            cls = adapter_info.get('cls', '')

            if not module or not cls:
                continue

            # Filter by module name if specified
            module_name = module.split('.')[-1]
            if filter_modules and module_name not in filter_modules:
                continue

            adapter_file = self.get_adapter_file_from_module(module)
            adapter_data = self.analyze_adapter_file(adapter_file)
            if not adapter_data:
                continue

            source_name = self.SOURCE_NAME_ALIASES.get(
                adapter_data['source_name'],
                adapter_data['source_name'],
            )

            # Filter by source name if specified
            if filter_sources and source_name not in filter_sources:
                continue

            by_source[source_name].append({
                'adapter_name': adapter_name,
                'adapter_file': adapter_file,
                'config': config,
                'adapter_data': adapter_data
            })

        if not by_source:
            print(f"\n⚠ No adapters matched the filter criteria")
            return

        for source_name, adapters in sorted(by_source.items()):
            print(f"\nGenerating schema for: {source_name}")
            self.generate_schema(source_name, adapters)

        print(f"\n✓ Schema generation complete! Output directory: {self.output_dir}")
        print(f"  Generated {len(by_source)} schema files")

    def generate_schema(self, source_name: str, adapters: List[Dict]):
        source_urls = [a['adapter_data']['source_url'] for a in adapters if a['adapter_data']['source_url']]
        website = self.extract_base_url(source_urls[0]) if source_urls else ''

        # Check if output file already exists
        filename_source = self.SOURCE_FILENAME_ALIASES.get(source_name, source_name)
        output_filename = filename_source.replace(' ', '_').replace('-', '_').replace('/', '_') + '.yaml'
        output_path = self.output_dir / output_filename

        # Load existing schema if it exists
        if output_path.exists():
            with open(output_path, 'r') as f:
                existing_schema = yaml.safe_load(f)
            nodes = existing_schema.get('nodes', {})
            relationships = existing_schema.get('relationships', {})
            schema = {
                'name': existing_schema.get('name', source_name),
                # Prefer the freshly computed value; only fall back to what's already on
                # disk when this run couldn't resolve a URL from any adapter (e.g. a
                # filtered run that didn't touch the adapter carrying it).
                'website': website if website else existing_schema.get('website', '')
            }
        else:
            schema = {
                'name': source_name,
            }
            if website:
                schema['website'] = website
            nodes = {}
            relationships = {}

        for adapter_info in adapters:
            adapter_name = adapter_info['adapter_name']
            adapter_cfg = adapter_info['config']
            adapter_data = adapter_info['adapter_data']
            analyzer = adapter_data['analyzer']

            writes_nodes = adapter_cfg.get('nodes', False)
            writes_edges = adapter_cfg.get('edges', False)
            is_ontology_adapter = analyzer.inherits_from_ontology_adapter()
            adapter_args = adapter_cfg.get('adapter', {}).get('args', {})

            if is_ontology_adapter:
                ontology_type = adapter_args.get('type', '')
                if ontology_type == 'node':
                    pass

            adapter_source_url = adapter_data['source_url'] if adapter_data['source_url'] else ''
            # Re-resolve self.attr = ClassName.DICT[ctor_arg]-style source_url precisely
            # for this specific block, when its args give the ctor_arg a literal value
            # the static analysis above couldn't see (it only had a naive first-value
            # fallback -- see AdapterAnalyzer.dict_lookups).
            lookup = analyzer.dict_lookups.get('source_url')
            if lookup:
                dict_name, arg_name = lookup
                key = adapter_args.get(arg_name)
                dict_value = analyzer.class_attributes.get(dict_name)
                if key is not None and isinstance(dict_value, dict) and key in dict_value:
                    adapter_source_url = dict_value[key]

            for label in self.get_labels_for_adapter_config(adapter_name, adapter_cfg, adapter_data):
                type_infos = self.get_schema_type_infos(label)
                if not type_infos:
                    print(f"  Warning: No schema config found for label: {label} (adapter: {adapter_name})")
                    continue

                for type_info in type_infos:
                    type_config = type_info['config']
                    type_name = type_info['name']
                    output_label = type_config.get('output_label')
                    is_edge = type_config.get('represented_as') == 'edge'

                    # Process nodes
                    if writes_nodes and not is_edge:
                        node_props = analyzer.get_properties_for_label('get_nodes', label, 1, 2)
                        if not node_props:
                            node_props = analyzer.get_node_properties()

                        if is_ontology_adapter:
                            ontology_adapter_path = self.adapters_dir / 'ontologies_adapter.py'
                            parent_props = analyzer.get_parent_class_properties('get_nodes', ontology_adapter_path)
                            node_props = node_props.union(parent_props)

                        valid_props = self.get_valid_properties(label, node_props, type_info=type_info)
                        description = type_config.get('description', '')

                        # Add or update node
                        if type_name not in nodes:
                            nodes[type_name] = {}
                        nodes[type_name]['url'] = adapter_source_url
                        nodes[type_name]['input_label'] = label
                        if output_label:
                            nodes[type_name]['output_label'] = output_label
                        else:
                            nodes[type_name].pop('output_label', None)
                        if description:
                            nodes[type_name]['description'] = description.strip()
                        # Merge properties: refresh types redetected this run, keep others
                        # learned from earlier (e.g. filtered) runs.
                        if valid_props:
                            if 'properties' not in nodes[type_name]:
                                nodes[type_name]['properties'] = {}
                            nodes[type_name]['properties'].update(valid_props)

                    # Process edges
                    elif writes_edges and is_edge:
                        edge_props = analyzer.get_properties_for_label('get_edges', label, 2, 3)
                        if not edge_props:
                            edge_props = analyzer.get_edge_properties()

                        if is_ontology_adapter:
                            ontology_adapter_path = self.adapters_dir / 'ontologies_adapter.py'
                            parent_props = analyzer.get_parent_class_properties('get_edges', ontology_adapter_path)
                            edge_props = edge_props.union(parent_props)

                        valid_props = self.get_valid_properties(label, edge_props, type_info=type_info)
                        description = type_config.get('description', '')
                        source = type_config.get('source')
                        target = type_config.get('target')

                        # Add or update relationship
                        if type_name not in relationships:
                            relationships[type_name] = {}
                        relationships[type_name]['url'] = adapter_source_url
                        relationships[type_name]['input_label'] = label
                        if output_label:
                            relationships[type_name]['output_label'] = output_label
                        else:
                            relationships[type_name].pop('output_label', None)
                        if description:
                            relationships[type_name]['description'] = description.strip()
                        if source:
                            relationships[type_name]['source'] = source
                        if target:
                            relationships[type_name]['target'] = target
                        # Merge properties: refresh types redetected this run, keep others
                        # learned from earlier (e.g. filtered) runs.
                        if valid_props:
                            if 'properties' not in relationships[type_name]:
                                relationships[type_name]['properties'] = {}
                            relationships[type_name]['properties'].update(valid_props)

        if nodes:
            schema['nodes'] = nodes
        if relationships:
            schema['relationships'] = relationships

        if nodes or relationships:
            with open(output_path, 'w') as f:
                yaml.dump(schema, f, Dumper=FlowStyleListDumper, default_flow_style=False, sort_keys=False, allow_unicode=True)

            print(f"  ✓ Generated: {output_filename}")
            print(f"    - Nodes: {len(nodes)}, Relationships: {len(relationships)}")
        else:
            print(f"  Warning: No nodes or relationships found for {source_name}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate data source schema YAML files from adapters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--species',
        choices=sorted(SchemaGenerator.SPECIES_CONFIGS),
        help='Generate all configured data source schemas for a species using repo defaults.',
    )
    parser.add_argument('--schema-config', help='Path to schema configuration YAML file')
    parser.add_argument('--adapter-config', help='Path to adapter configuration YAML file')
    parser.add_argument('--adapters-dir', default='biocypher_metta/adapters', help='Directory containing adapter Python files')
    parser.add_argument('--output-dir', help='Output directory for generated schema files')
    parser.add_argument('--adapter', action='append', help='Generate schema only for specific adapter(s). Can be used multiple times.')
    parser.add_argument('--module', action='append', help='Generate schema for all adapters using specific module(s). Example: candidate_cis_regulatory_promoter_adapter. Can be used multiple times.')
    parser.add_argument('--source', action='append', help='Generate schema only for specific data source(s). Can be used multiple times.')
    parser.add_argument(
        '--include-inactive-adapters',
        action=argparse.BooleanOptionalAction,
        default=None,
        help='Also include commented-out adapter blocks from the adapter config. Defaults to on for --species.',
    )

    args = parser.parse_args()

    if args.species:
        defaults = SchemaGenerator.SPECIES_CONFIGS[args.species]
        args.schema_config = args.schema_config or defaults['schema_config']
        args.adapter_config = args.adapter_config or defaults['adapter_config']
        args.output_dir = args.output_dir or defaults['output_dir']

    missing_args = [
        name
        for name in ('schema_config', 'adapter_config', 'adapters_dir', 'output_dir')
        if not getattr(args, name)
    ]
    if missing_args:
        parser.error(
            "the following arguments are required unless --species is used: "
            + ", ".join(f"--{name.replace('_', '-')}" for name in missing_args)
        )

    if not Path(args.schema_config).exists():
        print(f"Error: Schema config not found: {args.schema_config}")
        sys.exit(1)
    if not Path(args.adapter_config).exists():
        print(f"Error: Adapter config not found: {args.adapter_config}")
        sys.exit(1)
    if not Path(args.adapters_dir).exists():
        print(f"Error: Adapters directory not found: {args.adapters_dir}")
        sys.exit(1)

    schema_config_data = None
    adapter_config_data = None
    if args.species:
        schema_config_data = SchemaGenerator.load_schema_config(
            args.schema_config,
            include_primer=True,
        )
        include_inactive = True if args.include_inactive_adapters is None else args.include_inactive_adapters
        if include_inactive:
            with open(args.adapter_config) as f:
                adapter_config_data = load_yaml_with_includes(f) or {}
            inactive_adapters = SchemaGenerator.load_commented_adapter_config(args.adapter_config)
            for adapter_name, adapter_cfg in inactive_adapters.items():
                adapter_config_data.setdefault(adapter_name, adapter_cfg)
            discovery_generator = SchemaGenerator(
                args.schema_config,
                args.adapter_config,
                args.adapters_dir,
                args.output_dir,
                adapter_config_data=adapter_config_data,
                schema_config_data=schema_config_data,
            )
            discovered_adapters = discovery_generator.discover_species_adapters(
                args.species,
                adapter_config_data,
            )
            for adapter_name, adapter_cfg in discovered_adapters.items():
                adapter_config_data.setdefault(adapter_name, adapter_cfg)

    generator = SchemaGenerator(
        args.schema_config,
        args.adapter_config,
        args.adapters_dir,
        args.output_dir,
        adapter_config_data=adapter_config_data,
        schema_config_data=schema_config_data,
    )

    generator.generate_all_schemas(
        filter_adapters=args.adapter,
        filter_modules=args.module,
        filter_sources=args.source
    )


if __name__ == '__main__':
    main()
