#!/usr/bin/env python3
"""
Scan python files to find call sites to genome assembly-dependent ID builders:
- build_regulatory_region_id
- build_variant_id
- build_variant_id_from_hgvs

Reports which callers explicitly pass an assembly parameter and which rely on the default.
"""

import ast
import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Any


BUILDER_PARAM_COUNTS = {
    'build_regulatory_region_id': {
        'assembly_pos_index': 3,  # (chr, pos_start, pos_end, assembly)
    },
    'build_variant_id': {
        'assembly_pos_index': 4,  # (chr, pos_first_ref_base, ref_seq, alt_seq, assembly)
    },
    'build_variant_id_from_hgvs': {
        'assembly_pos_index': 2,  # (hgvs_id, validate, assembly)
    },
}


class AssemblyCallVisitor(ast.NodeVisitor):
    def __init__(self, filepath: str, rel_path: str):
        self.filepath = filepath
        self.rel_path = rel_path
        self.scope_stack: List[str] = []
        self.call_sites: List[Dict[str, Any]] = []

    def visit_ClassDef(self, node: ast.ClassDef):
        self.scope_stack.append(node.name)
        self.generic_visit(node)
        self.scope_stack.pop()

    def visit_FunctionDef(self, node: ast.FunctionDef):
        self.scope_stack.append(node.name)
        self.generic_visit(node)
        self.scope_stack.pop()

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef):
        self.scope_stack.append(node.name)
        self.generic_visit(node)
        self.scope_stack.pop()

    def visit_Call(self, node: ast.Call):
        func_name = None
        if isinstance(node.func, ast.Name):
            func_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            func_name = node.func.attr

        if func_name in BUILDER_PARAM_COUNTS:
            spec = BUILDER_PARAM_COUNTS[func_name]
            passes_assembly = False
            assembly_expr = None

            # Check keyword arguments
            for kw in node.keywords:
                if kw.arg == 'assembly':
                    passes_assembly = True
                    try:
                        assembly_expr = ast.unparse(kw.value)
                    except Exception:
                        assembly_expr = "<kwarg>"
                    break

            # Check positional arguments
            if not passes_assembly:
                if len(node.args) > spec['assembly_pos_index']:
                    passes_assembly = True
                    try:
                        assembly_expr = ast.unparse(node.args[spec['assembly_pos_index']])
                    except Exception:
                        assembly_expr = "<pos_arg>"

            scope = '.'.join(self.scope_stack) if self.scope_stack else '<module>'
            self.call_sites.append({
                'file': self.rel_path,
                'line': node.lineno,
                'func_called': func_name,
                'scope': scope,
                'passes_assembly': passes_assembly,
                'assembly_expr': assembly_expr if passes_assembly else 'None (default GRCh38)',
            })

        self.generic_visit(node)


def scan_directory(target_dir: Path) -> List[Dict[str, Any]]:
    all_calls = []
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.endswith('.py'):
                full_path = Path(root) / file
                try:
                    rel_path = full_path.relative_to(target_dir)
                except ValueError:
                    rel_path = full_path

                try:
                    with open(full_path, 'r', encoding='utf-8') as f:
                        tree = ast.parse(f.read(), filename=str(full_path))
                    visitor = AssemblyCallVisitor(str(full_path), str(rel_path))
                    visitor.visit(tree)
                    all_calls.extend(visitor.call_sites)
                except Exception as e:
                    print(f"Warning: could not parse {full_path}: {e}", file=sys.stderr)

    return all_calls


def main():
    parser = argparse.ArgumentParser(
        description="Find all callers of genome assembly-dependent ID builders in BioCypher adapters."
    )
    parser.add_argument(
        'target_dir',
        nargs='?',
        default='biocypher_metta',
        help="Target directory to scan (default: biocypher_metta)"
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help="List each caller in detail"
    )
    args = parser.parse_args()

    target_path = Path(args.target_dir)
    if not target_path.exists():
        # Try finding it relative to current working dir or repo root
        candidates = [
            Path('.') / args.target_dir,
            Path(__file__).resolve().parent.parent.parent / args.target_dir,
        ]
        found = False
        for cand in candidates:
            if cand.exists():
                target_path = cand
                found = True
                break
        if not found:
            print(f"Error: Target directory '{args.target_dir}' not found.", file=sys.stderr)
            sys.exit(1)

    calls = scan_directory(target_path)

    # Sort calls by file, line
    calls.sort(key=lambda x: (x['file'], x['line']))

    explicit_calls = [c for c in calls if c['passes_assembly']]
    default_calls = [c for c in calls if not c['passes_assembly']]

    print(f"\n========================================================")
    print(f"  Assembly Builder Call Sites Scan: {target_path}")
    print(f"========================================================\n")
    print(f"Total call sites found: {len(calls)}")
    print(f"  - Passing assembly explicitly: {len(explicit_calls)}")
    print(f"  - Relying on default (GRCh38): {len(default_calls)}\n")

    if args.list:
        header = f"{'File':<40} {'Line':<6} {'Function Called':<28} {'Scope':<32} {'Assembly Arg'}"
        print(header)
        print("-" * len(header))
        for c in calls:
            print(
                f"{c['file']:<40} "
                f"{c['line']:<6} "
                f"{c['func_called']:<28} "
                f"{c['scope']:<32} "
                f"{c['assembly_expr']}"
            )
        print("-" * len(header) + "\n")

    sys.exit(0)


if __name__ == '__main__':
    main()

