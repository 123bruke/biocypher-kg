#!/usr/bin/env python3
"""
Verify genome assembly configurations across all species in BioCypher.

Ensures that:
1. Every species has the correct, verified genome assembly configured.
2. Non-human species (dmel, cel, mmu, rno) do NOT receive human assembly tags (GRCh38/hg38).
3. All configured assemblies exist in ALLOWED_ASSEMBLIES in helpers.py.
4. EPD adapter generates identifiers matching the verified assembly tag for each species.
5. Both production and sample configurations match.

Exits with code 0 on success, or code 1 if any assembly issue is found.
"""

import sys
import os
import argparse
from pathlib import Path
import yaml

# Expected verified assemblies per species
VERIFIED_SPECIES_ASSEMBLIES = {
    'hsa': 'GRCh38',
    'dmel': 'dm6',
    'cel': 'ce6',
    'mmu': 'mm10',
    'rno': 'rn6',
}

HUMAN_ASSEMBLIES = {'GRCh38', 'hg38', 'GRCh37', 'hg19'}


def load_yaml_config(file_path: Path) -> dict:
    with open(file_path, 'r', encoding='utf-8') as f:
        # Ignore custom yaml tags if present (e.g. !include)
        class SafeLoaderIgnoreUnknown(yaml.SafeLoader):
            pass

        def unknown_constructor(loader, tag_suffix, node):
            return loader.construct_scalar(node)

        SafeLoaderIgnoreUnknown.add_multi_constructor('', unknown_constructor)
        return yaml.load(f, Loader=SafeLoaderIgnoreUnknown) or {}


def check_assembly_configurations(repo_root: Path) -> bool:
    errors = []
    checks_passed = 0

    print("========================================================")
    print("  BioCypher Species Genome Assembly Verification")
    print(f"  Target Root: {repo_root.resolve()}")
    print("========================================================\n")

    # 1. Check ALLOWED_ASSEMBLIES in helpers.py
    helpers_file = repo_root / "biocypher_metta" / "adapters" / "helpers.py"
    allowed_assemblies = []
    if helpers_file.exists():
        # Add repo root to sys.path to import helpers
        sys.path.insert(0, str(repo_root))
        try:
            from biocypher_metta.adapters.helpers import (
                ALLOWED_ASSEMBLIES,
                build_regulatory_region_id,
            )
            allowed_assemblies = list(ALLOWED_ASSEMBLIES)
            print(f"[OK] helpers.ALLOWED_ASSEMBLIES loaded: {allowed_assemblies}")
            checks_passed += 1
        except Exception as e:
            errors.append(f"Failed to import ALLOWED_ASSEMBLIES from helpers: {e}")
    else:
        errors.append(f"helpers.py not found at {helpers_file}")

    for species, expected_assembly in VERIFIED_SPECIES_ASSEMBLIES.items():
        if expected_assembly not in allowed_assemblies:
            errors.append(
                f"Verified assembly '{expected_assembly}' for species '{species}' "
                f"is missing from ALLOWED_ASSEMBLIES ({allowed_assemblies})"
            )
        else:
            checks_passed += 1

    # 2. Check each species config and sample config
    config_dir = repo_root / "config"
    if not config_dir.exists():
        errors.append(f"Config directory not found: {config_dir}")
        for err in errors:
            print(f"[ERROR] {err}")
        return False

    for species, expected_assembly in VERIFIED_SPECIES_ASSEMBLIES.items():
        species_dir = config_dir / species
        if not species_dir.exists():
            print(f"[WARN] Species directory {species_dir} does not exist, skipping.")
            continue

        config_files = [
            species_dir / f"{species}_adapters_config.yaml",
            species_dir / f"{species}_adapters_config_sample.yaml",
        ]

        for cfg_file in config_files:
            if not cfg_file.exists():
                print(f"[WARN] Config file {cfg_file.name} does not exist, skipping.")
                continue

            rel_cfg = cfg_file.relative_to(repo_root)
            try:
                data = load_yaml_config(cfg_file)
            except Exception as e:
                errors.append(f"Failed to parse YAML {rel_cfg}: {e}")
                continue

            # Check epd entries
            for epd_key in ['epd_promoter', 'epd_promoter_regulates_gene']:
                if epd_key not in data:
                    continue

                entry = data[epd_key]
                args = (entry.get('adapter') or {}).get('args') or {}
                configured_assembly = args.get('assembly')

                # Check 1: Non-human species MUST configure an assembly
                if species != 'hsa':
                    if configured_assembly is None:
                        errors.append(
                            f"[{rel_cfg}] {epd_key}: missing 'assembly' argument. "
                            f"Non-human species '{species}' would fall back to default human GRCh38!"
                        )
                        continue

                    # Check 2: Non-human species MUST NOT use human assembly
                    if configured_assembly in HUMAN_ASSEMBLIES:
                        errors.append(
                            f"[{rel_cfg}] {epd_key}: invalid assembly '{configured_assembly}' for non-human species '{species}'. "
                            f"Non-human elements must NOT receive a human assembly tag!"
                        )
                        continue

                    # Check 3: Configured assembly must match verified assembly
                    if configured_assembly != expected_assembly:
                        errors.append(
                            f"[{rel_cfg}] {epd_key}: expected assembly '{expected_assembly}' for '{species}', "
                            f"but found '{configured_assembly}'."
                        )
                        continue

                else:
                    # For human, if configured, it must be GRCh38
                    if configured_assembly is not None and configured_assembly != 'GRCh38':
                        errors.append(
                            f"[{rel_cfg}] {epd_key}: human assembly should be GRCh38, found '{configured_assembly}'."
                        )
                        continue

                # Check 4: Configured assembly must be in ALLOWED_ASSEMBLIES
                if configured_assembly and configured_assembly not in allowed_assemblies:
                    errors.append(
                        f"[{rel_cfg}] {epd_key}: assembly '{configured_assembly}' is not in ALLOWED_ASSEMBLIES."
                    )
                    continue

                print(f"[OK] {rel_cfg} -> {epd_key}: assembly='{configured_assembly}'")
                checks_passed += 1

    # 3. Verify ID generation logic for all species
    if 'build_regulatory_region_id' in locals():
        print("\nVerifying ID Generation for each species:")
        for species, expected_assembly in VERIFIED_SPECIES_ASSEMBLIES.items():
            test_chr = "chr1"
            test_start = 1000000
            test_end = 1000060

            if species == 'hsa':
                # Default (no arg)
                id_default = build_regulatory_region_id(test_chr, test_start, test_end)
                id_explicit = build_regulatory_region_id(test_chr, test_start, test_end, assembly='GRCh38')
                if id_default != f"{test_chr}_{test_start}_{test_end}_GRCh38":
                    errors.append(f"Human default ID builder returned unexpected value: {id_default}")
                elif id_explicit != id_default:
                    errors.append(f"Human explicit ID {id_explicit} does not match default ID {id_default}")
                else:
                    print(f"[OK] Species '{species}' ID generation: {id_default}")
                    checks_passed += 1
            else:
                id_species = build_regulatory_region_id(test_chr, test_start, test_end, assembly=expected_assembly)
                if not id_species.endswith(f"_{expected_assembly}"):
                    errors.append(f"Species '{species}' ID {id_species} does not end with '_{expected_assembly}'")
                elif any(h in id_species for h in HUMAN_ASSEMBLIES):
                    errors.append(f"Species '{species}' ID {id_species} incorrectly contains a human assembly tag!")
                else:
                    print(f"[OK] Species '{species}' ID generation: {id_species}")
                    checks_passed += 1

        # Test unsupported assembly rejection
        try:
            build_regulatory_region_id(test_chr, test_start, test_end, assembly='INVALID_ASM')
            errors.append("build_regulatory_region_id failed to reject unsupported assembly 'INVALID_ASM'")
        except ValueError:
            print("[OK] Unsupported assembly correctly raises ValueError")
            checks_passed += 1

    print("\n--------------------------------------------------------")
    if errors:
        print(f"FAILED: {len(errors)} error(s) detected:")
        for err in errors:
            print(f"  [ERROR] {err}")
        print("--------------------------------------------------------\n")
        return False
    else:
        print(f"SUCCESS: All {checks_passed} assembly checks passed successfully.")
        print("--------------------------------------------------------\n")
        return True


def main():
    parser = argparse.ArgumentParser(
        description="Verify genome assembly configuration across BioCypher species."
    )
    parser.add_argument(
        'root_dir',
        nargs='?',
        default='.',
        help="Root directory of the repository (default: .)"
    )
    args = parser.parse_args()

    repo_root = Path(args.root_dir).resolve()
    if not (repo_root / "config").exists() and (repo_root / "biocypher-kg-main" / "config").exists():
        repo_root = repo_root / "biocypher-kg-main"

    success = check_assembly_configurations(repo_root)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()

