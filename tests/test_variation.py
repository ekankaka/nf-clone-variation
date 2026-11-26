#!/usr/bin/env python3
"""Run from the repository root: python tests/test_variation.py."""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data"
SCRIPT = ROOT / "bin/assess_clone_variation.py"


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def check_results(folder):
    summary = read(folder / "clone_variation_summary.csv")
    assert len(summary) == 2, summary
    varied, stable = summary
    expected = {"n_callable_consensus_bases": 11, "n_evaluable_positions": 11,
                "n_callable_family_bases": 54, "n_substitution_sites": 1,
                "n_substitutions": 2, "n_insertion_sites": 1, "n_insertion_events": 2,
                "n_inserted_bases": 4, "n_deletion_sites": 1, "n_deletion_base_sites": 2,
                "n_deletion_events": 2, "n_deleted_bases": 4, "n_families_with_variants": 2,
                "n_multifamily_variants": 3, "n_multisample_variants": 0,
                "n_multitimepoint_variants": 0, "n_evaluable_insertion_junctions": 9,
                "n_callable_family_insertion_junctions": 38}
    for field, value in expected.items():
        assert int(varied[field]) == value, (field, varied[field], value)
    assert abs(float(varied["variant_bases_per_kb"]) - 10000 / 54) < 1e-6
    assert int(stable["n_callable_family_bases"]) == 12
    assert float(stable["variant_bases_per_kb"]) == 0
    variants = read(folder / "variants.csv")
    assert len(variants) == 3
    by_type = {v["variant_type"]: v for v in variants}
    for v in variants:
        assert int(v["n_supporting_families"]) == 2
        assert int(v["n_covering_families"]) == 5
        assert float(v["variant_frequency"]) == .4
    assert by_type["substitution"]["consensus_start"] == "2"
    assert by_type["insertion"]["alternative_allele"] == "GG"
    assert by_type["deletion"]["consensus_allele"] == "TA"
    timepoints = read(folder / "variants_by_timepoint.csv")
    assert len(timepoints) == 6
    for v in timepoints:
        expected = (0, 2) if v["timepoint"] == "1" else (2, 3)
        assert (int(v["n_supporting_families"]), int(v["n_covering_families"])) == expected
    skip = folder / "review/skipped_clones.csv"
    if not skip.exists():
        skip = folder / "skipped_clones.csv"
    skipped = read(skip)
    assert len(skipped) == 1
    assert skipped[0]["skip_reason"] == "Alignment file not found"
    print("Synthetic counts, rates, longitudinal frequencies and skipped clones: PASS")


def main():
    if len(sys.argv) > 1:
        check_results(Path(sys.argv[1]))
        return
    clones = read(DATA / "clone_summary.csv")
    assignments = read(DATA / "read_assignments.csv")
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        for i, row in enumerate(clones):
            if int(row["n_unique_reads"]) < 3:
                continue
            task = folder / f"clone{i}"
            task.mkdir()
            header = row["consensus_header"]
            alignment = DATA / "clone_alignments" / (header + "_alignment.fasta")
            if not alignment.exists():
                alignment = DATA / "clone_alignments" / (header + ".fasta")
            metadata = {"clone": row, "missing_alignment": not alignment.exists(),
                        "assignments": [a for a in assignments if
                                        a["patient_id"] == row["patient_id"] and a["clone_id"] == row["clone_id"]]}
            (task / "metadata.json").write_text(json.dumps(metadata))
            subprocess.run([sys.executable, str(SCRIPT), "analyze",
                            "--metadata", "metadata.json", "--alignment", str(alignment)],
                           cwd=task, check=True)
            if alignment.exists():
                for ext in ("png", "pdf"):
                    assert (task / (header + "." + ext)).stat().st_size > 0
        subprocess.run([sys.executable, str(SCRIPT), "merge", "--input", str(DATA / "clone_summary.csv")],
                       cwd=folder, check=True)
        check_results(folder)


if __name__ == "__main__":
    main()
