# nf-clone-variation

A Nextflow pipeline for visualizing and quantifying within-clone HIV sequence variation relative to a derived clone consensus, including longitudinal changes in allele frequency.

It uses outputs from [nf-clone-consensus](https://github.com/ekankaka/nf-clone-consensus). No realignment is performed.

## Requirements

- Nextflow and Java
- Docker, Singularity or Apptainer (or Python 3.10+ and Matplotlib for native execution)
- SLURM for HPC execution

## Inputs

| Parameter | Description |
|---|---|
| --input | clone_summary.csv |
| --clone_alignment_folder | Directory containing clone alignments |
| --read_assignments | review/read_assignments.csv |
| --min_unique_reads | Minimum distinct PCR families per clone (default: 3) |
| --min_family_fraction | Minimum within-family agreement (default: 0.60; match the upstream setting if changed) |
| --outdir | Output directory (default: results) |

An alignment is named <consensus_header>_alignment.fasta or <consensus_header>.fasta. The first sequence must be the aligned consensus. Read assignments provide sample IDs, PCR IDs and timepoints. The clone summary must contain the upstream consensus_sequence, n_reads, n_unique_reads, patient_id, clone_id and consensus_header fields.

## Run

Local synthetic example (requires Docker):

```bash
git clone https://github.com/ekankaka/nf-clone-variation.git
cd nf-clone-variation

nextflow run main.nf \
    --input tests/data/clone_summary.csv \
    --clone_alignment_folder tests/data/clone_alignments \
    --read_assignments tests/data/read_assignments.csv \
    --outdir results_test \
    -profile local,docker
```

HPC example:

```bash
INFILE=/path/to/nf-clone-consensus/results/clone_summary.csv
CLONE_ALIGNMENT_FOLDER=/path/to/nf-clone-consensus/results/clone_alignments
READ_ASSIGNMENTS=/path/to/nf-clone-consensus/results/review/read_assignments.csv

nextflow run main.nf \
    --input "$INFILE" \
    --clone_alignment_folder "$CLONE_ALIGNMENT_FOLDER" \
    --read_assignments "$READ_ASSIGNMENTS" \
    --outdir results \
    -profile slurm,singularity \
    -c /path/to/hpc.config \
    -resume
```

Use local or slurm with docker, singularity or apptainer. HPC queue, account, work directory and cache settings belong in an external HPC configuration or command line. Use an absolute output path on HPC. The container tag is a Python 3.12 Jupyter SciPy image providing Matplotlib; for a frozen production environment, override the container with a verified digest.

## Method

1. Select clones with at least three independent PCR families.
2. Verify the alignment and map each read to its sample, timepoint and PCR family.
3. Treat terminal gaps and ambiguous read bases as missing information. Use A/C/G/T and internal gaps as informative calls. Mask ambiguous consensus bases, including N.
4. Collapse duplicate reads to one call per (sample_id, PCR_ID) family per column using the minimum agreement fraction. Missing calls do not vote.
5. Quantify substitutions and deleted bases at unambiguous consensus positions with at least two informative families. Count insertion alleles at internal consensus gap blocks with unambiguous flanks and at least two families spanning the junction.
6. Report fully resolved deletion events separately from per-base deletion burden. Exact event frequencies require at least two families with informative calls across the event and base calls at both flanks. Count contiguous deleted bases as one event; alternative inserted sequences are separate alleles.
7. Compute frequencies using locus-specific family denominators, both pooled and per timepoint. No coverage gives NA, not zero. Keep singleton variants and report independent support.
8. Plot all input reads ordered by timepoint (including duplicates), with raw substitutions, internal gaps and candidate insertions, followed by family-level variant frequencies and coverage.

Plots use 1-based alignment columns. Tables also include ungapped consensus positions. The family frequency panel shows per-column difference frequencies, not exact multi-base indel allele frequencies. Raw read differences may not survive family-level filters. Numeric timepoints sort numerically; other labels sort lexically.

## Outputs

```text
results/
├── clone_variation_summary.csv
├── variants.csv
├── variants_by_timepoint.csv
├── plots/
│   ├── <consensus_header>.png
│   └── <consensus_header>.pdf
└── review/
    └── skipped_clones.csv
```

The clone summary preserves upstream columns and appends coverage-adjusted variation statistics. Clones with zero differences remain in the summary. Eligible clones with missing or inconsistent alignments are listed in skipped_clones.csv; clones below the threshold are not analyzed. See [data dictionary](docs/data_dictionary.md).

## Tests

The synthetic input includes a six-read clone with five independent families, two timepoints, a PCR duplicate, a terminal fragment, a masked base, one substitution, one two-base insertion and one two-base deletion supported by two families. A second clone is invariant; a third eligible clone is missing its alignment; a fourth is below the threshold.

Run the analysis and merge tests with Python and Matplotlib:

```bash
python tests/test_variation.py
```

After running the synthetic Nextflow example, verify its outputs:

```bash
python tests/test_variation.py results_test
```

## Interpretation and limitations

The pooled consensus is a comparison sequence, not an ancestral sequence. Variants are differences relative to that consensus, not necessarily new mutations. Ambiguous consensus positions are deliberately excluded, even if they contain true polymorphism. Internal gaps may reflect alignment artifacts; indels are not left-normalized across repetitive sequences. Exact deletion events that cannot be bounded are excluded from the event table but their evaluable deleted bases contribute to base burden. FASTA lacks quality scores; technical artifacts cannot be ruled out. PCR-family support is not proof of biological variation. Longitudinal comparisons require adequate coverage at each timepoint and confirmation of ART suppression using clinical data. No inference of viral replication or mutation rate is made by this pipeline.

License: MIT.
