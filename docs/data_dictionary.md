# Data dictionary

All tables are CSV. Counts and rates refer to independent (sample_id, PCR_ID) families, not raw reads. Frequencies are fractions from 0 to 1; percent fields are from 0 to 100. NA means no informative denominator. Coordinates are 1-based.

## clone_variation_summary.csv

All original clone_summary.csv columns are preserved. The following are appended.

| Field | Definition |
|---|---|
| n_callable_consensus_bases | Consensus A/C/G/T positions, excluding gaps and ambiguity |
| n_evaluable_positions | Callable consensus positions with at least two informative families |
| percent_evaluable_consensus | 100 × evaluable positions / callable consensus bases |
| mean_family_depth | Mean informative depth over all callable consensus positions |
| n_callable_family_bases | Sum of informative family calls (including internal deletions) across evaluable positions |
| n_evaluable_insertion_junctions | Adjacent unambiguous consensus positions with at least two families spanning the junction, including ordinary non-gap junctions |
| n_callable_family_insertion_junctions | Sum of family coverage over evaluable insertion junctions |
| n_evaluable_families | Families with at least one evaluable base or insertion block |
| n_substitution_sites | Distinct evaluable consensus bases with a substitution |
| n_insertion_sites | Distinct insertion blocks with a reported insertion allele |
| n_deletion_sites | Distinct start/end intervals with a reported, bounded deletion allele |
| n_deletion_base_sites | Distinct evaluable consensus bases with an internal gap in any family, including unresolved events |
| n_substitutions | Number of substitution observations summed across families |
| n_insertion_events | Number of insertion allele observations summed across families |
| n_deletion_events | Number of fully resolved, reported deletion allele observations summed across families |
| n_inserted_bases | Number of inserted nucleotides summed across families |
| n_deleted_bases | Number of evaluable deleted nucleotides summed across families, including unresolved deletions |
| n_families_with_variants | Families with an evaluable substitution, inserted base or deleted base |
| n_multifamily_variants | Distinct reported alleles supported by at least two families |
| n_multisample_variants | Distinct reported alleles supported by at least two samples |
| n_multitimepoint_variants | Distinct reported alleles supported at two or more timepoints |
| max_variant_frequency | Largest pooled frequency among reported alleles; 0 if none reported and there are evaluable families |
| substitution_sites_per_kb | 1000 × substitution sites / evaluable positions |
| substitutions_per_kb | 1000 × substitution observations / callable family bases |
| inserted_bases_per_kb | 1000 × inserted bases / callable family bases |
| deleted_bases_per_kb | 1000 × deleted bases / callable family bases |
| variant_bases_per_kb | 1000 × (substitutions + inserted bases + deleted bases) / callable family bases |
| percent_families_with_variants | 100 × families with variants / evaluable families |

The base-normalized insertion burden measures inserted bases relative to interrogated consensus length, not insertion probability per junction. The separate junction denominator allows junction-normalized analyses. Multiple family observations at the same site contribute separately to base burden. All event counts and allele-based summaries exclude unresolved deletion boundaries; consequently deleted-base counts need not equal the sum of lengths of reported deletion events.

## variants.csv

One row per distinct reported alternative allele per clone.

| Field | Definition |
|---|---|
| patient_id, clone_id, consensus_header | Upstream clone identifiers |
| variant_type | substitution, insertion or deletion |
| alignment_start, alignment_end | Inclusive alignment interval; an insertion uses the full consensus-gap block |
| consensus_start, consensus_end | Inclusive ungapped consensus coordinates; an insertion is anchored between these two adjacent bases |
| consensus_allele | A/C/G/T for substitutions, deleted string for deletions, "-" for insertions |
| alternative_allele | Alternative nucleotide or inserted string; "-" for deletions |
| variant_length | Number of affected nucleotides |
| n_supporting_families | Distinct families with the exact allele |
| n_covering_families | Informative families at the locus (at least two pooled) |
| variant_frequency | Supporting / covering families |
| n_supporting_samples, n_supporting_timepoints | Number of distinct samples or timepoints supporting the exact allele |
| supporting_samples, supporting_timepoints | Semicolon-separated sample/timepoint labels |
| supporting_families | JSON list of [sample_id, PCR_ID] pairs |

For substitutions, informative families may have a different nucleotide or internal deletion at that position. For insertions, all columns in the gap block must be informative and both flanking family calls must be bases. For a bounded deletion, all columns between its flanks must be informative and both flanks must have family base calls. Overlapping longer deletions that remove a flank do not contribute to that exact-event denominator.

## variants_by_timepoint.csv

Repeats the variant identifiers, coordinates and alleles from variants.csv for every timepoint represented in the clone, adding timepoint, n_supporting_families, n_covering_families and variant_frequency. Zero support with positive coverage is 0; absent informative coverage is NA. Only pooled reported alleles are included. No temporal significance test is performed.

## review/skipped_clones.csv

Original clone fields plus skip_reason. Only eligible clones that could not be evaluated are listed. Missing required CSV fields, unsafe filenames and unexpected runtime failures stop the workflow rather than silently generating incomplete output.

## Plots

The read panel displays raw alignment differences, with gray indicating informative coverage, lighter shading indicating masked consensus columns and blank stretches indicating missing read coverage. Reads are sorted by timepoint, sample and PCR family; all raw reads are shown, including PCR duplicates. Labels are omitted for large clones to keep plots legible. Family-level frequency and depth panels use the same alignment coordinates; the depth panel shows informative coverage at callable consensus positions only. Candidate insertion marks in raw reads do not guarantee an independently supported, bounded insertion allele.
