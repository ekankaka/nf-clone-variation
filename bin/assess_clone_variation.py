#!/usr/bin/env python3
"""Assess aligned reads relative to the pooled clone consensus."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from itertools import groupby
from pathlib import Path

BASES = "ACGT"
STATES = BASES + "-"
IDS = ["patient_id", "clone_id", "consensus_header"]
COUNTS = """n_callable_consensus_bases n_evaluable_positions percent_evaluable_consensus
mean_family_depth n_callable_family_bases n_evaluable_insertion_junctions
n_callable_family_insertion_junctions n_evaluable_families n_substitution_sites
n_insertion_sites n_deletion_sites n_deletion_base_sites n_substitutions n_insertion_events n_deletion_events
n_inserted_bases n_deleted_bases n_families_with_variants n_multifamily_variants
n_multisample_variants n_multitimepoint_variants max_variant_frequency
substitution_sites_per_kb substitutions_per_kb inserted_bases_per_kb
deleted_bases_per_kb variant_bases_per_kb percent_families_with_variants""".split()
DETAIL = IDS + """variant_type alignment_start alignment_end consensus_start
consensus_end consensus_allele alternative_allele variant_length
n_supporting_families n_covering_families variant_frequency
n_supporting_samples n_supporting_timepoints supporting_samples supporting_timepoints
supporting_families""".split()
TIME = IDS + """variant_type alignment_start alignment_end consensus_start consensus_end
consensus_allele alternative_allele variant_length timepoint n_supporting_families
n_covering_families variant_frequency""".split()


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, fields, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif line:
            if not records:
                raise ValueError("Invalid FASTA")
            records[-1][1] += line.upper().replace(".", "-")
    return records


def read_states(seq):
    # Terminal gaps and ambiguous bases are missing, not deletions.
    seq = seq.strip()  # preserve alignment width
    left = len(seq) - len(seq.lstrip("-"))
    right = len(seq.rstrip("-"))
    return [b if left <= i < right and b in STATES else "N"
            for i, b in enumerate(seq)]


def family_call(votes, fraction):
    counts = Counter(b for b in votes if b in STATES)
    if not counts:
        return "N"
    base, count = counts.most_common(1)[0]
    return base if count / sum(counts.values()) >= fraction else "N"


def runs(indices):
    """Return consecutive runs of integer indices, inclusive."""
    for _, group in groupby(enumerate(indices), lambda x: x[1] - x[0]):
        values = [v for _, v in group]
        yield values[0], values[-1]


def ratio(n, d, scale=1):
    return round(scale * n / d, 8) if d else "NA"


def time_order(value):
    try:
        return (0, float(value), value)
    except ValueError:
        return (1, 0, value)


def analyze(row, assignments, alignment, fraction):
    records = fasta(alignment)
    if len(records) < 2 or records[0][0] != row["consensus_header"]:
        raise ValueError("Consensus header is missing or is not first")
    ref = records[0][1]
    if not ref or any(len(seq) != len(ref) for _, seq in records):
        raise ValueError("Empty or unequal alignment widths")
    if ref.replace("-", "") != row["consensus_sequence"].upper():
        raise ValueError("Consensus sequence differs from clone summary")
    metadata = {r["fasta_header"]: r for r in assignments}
    if len(metadata) != len(assignments):
        raise ValueError("Duplicate read-assignment headers")
    families, info, reads = defaultdict(list), {}, []
    seen = set()
    for header, seq in records[1:]:
        # The upstream pipeline appends this label after reverse complementation.
        name = header.removesuffix(" [reverse_complemented]")
        if name not in metadata or name in seen:
            raise ValueError("Unmatched or duplicate alignment read: " + name)
        seen.add(name)
        m = metadata[name]
        key = (m["sample_id"], m["PCR_ID"])
        if not all(key) or not m["timepoint"]:
            raise ValueError("Missing PCR family or timepoint")
        if key in info and info[key] != m["timepoint"]:
            raise ValueError("Conflicting family timepoints")
        info[key] = m["timepoint"]
        states = read_states(seq)
        families[key].append(states)
        reads.append((key, states, name))
    if (len(seen) != int(row["n_reads"]) or seen != set(metadata)
            or len(families) != int(row["n_unique_reads"])):
        raise ValueError("Alignment, assignments and clone counts disagree")

    calls = {key: [family_call(v, fraction) for v in zip(*seqs)]
             for key, seqs in families.items()}
    width = len(ref)
    depth = [sum(s[i] in STATES for s in calls.values()) for i in range(width)]
    ungapped = [i for i, b in enumerate(ref) if b != "-"]
    pos = {i: j + 1 for j, i in enumerate(ungapped)}
    callable_cols = [i for i in ungapped if ref[i] in BASES]
    eligible = {i for i in callable_cols if depth[i] >= 2}

    # Only internal insertion blocks with unambiguous consensus flanks.
    blocks = [(a, b) for a, b in runs([i for i, c in enumerate(ref) if c == "-"])
              if a > 0 and b + 1 < width and ref[a - 1] in BASES and ref[b + 1] in BASES]

    def insertion_call(s, a, b):
        if s[a - 1] not in BASES or s[b + 1] not in BASES:
            return None
        if any(c not in STATES for c in s[a:b + 1]):
            return None
        return "".join(c for c in s[a:b + 1] if c in BASES)

    # Include ordinary adjacent consensus-base junctions in the insertion
    # opportunity denominator, not just junctions with an observed alignment gap.
    junctions = [(a, b) for a, b in zip(ungapped, ungapped[1:])
                 if ref[a] in BASES and ref[b] in BASES]
    junction_depth = {(a, b): sum(insertion_call(s, a + 1, b - 1) is not None
                                  for s in calls.values()) for a, b in junctions}
    eligible_blocks = [(a, b) for a, b in blocks if junction_depth[a - 1, b + 1] >= 2]
    events = defaultdict(set)
    sub_sites, deleted_sites, affected = set(), set(), set()
    n_sub = n_del = n_ins = 0
    for key, s in calls.items():
        for i in sorted(eligible):
            if s[i] in BASES and s[i] != ref[i]:
                events[("substitution", i, i, ref[i], s[i])].add(key)
                n_sub += 1
                sub_sites.add(i)
                affected.add(key)
            elif s[i] == "-":
                n_del += 1
                deleted_sites.add(i)
                affected.add(key)
        for a, b in eligible_blocks:
            alt = insertion_call(s, a, b)
            if alt:
                events[("insertion", a, b, "-", alt)].add(key)
                n_ins += len(alt)
                affected.add(key)
        # Count complete, bounded deletion alleles on ungapped coordinates.
        # Per-base deletion burden above also includes partially resolved events.
        gaps = [j for j, i in enumerate(ungapped) if s[i] == "-" and ref[i] in BASES]
        for first, last in runs(gaps):
            if first == 0 or last == len(ungapped) - 1:
                continue
            a, b = ungapped[first], ungapped[last]
            left, right = ungapped[first - 1], ungapped[last + 1]
            if (ref[left] not in BASES or ref[right] not in BASES
                    or s[left] not in BASES or s[right] not in BASES
                    or any(i not in eligible for i in ungapped[first:last + 1])
                    or any(s[i] == "N" for i in range(left, right + 1))):
                continue
            allele = "".join(ref[i] for i in ungapped[first:last + 1])
            events[("deletion", a, b, allele, "-")].add(key)

    details, longitudinal = [], []
    ids = {k: row[k] for k in IDS}
    for (kind, a, b, original, alt), support in sorted(events.items()):
        if kind == "insertion":
            covering = {k for k, s in calls.items() if insertion_call(s, a, b) is not None}
            start, end = pos[a - 1], pos[b + 1]
        elif kind == "deletion":
            left, right = ungapped[pos[a] - 2], ungapped[pos[b]]
            covering = {k for k, s in calls.items()
                        if s[left] in BASES and s[right] in BASES
                        and all(v in STATES for v in s[left:right + 1])}
            start, end = pos[a], pos[b]
        else:
            covering = {k for k, s in calls.items() if s[a] in STATES}
            start = end = pos[a]
        if len(covering) < 2:
            continue
        common = dict(ids, variant_type=kind, alignment_start=a + 1, alignment_end=b + 1,
                      consensus_start=start, consensus_end=end, consensus_allele=original,
                      alternative_allele=alt, variant_length=len(original if kind == "deletion" else alt))
        samples = sorted({k[0] for k in support})
        times = sorted({info[k] for k in support}, key=time_order)
        details.append(dict(common, n_supporting_families=len(support),
                            n_covering_families=len(covering),
                            variant_frequency=ratio(len(support), len(covering)),
                            n_supporting_samples=len(samples), n_supporting_timepoints=len(times),
                            supporting_samples=";".join(samples), supporting_timepoints=";".join(times),
                            supporting_families=json.dumps(sorted(support))))
        for t in sorted(set(info.values()), key=time_order):
            n = sum(info[k] == t for k in support)
            d = sum(info[k] == t for k in covering)
            longitudinal.append(dict(common, timepoint=t, n_supporting_families=n,
                                     n_covering_families=d, variant_frequency=ratio(n, d)))

    ins = [v for v in details if v["variant_type"] == "insertion"]
    dels = [v for v in details if v["variant_type"] == "deletion"]
    family_bases = sum(depth[i] for i in eligible)
    evaluable_families = {k for k, s in calls.items()
                          if any(s[i] in STATES for i in eligible)
                          or any(insertion_call(s, a, b) is not None for a, b in eligible_blocks)}
    summary = dict(row, n_callable_consensus_bases=len(callable_cols),
                   n_evaluable_positions=len(eligible),
                   percent_evaluable_consensus=ratio(len(eligible), len(callable_cols), 100),
                   mean_family_depth=ratio(sum(depth[i] for i in callable_cols), len(callable_cols)),
                   n_callable_family_bases=family_bases,
                   n_evaluable_insertion_junctions=sum(d >= 2 for d in junction_depth.values()),
                   n_callable_family_insertion_junctions=sum(d for d in junction_depth.values() if d >= 2),
                   n_evaluable_families=len(evaluable_families),
                   n_substitution_sites=len(sub_sites),
                   n_insertion_sites=len({v["alignment_start"] for v in ins}),
                   n_deletion_sites=len({(v["alignment_start"], v["alignment_end"]) for v in dels}),
                   n_deletion_base_sites=len(deleted_sites),
                   n_substitutions=n_sub, n_insertion_events=sum(v["n_supporting_families"] for v in ins),
                   n_deletion_events=sum(v["n_supporting_families"] for v in dels),
                   n_inserted_bases=n_ins, n_deleted_bases=n_del,
                   n_families_with_variants=len(affected),
                   n_multifamily_variants=sum(v["n_supporting_families"] >= 2 for v in details),
                   n_multisample_variants=sum(v["n_supporting_samples"] >= 2 for v in details),
                   n_multitimepoint_variants=sum(v["n_supporting_timepoints"] >= 2 for v in details),
                   max_variant_frequency=max((v["variant_frequency"] for v in details), default=0)
                   if evaluable_families else "NA",
                   substitution_sites_per_kb=ratio(len(sub_sites), len(eligible), 1000),
                   substitutions_per_kb=ratio(n_sub, family_bases, 1000),
                   inserted_bases_per_kb=ratio(n_ins, family_bases, 1000),
                   deleted_bases_per_kb=ratio(n_del, family_bases, 1000),
                   variant_bases_per_kb=ratio(n_sub + n_ins + n_del, family_bases, 1000),
                   percent_families_with_variants=ratio(len(affected), len(evaluable_families), 100))
    plot(row, ref, reads, calls, info, depth, blocks, insertion_call)
    return summary, details, longitudinal


def plot(row, ref, reads, calls, info, depth, blocks, insertion_call):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {"covered": "#9ca3af", "substitution": "#c83239",
              "insertion": "#2563eb", "deletion": "#d97706", "masked": "#dedede"}
    reads = sorted(reads, key=lambda r: (time_order(info[r[0]]), r[0], r[2]))
    fig, (ax, freq, cov) = plt.subplots(3, 1, sharex=True,
        figsize=(13, min(50, max(6, 4 + len(reads) * 0.20))),
        gridspec_kw={"height_ratios": [max(2, len(reads) * 0.17), 1.4, 1.4]})
    width = len(ref)
    insertion_cols = {i for a, b in blocks for i in range(a, b + 1)}
    masked = [i for i, c in enumerate(ref) if c not in STATES]
    for a, b in runs(masked):
        ax.axvspan(a + .5, b + 1.5, color=colors["masked"], alpha=.6)
    for y, (key, s, name) in enumerate(reads):
        for a, b in runs([i for i, v in enumerate(s) if v in STATES]):
            ax.hlines(y, a + 1, b + 1, color=colors["covered"], linewidth=2)
        for kind, indices in [
            ("substitution", [i for i, v in enumerate(s) if ref[i] in BASES and v in BASES and v != ref[i]]),
            ("deletion", [i for i, v in enumerate(s) if ref[i] in BASES and v == "-"]),
            ("insertion", [i for i in sorted(insertion_cols) if s[i] in BASES])]:
            if indices:
                ax.scatter([i + 1 for i in indices], [y] * len(indices),
                           color=colors[kind], marker="|", s=28, linewidths=1.3)
    ax.set_ylim(len(reads) - .5, -.5)
    if len(reads) <= 80:
        ax.set_yticks(range(len(reads)), [f"{info[k]} | {k[0]} | {k[1]} | {i+1}"
                                        for i, (k, _, _) in enumerate(reads)], fontsize=6)
    else:
        ax.set_ylabel("Reads (ordered by timepoint)")
    ax.set_title(row["consensus_header"], loc="left", fontsize=10)
    ax.legend(handles=[Line2D([0], [0], color=c, label=k, lw=3)
                       for k, c in colors.items()], loc="upper center",
              bbox_to_anchor=(.5, 1.20), ncol=5, fontsize=8, frameon=False)
    for kind, color in colors.items():
        if kind not in ("substitution", "deletion", "insertion"):
            continue
        xs, ys = [], []
        for i in range(width):
            if kind == "insertion":
                if i not in insertion_cols:
                    continue
                a, b = next((a, b) for a, b in blocks if a <= i <= b)
                covered = [s for s in calls.values() if insertion_call(s, a, b) is not None]
                n = sum(s[i] in BASES for s in covered)
            else:
                if ref[i] not in BASES:
                    continue
                covered = [s for s in calls.values() if s[i] in STATES]
                n = sum((s[i] == "-" if kind == "deletion" else
                         s[i] in BASES and s[i] != ref[i]) for s in covered)
            if n and len(covered) >= 2:
                xs.append(i + 1)
                ys.append(n / len(covered))
        freq.vlines(xs, 0, ys, color=color, linewidth=1)
    freq.set_ylim(0, 1)
    freq.set_ylabel("Family variant\nfrequency")
    cov.fill_between(range(1, width + 1), [depth[i] if ref[i] in BASES else 0
                                          for i in range(width)], color="#25856b", alpha=.6)
    cov.set_ylabel("Family depth")
    cov.set_xlabel("Alignment column (1-based)")
    cov.set_xlim(.5, width + .5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(row["consensus_header"] + "." + ext, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["analyze", "merge"])
    p.add_argument("--metadata")
    p.add_argument("--alignment")
    p.add_argument("--input")
    p.add_argument("--family-fraction", type=float, default=.60)
    args = p.parse_args()
    if args.mode == "analyze":
        data = json.loads(Path(args.metadata).read_text())
        row = data["clone"]
        Path("report").mkdir()
        summary, details, time, skipped = [], [], [], []
        try:
            if data["missing_alignment"]:
                raise ValueError("Alignment file not found")
            s, details, time = analyze(row, data["assignments"], args.alignment, args.family_fraction)
            summary = [s]
        except ValueError as e:
            skipped = [dict(row, skip_reason=str(e))]
        write_csv("report/summary.csv", list(row) + COUNTS, summary)
        write_csv("report/variants.csv", DETAIL, details)
        write_csv("report/timepoints.csv", TIME, time)
        write_csv("report/skipped.csv", list(row) + ["skip_reason"], skipped)
    else:
        with open(args.input, encoding="utf-8-sig", newline="") as f:
            fields = csv.DictReader(f).fieldnames
        for source, target, cols in [
            ("summary.csv", "clone_variation_summary.csv", fields + COUNTS),
            ("variants.csv", "variants.csv", DETAIL),
            ("timepoints.csv", "variants_by_timepoint.csv", TIME),
            ("skipped.csv", "skipped_clones.csv", fields + ["skip_reason"])]:
            rows = [row for path in sorted(Path(".").glob("clone*/report/" + source))
                    for row in read_csv(path)]
            rows.sort(key=lambda r: tuple(r.get(k, "") for k in IDS) +
                      (int(r.get("alignment_start", 0)), r.get("variant_type", ""),
                       r.get("alternative_allele", ""), time_order(r.get("timepoint", ""))))
            write_csv(target, cols, rows)


if __name__ == "__main__":
    main()
