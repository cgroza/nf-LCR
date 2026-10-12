#!/usr/bin/env python3
"""Extract and CHM13-strand-normalize chr22q haplotypes using LCR flanks."""

import argparse
import csv
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LCRS = {
    "A": (846123, 1410825),
    "B": (2643128, 2766966),
    "C": (3076285, 3146770),
    "D": (3520158, 3976302),
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, header, sequence = [], None, []
    with path.open() as fasta:
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header, sequence = line[1:].split()[0], []
            else:
                sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as output:
        for name, sequence in records:
            output.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                output.write(f"{sequence[start:start + 80]}\n")


def paf_hits(path: Path) -> dict[tuple[str, str], tuple[float, list[str]]]:
    best: dict[tuple[str, str], tuple[float, list[str]]] = {}
    with path.open() as paf:
        for line in paf:
            fields = line.rstrip("\n").split("\t")
            query_length, query_start, query_end = map(int, fields[1:4])
            matches, alignment_length = map(int, fields[9:11])
            coverage = (query_end - query_start) / query_length
            identity = matches / alignment_length if alignment_length else 0
            score = coverage * identity
            key = (fields[0], fields[5])
            if key not in best or score > best[key][0]:
                best[key] = (score, fields)
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-roi", type=Path, required=True)
    parser.add_argument("--assemblies", type=Path, nargs="+", required=True)
    parser.add_argument("--flank-paf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--flank-bp", type=int, default=30000)
    parser.add_argument("--min-coverage", type=float, default=0.90)
    parser.add_argument("--min-identity", type=float, default=0.90)
    args = parser.parse_args()

    reference = read_fasta(args.reference_roi)
    if len(reference) != 1:
        raise SystemExit("--reference-roi must contain exactly one 18-23 Mb sequence.")
    reference_name, reference_sequence = reference[0]
    haplotypes: list[tuple[str, str]] = []
    names = {reference_name}
    for assembly in args.assemblies:
        for name, sequence in read_fasta(assembly):
            if name in names:
                raise SystemExit(f"Duplicate FASTA record name: {name}")
            names.add(name)
            haplotypes.append((name, sequence))
    if not haplotypes:
        raise SystemExit("No input FASTA records were found.")

    hits = paf_hits(args.flank_paf)
    lengths = dict(haplotypes)
    selected: list[tuple[str, str]] = [(reference_name, reference_sequence)]
    orientations = [{
        "haplotype": reference_name,
        "orientation_relative_to_chm13": "+",
        "reverse_complement": "false",
        "status": "reference",
    }]
    windows = {lcr: [{
        "haplotype": reference_name,
        "lcr_start_0based": start,
        "lcr_end_0based": end,
    }] for lcr, (start, end) in LCRS.items()}
    anchor_rows = []
    excluded = []
    boundary_rows = [
        {
            "haplotype": reference_name,
            "lcr": f"LCR-{lcr}",
            "boundary_source": "reference",
        }
        for lcr in LCRS
    ]

    for name, sequence in haplotypes:
        chosen = {}
        for lcr in LCRS:
            for side in ("left", "right"):
                query = f"{lcr}_{side}"
                record = hits.get((query, name))
                if record is None:
                    continue
                _, fields = record
                qlen, qstart, qend = map(int, fields[1:4])
                matches, alen = map(int, fields[9:11])
                coverage = (qend - qstart) / qlen
                identity = matches / alen if alen else 0
                if coverage >= args.min_coverage and identity >= args.min_identity:
                    chosen[query] = fields
        anchored_lcrs = {
            lcr
            for lcr in LCRS
            if any(f"{lcr}_{side}" in chosen for side in ("left", "right"))
        }
        if len(anchored_lcrs) < 2:
            excluded.append({
                "haplotype": name,
                "reason": "fewer_than_two_localized_lcrs",
                "anchors_found": ",".join(sorted(chosen)),
            })
            continue
        strands = {fields[4] for fields in chosen.values()}
        if len(strands) != 1:
            excluded.append({
                "haplotype": name,
                "reason": "inconsistent_anchor_orientations",
                "anchors_found": ",".join(sorted(chosen)),
            })
            continue
        reverse = strands == {"-"}
        length = len(sequence)
        transformed_windows = {}
        for lcr, (reference_start, reference_end) in LCRS.items():
            left, right = chosen.get(f"{lcr}_left"), chosen.get(f"{lcr}_right")
            if left is None and right is None:
                boundary_rows.append({
                    "haplotype": name,
                    "lcr": f"LCR-{lcr}",
                    "boundary_source": "not_localized",
                })
                continue
            expected_length = reference_end - reference_start
            if left is not None and right is not None:
                if not reverse:
                    start, end = int(left[8]), int(right[7])
                else:
                    start, end = length - int(left[8]), length - int(right[7])
                boundary_source = "both_flanks"
            elif left is not None:
                start = int(left[8]) if not reverse else length - int(left[8])
                end = start + expected_length
                boundary_source = "left_flank_plus_chm13_length"
            else:
                end = int(right[7]) if not reverse else length - int(right[7])
                start = end - expected_length
                boundary_source = "right_flank_minus_chm13_length"
            if start < 0 or end > length or start >= end:
                boundary_rows.append({
                    "haplotype": name,
                    "lcr": f"LCR-{lcr}",
                    "boundary_source": "inferred_boundary_out_of_bounds",
                })
                continue
            transformed_windows[lcr] = (start, end)
            boundary_rows.append({
                "haplotype": name,
                "lcr": f"LCR-{lcr}",
                "boundary_source": boundary_source,
            })
        normalized = sequence.translate(COMPLEMENT)[::-1] if reverse else sequence
        selected.append((name, normalized))
        orientations.append({
            "haplotype": name,
            "orientation_relative_to_chm13": "-" if reverse else "+",
            "reverse_complement": str(reverse).lower(),
            "status": "selected",
        })
        for lcr, (start, end) in transformed_windows.items():
            windows[lcr].append({
                "haplotype": name,
                "lcr_start_0based": start,
                "lcr_end_0based": end,
            })
        for query, fields in chosen.items():
            anchor_rows.append({
                "haplotype": name,
                "anchor": query,
                "source_strand": fields[4],
                "target_start_0based": fields[7],
                "target_end_0based": fields[8],
            })

    if len(selected) == 1:
        raise SystemExit(
            "No input records localized at least two LCRs with consistent flank-anchor "
            "orientation. Provide chr22q haplotypes containing at least two LCRs."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(args.output_dir / "normalized_chr22q_haplotypes.fa", selected)
    with (args.output_dir / "orientations.tsv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(orientations[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(orientations)
    for lcr, rows in windows.items():
        with (args.output_dir / f"lcr{lcr.lower()}_display_windows.tsv").open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    with (args.output_dir / "flank_anchor_placements.tsv").open("w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=[
                "haplotype", "anchor", "source_strand",
                "target_start_0based", "target_end_0based",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(anchor_rows)
    with (args.output_dir / "lcr_boundary_sources.tsv").open("w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=["haplotype", "lcr", "boundary_source"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(boundary_rows)
    with (args.output_dir / "excluded_haplotypes.tsv").open("w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=["haplotype", "reason", "anchors_found"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(excluded)


if __name__ == "__main__":
    main()
