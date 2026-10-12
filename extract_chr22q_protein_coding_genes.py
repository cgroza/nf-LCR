#!/usr/bin/env python3
"""Extract CHM13 protein-coding gene loci from a chr22q interval."""

import argparse
import csv
import gzip
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    sequence: list[str] = []
    with path.open() as fasta:
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence)
    return records


def attributes(text: str) -> dict[str, str]:
    return {
        key: value
        for field in text.split(";")
        if "=" in field
        for key, value in [field.split("=", 1)]
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotation_gff", type=Path)
    parser.add_argument("roi_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--chromosome", default="NC_060946.1")
    parser.add_argument("--roi-start-1based", type=int, default=18_000_000)
    parser.add_argument("--roi-end-1based", type=int, default=23_000_000)
    args = parser.parse_args()

    sequences = read_fasta(args.roi_fasta)
    if len(sequences) != 1:
        raise SystemExit("The ROI FASTA must contain exactly one CHM13 sequence.")
    roi_sequence = next(iter(sequences.values()))
    genes: list[dict[str, str | int]] = []

    with gzip.open(args.annotation_gff, "rt") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != args.chromosome or fields[2] != "gene":
                continue
            start_1based, end_1based = int(fields[3]), int(fields[4])
            if start_1based > args.roi_end_1based or end_1based < args.roi_start_1based:
                continue
            attrs = attributes(fields[8])
            if attrs.get("gene_biotype") != "protein_coding":
                continue
            start_0based = start_1based - args.roi_start_1based
            end_0based = end_1based - args.roi_start_1based + 1
            genes.append(
                {
                    "gene": attrs["gene"],
                    "gene_id": attrs["Dbxref"].split(":")[-1],
                    "start_1based": start_1based,
                    "end_1based": end_1based,
                    "strand": fields[6],
                    "roi_start_0based": start_0based,
                    "roi_end_0based": end_0based,
                    "length_bp": end_0based - start_0based,
                }
            )

    genes.sort(key=lambda gene: int(gene["roi_start_0based"]))
    with args.output_tsv.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(genes[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(genes)
    with args.output_fasta.open("w") as output:
        for gene in genes:
            output.write(
                f">{gene['gene']}|{gene['start_1based']}-{gene['end_1based']}|{gene['strand']}\n"
            )
            sequence = roi_sequence[
                int(gene["roi_start_0based"]):int(gene["roi_end_0based"])
            ]
            for start in range(0, len(sequence), 80):
                output.write(f"{sequence[start:start + 80]}\n")


if __name__ == "__main__":
    main()
