#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_chr22q_lcr_analysis.sh \
    --fasta assemblies.fa \
    --hor-consensuses 22q11.2_25-28-41-105-120_HORs.fasta \
    --annotation GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz \
    --output results

Optional:
  --reference-record CHM13#0#chr22:18000000-23000000

The FASTA must contain a CHM13 chr22 record and additional PanSN-formatted
haplotype contigs. The expected CHM13 record name is
CHM13#0#chr22:18000001-23000000, representing the already subset 18-23 Mb
interval. Use --reference-record only when using a different identifier.

Non-reference records are retained when 30-kb CHM13 flanks localize at least
two of LCR-A/B/C/D with
>=90% coverage and >=90% identity. Input records must therefore be contiguous
chr22q haplotypes spanning at least two LCRs. LCRs with no flanking anchor are
shown as absent rather than causing the entire haplotype to be dropped.

Outputs are written to OUTPUT/strand_normalized/. Existing output directories
are never overwritten.

Requirements: bash, python3, minimap2, Rscript, and R package ggplot2.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fasta=
reference_record=
hor_consensuses=
annotation=
output=

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) fasta=$2; shift 2 ;;
    --reference-record) reference_record=$2; shift 2 ;;
    --hor-consensuses) hor_consensuses=$2; shift 2 ;;
    --annotation) annotation=$2; shift 2 ;;
    --output) output=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --*) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    *) break ;;
  esac
done

[[ -n "$fasta" && -n "$hor_consensuses" && -n "$annotation" && -n "$output" && $# -eq 0 ]] || {
  usage >&2
  exit 2
}

for command in python3 minimap2 Rscript; do
  command -v "$command" >/dev/null || { printf 'Missing required command: %s\n' "$command" >&2; exit 1; }
done
Rscript -e 'quit(status = if (requireNamespace("ggplot2", quietly = TRUE)) 0 else 1)' ||
  { printf 'Missing required R package: ggplot2\n' >&2; exit 1; }
for path in "$fasta" "$hor_consensuses" "$annotation"; do
  [[ -f "$path" ]] || { printf 'Input file does not exist: %s\n' "$path" >&2; exit 1; }
done
[[ ! -e "$output" ]] || { printf 'Refusing to overwrite existing output: %s\n' "$output" >&2; exit 1; }

mkdir -p "$output/strand_normalized"
work="$output/strand_normalized"
trap 'printf "Analysis failed; partial results retained in %s\n" "$output" >&2' ERR

python3 - "$fasta" "$reference_record" "$work/chm13_chr22_18m_23m.fa" "$work/input_assemblies.fa" <<'PY'
import sys
fasta, reference_name, reference_output, assemblies_output = sys.argv[1:]
records, header, sequence = [], None, []
for line in open(fasta):
    line = line.strip()
    if line.startswith(">"):
        if header is not None:
            records.append((header, "".join(sequence)))
        header, sequence = line[1:].split()[0], []
    else:
        sequence.append(line)
if header is not None:
    records.append((header, "".join(sequence)))
if reference_name:
    reference_records = [(name, sequence) for name, sequence in records if name == reference_name]
else:
    reference_records = [
        (name, sequence)
        for name, sequence in records
        if name == "CHM13#0#chr22:18000001-23000000"
    ]
if len(reference_records) != 1:
    raise SystemExit(
        "Could not identify CHM13#0#chr22:18000001-23000000. "
        "Provide --reference-record with its exact PanSN name."
    )
reference_name, reference_sequence = reference_records[0]
roi = reference_sequence
with open(reference_output, "w") as output:
    output.write(f">{reference_name}\n")
    for start in range(0, len(roi), 80):
        output.write(f"{roi[start:start + 80]}\n")
with open(assemblies_output, "w") as output:
    for name, sequence in records:
        if name == reference_name:
            continue
        output.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            output.write(f"{sequence[start:start + 80]}\n")
PY

python3 - "$work/chm13_chr22_18m_23m.fa" "$work/chm13_lcr_flanks.fa" <<'PY'
import sys
reference, output = sys.argv[1:]
lines = [line.strip() for line in open(reference) if not line.startswith(">")]
sequence = "".join(lines)
lcrs = {"A": (846123, 1410825), "B": (2643128, 2766966), "C": (3076285, 3146770), "D": (3520158, 3976302)}
with open(output, "w") as out:
    for lcr, (start, end) in lcrs.items():
        for side, begin, finish in (("left", start - 30000, start), ("right", end, end + 30000)):
            out.write(f">{lcr}_{side}\n{sequence[begin:finish]}\n")
PY

minimap2 -x asm20 -c --eqx --secondary=yes -p 0.50 -N 50 \
  "$work/input_assemblies.fa" "$work/chm13_lcr_flanks.fa" > "$work/lcr_flanks_to_inputs.paf"

python3 "$script_dir/prepare_chr22q_haplotypes.py" \
  --reference-roi "$work/chm13_chr22_18m_23m.fa" \
  --assemblies "$work/input_assemblies.fa" \
  --flank-paf "$work/lcr_flanks_to_inputs.paf" \
  --output-dir "$work"
if [[ $(wc -l < "$work/excluded_haplotypes.tsv") -gt 1 ]]; then
  printf 'Warning: some input records were excluded; see %s\n' "$work/excluded_haplotypes.tsv" >&2
fi

python3 "$script_dir/extract_chr22q_protein_coding_genes.py" \
  "$annotation" "$work/chm13_chr22_18m_23m.fa" \
  "$work/chm13_chr22q_protein_coding_genes.fa" \
  "$work/chm13_chr22q_protein_coding_genes.tsv"

minimap2 -x asm20 -c --eqx --secondary=yes -E1,0 -p 0.25 -N 100 \
  "$work/normalized_chr22q_haplotypes.fa" "$hor_consensuses" > "$work/hor_to_haplotypes.sensitive.paf"
minimap2 -x asm20 -c --eqx --secondary=yes -p 0.50 -N 20 \
  "$work/normalized_chr22q_haplotypes.fa" "$work/chm13_chr22q_protein_coding_genes.fa" > "$work/chm13_protein_coding_genes_to_haplotypes.paf"

(
  cd "$output"
  Rscript "$script_dir/plot_chr22q_full_lcr_layout.R" \
    strand_normalized/hor_to_haplotypes.sensitive.paf \
    strand_normalized \
    strand_normalized/full_chr22q_lcr_layout
  Rscript "$script_dir/plot_chr22q_full_lcr_gene_layout.R"
  Rscript "$script_dir/plot_chr22q_lcr_gene_configurations.R"
)

printf 'Completed chr22q LCR analysis: %s\n' "$work"
