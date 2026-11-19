# chr22q family workflow

## SLURM Nextflow workflow

`main.nf` creates one analysis input per family, aligns each family, and
identifies chr22q intervals:

1. Reads a CSV sample sheet with exactly the required `ASSEMBLY,FAMILY`
   columns.
2. Decompresses the supplied full CHM13 genome when gzipped and concatenates
   it with all assembly FASTAs from a family.
3. Runs:

   ```bash
   wfmash -t 64 -p 95 -Y '#' -S 20k -j 200k --scaffold-out
   ```

   and writes the PAF alignment.
4. Retrieves homologous chr22q intervals using:

   ```bash
   impg query -a ALIGNMENTS.paf \
     -r NC_060946.1:18000000-23000000 -d 1000000 |
     awk '$3-$2 > 100e3 {print $1":"$2+1"-"$3}'
   ```

5. Runs `samtools faidx -r` to extract those intervals from the corresponding
   family FASTA into a chr22q-only FASTA.
6. Renames the extracted CHM13 chromosome 22 interval to
   `CHM13#0#chr22:18000001-23000000`, then runs `run_chr22q_lcr_analysis.sh`
   using the supplied HOR consensus FASTA and CHM13 GFF.

The default Nextflow configuration runs every process in
`cgroza/nf-chr22q:latest` using Docker. Docker must be available on the SLURM
compute nodes; the image must be pullable there before jobs start.

The sample-sheet `ASSEMBLY` values are paths to individual FASTA files. They
may be absolute paths or paths relative to the directory from which Nextflow
is launched. For example:

```csv
ASSEMBLY,FAMILY
/data/family_01/proband.fa,family_01
/data/family_01/mother.fa,family_01
/data/family_02/proband.fa,family_02
```

```bash
nextflow run main.nf -profile slurm \
  --samplesheet families.csv \
  --chm13_fasta /path/to/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
  --hor_consensuses /path/to/22q11.2_25-28-41-105-120_HORs.fasta \
  --annotation /path/to/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz \
  --outdir results
```

Results are published into:

- `results/family_fastas/FAMILY/family.fa`
- `results/alignments/FAMILY/alignments.paf`
- `results/chr22q_intervals/FAMILY/chr22q_intervals.txt`
- `results/chr22q_fastas/FAMILY/chr22q.fa`
- `results/lcr_analysis/FAMILY/chr22q_lcr_analysis/`

The default query target is the NCBI RefSeq chromosome 22 record
`NC_060946.1:18000000-23000000`. The LCR analysis reference record is named
`CHM13#0#chr22:18000001-23000000`.

```bash
--reference_record 'CHM13#0#chr22:18000000-23000000'
```

The wfmash task requests 64 CPUs, 128 GB RAM, and 24 hours. The LCR-analysis
task requests 8 CPUs, 32 GB RAM, and 8 hours; extraction tasks request 4
CPUs, 16 GB RAM, and 8 hours.

## Container environment

`Dockerfile` uses Ubuntu and compiles `minimap2`, `samtools`, and IMPG
(including `wfmash`) from source. It also installs Nextflow, OpenJDK, R, and
builds `ggplot2` from CRAN source. The image copies only the workflow and
scripts; mount input FASTAs, the sample sheet, the HOR consensus FASTA, and
the CHM13 GFF at runtime.

```bash
./build.sh
docker run --rm -v "$PWD:/work" cgroza/nf-chr22q:latest \
  nextflow run /opt/chr22q/main.nf \
  --samplesheet /work/families.csv \
  --chm13_fasta /work/CHM13.fa \
  --hor_consensuses /work/22q11.2_25-28-41-105-120_HORs.fasta \
  --annotation /work/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz \
  --outdir /work/results
```

`build.sh` creates a temporary PEM bundle from the macOS system keychains and
passes it only as the BuildKit `mac_ca` secret. TLS verification remains
enabled; the certificate bundle is removed after the build and is not copied
into the final image. Pass an alternate image tag as its sole argument:

```bash
./build.sh cgroza/nf-chr22q:latest
```

For SLURM execution, run Nextflow on the cluster host and ensure the compute
nodes can access Docker and pull `cgroza/nf-chr22q:latest`.
