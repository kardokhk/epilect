# Epilect

**Assembly-free, per-read resolution of co-resident bacterial strains by DNA-methylation fingerprinting.**

Closely related bacterial strains can be near-identical at the sequence level
(more than 99% average nucleotide identity) yet carry distinct DNA-methylation
"fingerprints" laid down by their restriction-modification systems. Epilect reads
those per-read methylation patterns directly off nanopore modified-base BAMs and
uses them to:

1. **classify** each read to its strain of origin by a per-motif likelihood-ratio test,
2. **phase** heterozygous sites genome-wide by orienting SNP-phasing blocks with the
   per-read methylation labels (joint methylation + SNP phasing), and
3. **predict resolvability** - for any conspecific strain pair, whether methylation
   alone can separate them at the single-read level.

Classification needs no genome assembly, no variant calling, and no deep
resequencing: the only required input is a standard Dorado modified-base BAM.

The numerical core is a small set of pure, unit-tested functions (`epilect.core`);
the command-line interface (`epilect.cli`) wires them onto real BAMs via `modkit`,
`samtools`, `bcftools`, and `whatshap`.

## Requirements

- Linux x86_64 (tested), Python 3.9 or newer.
- CPU only. A GPU is needed only for upstream Dorado basecalling / modified-base
  calling, which produces the input BAM and is out of scope for this package.
- External command-line tools (installed via the provided conda environment):
  `modkit`, `samtools`, `bcftools`, `whatshap`, and (for self-assembly workflows)
  `minimap2`, `flye`, `nanomotif`.

## Installation

Expected install time is about 5-10 minutes.

```bash
# 1. Get the code
git clone https://github.com/kardokhk/epilect.git
cd epilect

# 2. Create the external-tools environment (pinned versions; see environment.yml)
conda env create -f environment.yml
conda activate epilect

# 3. Install the Python package (editable install; add dev extras for tests)
pip install -e ".[dev]"
```

Epilect finds `modkit` via the `MODKIT` environment variable if set, otherwise on
your `PATH` (the conda environment provides it). `samtools`, `bcftools`, and
`whatshap` are taken from your `PATH`.

## Quickstart

```bash
epilect --help          # top-level help
epilect classify --help # per-subcommand help

make test               # run the unit test suite
make demo               # classify the bundled example data and write a figure panel
```

`make demo` runs `epilect classify` on a small bundled example
(`demo/sdse_demo.bam`, about 1,400 reads, a 1:8 minor:major mix of two
Streptococcus dysgalactiae strains) and writes `demo/demo_panel.png`, the per-read
score distribution coloured by inferred strain. It also prints the expected vs
observed minor-strain fraction. Runtime is under about 2 minutes on a single CPU.

## Usage

### `classify` - per-read strain assignment

Extracts per-read modified-base calls from a modified-base BAM with
`modkit extract calls`, splits them into motif channels, and assigns each read to
one of two strains by the log-likelihood ratio of its methylation pattern.

```bash
epilect classify \
    --bam sample.modbam.bam \
    --reference ref.fna \
    --motifs "GATC:1:A:a,CCWGG:0:C:m" \
    --emission model.json \
    --out reads.tsv
```

- `--motifs` is a comma-separated list of `MOTIF:OFFSET:CANONICAL_BASE:MOD_CODE`
  (default: the Streptococcus dysgalactiae model `GATC:1:A:a,CCWGG:0:C:m`).
- `--emission model.json` gives per-strain methylation rates, for example
  `{"StrainA": {"GATC": 0.9, ...}, "StrainB": {...}}`. Alternatively pass
  `--estimate` to derive a two-cluster model directly from the data, or omit both
  to use the built-in model.
- Output is a TSV with columns: `read_id  strain  llr  n_informative_sites`.

### `phase` - methylation-only or joint methylation + SNP phasing

```bash
# Joint: whatshap phase blocks oriented genome-wide by per-read methylation labels
epilect phase --mode joint \
    --bam sample.bam --reference ref.fna \
    --vcf het.vcf.gz --labels reads.tsv --out phased.tsv

# Methylation-only baseline (no whatshap / SNP linkage)
epilect phase --mode meth-only \
    --bam sample.bam --vcf het.vcf.gz --labels reads.tsv --out phased.tsv
```

`--labels` is the `classify` output (or any `read_id<TAB>strain` TSV). Output is a
per-site TSV giving the inferred ALT-carrying strain.

### `resolvable` - a-priori resolvability tiering

```bash
# Tier every conspecific pair in a REBASE-style motif table
epilect resolvable --motif-table genome_motifs.tsv --min-motifs 3 --out tiers.json

# Tier a single pair from two motif-set files (one motif per line)
epilect resolvable --pair strainA.motifs strainB.motifs
```

Tiers:
- `T1_bidirectional` - each strain has a private motif; both are recoverable from single reads,
- `T1_onesided` - only one strain has a private motif,
- `T2_coverage` - strains differ but only via sparse motifs (needs coverage),
- `T0_none` - identical methylated-motif sets, so the pair is unresolvable by methylation.

## Input: how to produce a modified-base BAM

Epilect consumes a standard nanopore modified-base BAM (MM/ML tags). Basecall
R10.4.1 raw signal with Dorado using a modified-base model, then map with
`minimap2 -ax map-ont -y` (the `-y` flag propagates the modification tags) and
sort/index with `samtools`. No special library preparation is required.

## Tests

```bash
make test        # or: pytest -q
```

The suite covers the core numerics (likelihood-ratio classifier, multi-strain
assignment, and motif-density tiering).

## Citation

If you use Epilect, please cite the accompanying manuscript (in preparation).
Citation details will be added here on publication.

## License

MIT. See [LICENSE](LICENSE).
