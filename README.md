# MethPhaser

MethPhaser is now a single-command Python package.

`methphaser` performs the full pipeline end-to-end:
1. infer cross-block relationships from methylation
2. rewrite phased VCF and BAM
3. merge per-chromosome BAMs into one final BAM (automatic)

## Install

```bash
python -m pip install -e .
```

## Usage

```bash
methphaser \
  --bam /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --output-dir /path/to/output
```

Optional inputs:

- `--gtf`: haploblocks GTF (auto-discovered when omitted)
- `--vcf`: phased SNP VCF (auto-discovered when omitted)
- `--threads`: worker count
- `--region-strategy`: `full-span` (default) or `boundary-window`
- `--boundary-window-bp`: flank size for `boundary-window` mode (default `1000`)

## Outputs

Inside `--output-dir`:

- `methphaser.vcf`
- `methphaser.bam`
- `methphaser.bam.bai`
- `phaseblock_parent_assignment.tsv` (phaseblock-level father/mother assignment status)
- `work/` (internal intermediate relationship/read-assignment files)

## Built-in Imprinting Database

MethPhaser now bundles a genome-wide imprinting BED database at:

- `src/methphaser/data/imprinting_regions.hg38.bed`

Current file is built from four HumanICR tracks (`known_ICR_25`,
`putative_ICRs_35_65_final_N1488`, `ICRs_2_or_more_method_val_N1348`,
`ICRs_2_or_more_method_N4239`) and merged into one repository dataset.
It also includes inferred parent-of-origin methylation direction
(`methylated_parent = maternal|paternal|unknown`) using HumanICR gamete tracks
(`Oocyte_100`, `Sperm_with_sue_100`).
Source links are documented in:

- `src/methphaser/data/imprinting_regions.hg38.SOURCES.md`

To refresh the database from source tracks:

```bash
python tools/build_imprinting_database.py
```

## Parent Assignment Logic

For each phaseblock in `methphaser.vcf`, MethPhaser:

1. finds overlapping directional imprinting regions (`maternal`/`paternal` methylated),
2. measures H1/H2 methylation in those regions from `methphaser.bam`,
3. votes whether H1 or H2 is paternal based on expected methylation direction,
4. reports `assigned` only when evidence is decisive; otherwise reports `random`.

This avoids forcing `H1=father` when directional methylation evidence is absent.

## Behavior Simplifications

- Truth-VCF benchmarking was removed.
- `full-span` keeps legacy Step1 assignment behavior.
- `boundary-window` is available as a faster experimental mode.
- Final BAM merge is always performed.

## Test Subset Data

A local subset dataset is bundled under `tests/data`:

- `tests/data/HG002_2025q1_PAW70337_30x.chr21_4.9M_6.2M.subset.bam`
- `tests/data/reference/GRCh38.chr21.fa`
- `tests/data/PAW70337_output/SAMPLE.wf_snp.haploblocks.gtf`
- `tests/data/PAW70337_output/SAMPLE.wf_snp.vcf.gz`
