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
- `work/` (internal intermediate relationship/read-assignment files)

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
