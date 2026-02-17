from __future__ import annotations

import argparse
import sys

from .pipeline import PipelineConfig, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="methphaser: one-command methylation-assisted phasing"
    )
    parser.add_argument("-b", "--bam", required=True, help="haplotagged BAM file")
    parser.add_argument("-r", "--reference", required=True, help="reference FASTA")
    parser.add_argument(
        "-g",
        "--gtf",
        default=None,
        help="haploblocks GTF file (auto-discovered when omitted)",
    )
    parser.add_argument(
        "-v",
        "--vcf",
        default=None,
        help="phased SNP VCF file (auto-discovered when omitted)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="methphaser_output",
        help="output directory",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="threads")
    parser.add_argument(
        "--region-strategy",
        choices=["boundary-window", "full-span"],
        default="full-span",
        help="phase-region strategy for step1 read assignment",
    )
    parser.add_argument(
        "--boundary-window-bp",
        type=int,
        default=1000,
        help="window size for boundary-window strategy",
    )
    return parser


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    args = build_parser().parse_args(argv)

    outputs = run_pipeline(
        PipelineConfig(
            bam=args.bam,
            reference=args.reference,
            output_dir=args.output_dir,
            gtf=args.gtf,
            vcf=args.vcf,
            threads=max(1, args.threads),
            region_strategy=args.region_strategy,
            boundary_window_bp=max(1, args.boundary_window_bp),
        )
    )

    print("methphaser finished")
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
