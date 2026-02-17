from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
from typing import Dict, List

from .discovery import PhaseInputs, discover_phase_inputs
from .imprinting import ParentAssignmentConfig, run_parent_assignment
from .step1_parallel import Step1Config, run_step1
from .step2_postprocess import Step2Config, run_step2


@dataclass(frozen=True)
class PipelineConfig:
    bam: str
    reference: str
    output_dir: str = "methphaser_output"
    gtf: str | None = None
    vcf: str | None = None
    threads: int = 1
    region_strategy: str = "full-span"
    boundary_window_bp: int = 1000


def _chrom_sort_key(path: Path) -> tuple[int, str]:
    parts = path.name.split(".")
    chrom = parts[-3] if len(parts) >= 4 else ""
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return (int(suffix), chrom)
        if suffix == "X":
            return (23, chrom)
        if suffix == "Y":
            return (24, chrom)
        if suffix in {"M", "MT"}:
            return (25, chrom)
    return (99, chrom)


def _merge_chromosome_bams(chrom_bams: List[str], final_bam: Path, threads: int) -> None:
    if not chrom_bams:
        raise RuntimeError("no chromosome BAM outputs found to merge")

    bam_paths = sorted((Path(path) for path in chrom_bams), key=_chrom_sort_key)
    unsorted_bam = final_bam.with_suffix(".unsorted.bam")
    if len(bam_paths) == 1:
        shutil.copyfile(bam_paths[0], unsorted_bam)
    else:
        cmd = ["samtools", "merge", "-f", "-@", str(max(1, threads)), str(unsorted_bam)]
        cmd.extend(str(path) for path in bam_paths)
        subprocess.run(cmd, check=True)

    subprocess.run(
        [
            "samtools",
            "sort",
            "-@",
            str(max(1, threads)),
            "-o",
            str(final_bam),
            str(unsorted_bam),
        ],
        check=True,
    )
    unsorted_bam.unlink(missing_ok=True)
    subprocess.run(["samtools", "index", str(final_bam)], check=True)


def run_pipeline(config: PipelineConfig) -> Dict[str, str]:
    output_dir = Path(config.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    discovered: PhaseInputs = discover_phase_inputs(
        bam_file=config.bam,
        gtf=config.gtf,
        vcf=config.vcf,
    )

    processed_chromosomes = run_step1(
        Step1Config(
            bam_file=str(Path(config.bam).resolve()),
            reference=str(Path(config.reference).resolve()),
            gtf=str(discovered.gtf),
            output_dir=str(work_dir),
            threads=max(1, config.threads),
            region_strategy=config.region_strategy,
            boundary_window_bp=max(1, config.boundary_window_bp),
        )
    )

    output_vcf = output_dir / "methphaser.vcf"
    output_bam_prefix = output_dir / "methphaser"
    chrom_bams = run_step2(
        Step2Config(
            input_bam_file=str(Path(config.bam).resolve()),
            meth_phasing_input_folder=str(work_dir),
            output_vcf=str(output_vcf),
            output_bam_prefix=str(output_bam_prefix),
            vcf_called=str(discovered.vcf),
            threads=max(1, config.threads),
            chromosomes=processed_chromosomes,
        )
    )

    merged_bam = output_dir / "methphaser.bam"
    _merge_chromosome_bams(chrom_bams, merged_bam, max(1, config.threads))

    parent_outputs = run_parent_assignment(
        ParentAssignmentConfig(
            output_vcf=str(output_vcf),
            phaseblock_assignment_output=str(output_dir / "phaseblock_parent_assignment.tsv"),
            input_bam=str(merged_bam),
        )
    )

    outputs = {
        "output_dir": str(output_dir),
        "work_dir": str(work_dir),
        "output_vcf": str(output_vcf),
        "output_bam": str(merged_bam),
        "output_bam_index": str(merged_bam) + ".bai",
        "gtf_used": str(discovered.gtf),
        "vcf_used": str(discovered.vcf),
    }
    outputs.update(parent_outputs)
    return outputs
