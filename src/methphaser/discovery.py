from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, List, Optional, Sequence


class DiscoveryError(RuntimeError):
    """Raised when companion files cannot be inferred from BAM path."""


@dataclass(frozen=True)
class PhaseInputs:
    gtf: Path
    vcf: Path


def infer_run_token(bam_path: Path) -> Optional[str]:
    match = re.search(r"(PAW\d+)", bam_path.name)
    if match:
        return match.group(1)
    return None


def _unique_existing_dirs(paths: Iterable[Path]) -> List[Path]:
    seen = set()
    ordered: List[Path] = []
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_dir():
            ordered.append(resolved)
    return ordered


def candidate_output_dirs(bam_path: Path, run_token: Optional[str]) -> List[Path]:
    bam_dir = bam_path.parent
    search_dirs: List[Path] = [bam_dir]

    for parent in [bam_dir, bam_dir.parent, bam_dir.parent.parent]:
        search_dirs.append(parent)
        if run_token:
            search_dirs.append(parent / f"{run_token}_output")
            search_dirs.append(parent / f"{run_token}" / f"{run_token}_output")

    # Common ONT layout: .../giab_2025.01_mapped/*.bam and sibling ../giab_2025.01/<RUN>_output/
    if run_token:
        search_dirs.extend(
            [
                bam_dir.parent / "giab_2025.01" / f"{run_token}_output",
                bam_dir.parent.parent / "giab_2025.01" / f"{run_token}_output",
            ]
        )

    # Also inspect sibling output directories directly (fast, non-recursive).
    for parent in [bam_dir, bam_dir.parent, bam_dir.parent.parent]:
        search_dirs.extend(parent.glob("*_output"))

    return _unique_existing_dirs(search_dirs)


def _pick_first(directory: Path, patterns: Sequence[str]) -> Optional[Path]:
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def _discover_in_directory(directory: Path) -> tuple[Optional[Path], Optional[Path]]:
    gtf = _pick_first(
        directory,
        [
            "*.wf_snp.haploblocks.gtf",
            "*.haploblocks.gtf",
            "*.gtf",
        ],
    )
    vcf = _pick_first(
        directory,
        [
            "*.wf_snp.vcf.gz",
            "*.snp.vcf.gz",
            "*.vcf.gz",
        ],
    )
    return gtf, vcf


def discover_phase_inputs(
    bam_file: str,
    gtf: Optional[str] = None,
    vcf: Optional[str] = None,
) -> PhaseInputs:
    bam_path = Path(bam_file)
    if not bam_path.exists():
        raise DiscoveryError(f"BAM file does not exist: {bam_file}")

    resolved_gtf = Path(gtf).resolve() if gtf else None
    resolved_vcf = Path(vcf).resolve() if vcf else None

    if resolved_gtf and not resolved_gtf.exists():
        raise DiscoveryError(f"GTF file does not exist: {resolved_gtf}")
    if resolved_vcf and not resolved_vcf.exists():
        raise DiscoveryError(f"VCF file does not exist: {resolved_vcf}")

    if resolved_gtf and resolved_vcf:
        return PhaseInputs(gtf=resolved_gtf, vcf=resolved_vcf)

    run_token = infer_run_token(bam_path)
    search_dirs = candidate_output_dirs(bam_path, run_token)

    found_gtf = resolved_gtf
    found_vcf = resolved_vcf

    # Pass 1: prefer same directory containing both files.
    if not (found_gtf and found_vcf):
        for directory in search_dirs:
            gtf_candidate, vcf_candidate = _discover_in_directory(directory)
            if not found_gtf and gtf_candidate:
                found_gtf = gtf_candidate
            if not found_vcf and vcf_candidate:
                found_vcf = vcf_candidate
            if found_gtf and found_vcf:
                break

    if not found_gtf or not found_vcf:
        checked = "\n  - ".join(str(path) for path in search_dirs)
        raise DiscoveryError(
            "Unable to auto-discover haploblocks GTF and SNP VCF. "
            "Provide --gtf and --vcf explicitly.\n"
            f"Checked directories:\n  - {checked}"
        )

    return PhaseInputs(gtf=found_gtf.resolve(), vcf=found_vcf.resolve())
