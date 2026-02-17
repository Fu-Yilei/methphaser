from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import pysam
from pysam import VariantFile


METHYLATION_IDENTIFIER_0 = ("C", 0, "m")
METHYLATION_IDENTIFIER_1 = ("C", 1, "m")

IMPRINTING_COLUMNS = [
    "chrom",
    "start",
    "end",
    "region_id",
    "source_track",
    "source_database",
    "oocyte_mean",
    "sperm_mean",
    "methylated_parent",
    "direction_confidence",
    "oocyte_cpg_sites",
    "sperm_cpg_sites",
    "paternal_haplotype_rule",
]

DIRECTIONAL_IMPRINTING_COLUMNS = [
    "chrom",
    "start",
    "end",
    "region_id",
    "source_track",
    "source_database",
    "methylated_parent",
    "direction_confidence",
]

PHASEBLOCK_COLUMNS = [
    "chrom",
    "phase_set",
    "start",
    "end",
    "variant_count",
]

ASSIGNMENT_COLUMNS = [
    "chrom",
    "phase_set",
    "phaseblock_start",
    "phaseblock_end",
    "variant_count",
    "parent_assignment_status",
    "paternal_haplotype",
    "h1_parent",
    "h2_parent",
    "assignment_basis",
    "directional_region_count",
    "decisive_region_count",
    "h1_vote_weight",
    "h2_vote_weight",
    "vote_margin",
    "overlap_bp",
    "overlap_source_tracks",
    "overlap_source_databases",
    "overlap_region_ids",
]


DEFAULT_IMPRINTING_BED = (
    Path(__file__).resolve().parent / "data" / "imprinting_regions.hg38.bed"
).resolve()


@dataclass(frozen=True)
class ParentAssignmentConfig:
    output_vcf: str
    phaseblock_assignment_output: str
    input_bam: str | None = None
    imprinting_bed: str | None = None
    min_sites_per_haplotype: int = 3
    min_methylation_diff: float = 0.10


def default_imprinting_bed_path() -> Path:
    return DEFAULT_IMPRINTING_BED


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _normalize_parent(value: object) -> str:
    if value is None:
        return "unknown"
    text = str(value).strip().lower()
    if text in {"maternal", "mother", "mat"}:
        return "maternal"
    if text in {"paternal", "father", "pat"}:
        return "paternal"
    return "unknown"


def load_imprinting_regions(imprinting_bed: str | None = None) -> pd.DataFrame:
    bed_path = Path(imprinting_bed).resolve() if imprinting_bed else DEFAULT_IMPRINTING_BED
    if not bed_path.exists():
        raise FileNotFoundError(f"imprinting BED database does not exist: {bed_path}")

    regions = pd.read_csv(
        bed_path,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
    )
    if regions.empty:
        return pd.DataFrame(columns=IMPRINTING_COLUMNS)

    if regions.shape[1] < 3:
        raise ValueError(f"imprinting BED is missing required columns: {bed_path}")

    # Support older 6-column files by filling directional metadata as unknown.
    while regions.shape[1] < len(IMPRINTING_COLUMNS):
        regions[regions.shape[1]] = ""

    regions = regions.iloc[:, : len(IMPRINTING_COLUMNS)].copy()
    regions.columns = IMPRINTING_COLUMNS
    regions["start"] = pd.to_numeric(regions["start"], errors="coerce")
    regions["end"] = pd.to_numeric(regions["end"], errors="coerce")
    regions = regions.dropna(subset=["start", "end"]).copy()
    regions["start"] = regions["start"].astype(int)
    regions["end"] = regions["end"].astype(int)
    regions = regions[regions["start"] < regions["end"]].copy()
    regions["chrom"] = regions["chrom"].astype(str)
    regions = regions[regions["chrom"] != ""].copy()
    regions["region_id"] = regions["region_id"].fillna("").astype(str)
    regions["source_track"] = regions["source_track"].fillna("").astype(str)
    regions["source_database"] = regions["source_database"].fillna("").astype(str)
    regions["methylated_parent"] = regions["methylated_parent"].apply(_normalize_parent)
    regions["direction_confidence"] = regions["direction_confidence"].apply(_float_or_none)
    regions["oocyte_mean"] = regions["oocyte_mean"].apply(_float_or_none)
    regions["sperm_mean"] = regions["sperm_mean"].apply(_float_or_none)

    return regions.sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def extract_directional_regions(regions: pd.DataFrame) -> pd.DataFrame:
    if regions.empty:
        return pd.DataFrame(columns=DIRECTIONAL_IMPRINTING_COLUMNS)
    directional = regions[regions["methylated_parent"].isin(["maternal", "paternal"])].copy()
    return directional[DIRECTIONAL_IMPRINTING_COLUMNS].sort_values(
        ["chrom", "start", "end"]
    ).reset_index(drop=True)


def extract_phaseblocks_from_vcf(vcf_path: str) -> pd.DataFrame:
    phaseblock_dict: Dict[Tuple[str, int], dict] = {}
    vcf_file = VariantFile(vcf_path)
    samples = list(vcf_file.header.samples)
    if not samples:
        return pd.DataFrame(columns=PHASEBLOCK_COLUMNS)
    sample_name = samples[0]

    for rec in vcf_file:
        sample = rec.samples[sample_name]
        gt = sample.get("GT")
        ps = sample.get("PS")
        if gt is None or len(gt) < 2 or gt[0] is None or gt[1] is None:
            continue
        if gt[0] == gt[1]:
            continue
        if not sample.phased:
            continue
        if ps is None:
            continue

        key = (str(rec.contig), int(ps))
        rec_start = int(rec.pos)
        rec_end = int(rec.stop if rec.stop else rec.pos)
        if key not in phaseblock_dict:
            phaseblock_dict[key] = {
                "chrom": str(rec.contig),
                "phase_set": int(ps),
                "start": rec_start,
                "end": rec_end,
                "variant_count": 1,
            }
            continue

        phaseblock_dict[key]["start"] = min(phaseblock_dict[key]["start"], rec_start)
        phaseblock_dict[key]["end"] = max(phaseblock_dict[key]["end"], rec_end)
        phaseblock_dict[key]["variant_count"] += 1

    if not phaseblock_dict:
        return pd.DataFrame(columns=PHASEBLOCK_COLUMNS)
    phaseblocks = pd.DataFrame(phaseblock_dict.values(), columns=PHASEBLOCK_COLUMNS)
    return phaseblocks.sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def _candidate_chromosomes(chrom: str) -> List[str]:
    if chrom.startswith("chr"):
        return [chrom, chrom[3:]]
    return [chrom, f"chr{chrom}"]


def _directional_region_index(directional_regions: pd.DataFrame) -> Dict[str, List[dict]]:
    if directional_regions.empty:
        return {}
    index: Dict[str, List[dict]] = {}
    for chrom, group_df in directional_regions.groupby("chrom", sort=False):
        index[str(chrom)] = (
            group_df.sort_values(["start", "end"])[DIRECTIONAL_IMPRINTING_COLUMNS]
            .to_dict("records")
        )
    return index


def _pick_methylation_identifier(mm: dict) -> tuple[str, int, str] | None:
    if METHYLATION_IDENTIFIER_0 in mm:
        return METHYLATION_IDENTIFIER_0
    if METHYLATION_IDENTIFIER_1 in mm:
        return METHYLATION_IDENTIFIER_1
    return None


def _collect_haplotype_methylation_in_region(
    bam: pysam.AlignmentFile,
    chrom: str,
    start0: int,
    end0: int,
) -> dict:
    sums = {1: 0.0, 2: 0.0}
    counts = {1: 0, 2: 0}

    for reads in bam.fetch(chrom, int(start0), int(end0)):
        if not reads.has_tag("HP"):
            continue
        hp = int(reads.get_tag("HP"))
        if hp not in {1, 2}:
            continue
        mm = reads.modified_bases
        if not mm:
            continue
        methylation_identifier = _pick_methylation_identifier(mm)
        if not methylation_identifier:
            continue

        read_base_ref_loc = reads.get_reference_positions(full_length=True)
        for read_loc, score in mm[methylation_identifier]:
            if read_loc >= len(read_base_ref_loc):
                continue
            ref_pos = read_base_ref_loc[read_loc]
            if ref_pos is None:
                continue
            if int(start0) <= int(ref_pos) < int(end0):
                sums[hp] += float(score) / 255.0
                counts[hp] += 1

    hp1_mean = (sums[1] / counts[1]) if counts[1] > 0 else None
    hp2_mean = (sums[2] / counts[2]) if counts[2] > 0 else None
    return {
        "hp1_mean": hp1_mean,
        "hp2_mean": hp2_mean,
        "hp1_sites": counts[1],
        "hp2_sites": counts[2],
    }


def _infer_paternal_haplotype_vote(
    methylated_parent: str,
    hp1_mean: float | None,
    hp2_mean: float | None,
    hp1_sites: int,
    hp2_sites: int,
    min_sites_per_haplotype: int,
    min_methylation_diff: float,
) -> tuple[int | None, float]:
    if hp1_mean is None or hp2_mean is None:
        return None, 0.0
    if hp1_sites < min_sites_per_haplotype or hp2_sites < min_sites_per_haplotype:
        return None, 0.0

    meth_diff = float(hp1_mean) - float(hp2_mean)
    if abs(meth_diff) < float(min_methylation_diff):
        return None, 0.0

    if methylated_parent == "maternal":
        # Maternal allele should be more methylated.
        paternal_hap = 2 if meth_diff > 0 else 1
    elif methylated_parent == "paternal":
        # Paternal allele should be more methylated.
        paternal_hap = 1 if meth_diff > 0 else 2
    else:
        return None, 0.0

    return paternal_hap, abs(meth_diff)


def _collect_overlapping_directional_regions(
    region_index: Dict[str, List[dict]],
    chrom: str,
    block_start1: int,
    block_end1: int,
) -> List[dict]:
    overlaps: List[dict] = []
    block_start0 = int(block_start1) - 1
    block_end0 = int(block_end1)

    for candidate in _candidate_chromosomes(chrom):
        regions = region_index.get(candidate, [])
        for row in regions:
            region_start = int(row["start"])
            region_end = int(row["end"])
            if region_end <= block_start0:
                continue
            if region_start >= block_end0:
                break
            overlaps.append(row)
        if overlaps:
            return overlaps
    return overlaps


def annotate_phaseblocks(
    phaseblocks: pd.DataFrame,
    directional_regions: pd.DataFrame,
    input_bam: str | None,
    min_sites_per_haplotype: int,
    min_methylation_diff: float,
) -> pd.DataFrame:
    if phaseblocks.empty:
        return pd.DataFrame(columns=ASSIGNMENT_COLUMNS)

    region_index = _directional_region_index(directional_regions)
    bam = pysam.AlignmentFile(input_bam, "rb") if input_bam else None
    rows = []

    try:
        for block in phaseblocks.itertuples(index=False):
            overlaps = _collect_overlapping_directional_regions(
                region_index,
                str(block.chrom),
                int(block.start),
                int(block.end),
            )
            overlap_bp = 0
            source_tracks = set()
            source_databases = set()
            overlap_region_ids = set()
            h1_votes = 0.0
            h2_votes = 0.0
            decisive_region_count = 0

            for overlap in overlaps:
                block_start0 = int(block.start) - 1
                block_end0 = int(block.end)
                region_start = int(overlap["start"])
                region_end = int(overlap["end"])
                overlap_bp += max(
                    0,
                    min(block_end0, region_end) - max(block_start0, region_start),
                )
                source_tracks.add(str(overlap["source_track"]))
                source_databases.add(str(overlap["source_database"]))
                overlap_region_ids.add(str(overlap["region_id"]))

                if bam is None:
                    continue
                meth = _collect_haplotype_methylation_in_region(
                    bam,
                    str(block.chrom),
                    region_start,
                    region_end,
                )
                paternal_hap, vote_weight = _infer_paternal_haplotype_vote(
                    methylated_parent=str(overlap["methylated_parent"]),
                    hp1_mean=meth["hp1_mean"],
                    hp2_mean=meth["hp2_mean"],
                    hp1_sites=int(meth["hp1_sites"]),
                    hp2_sites=int(meth["hp2_sites"]),
                    min_sites_per_haplotype=max(1, int(min_sites_per_haplotype)),
                    min_methylation_diff=max(0.0, float(min_methylation_diff)),
                )
                if paternal_hap is None:
                    continue

                decisive_region_count += 1
                if paternal_hap == 1:
                    h1_votes += vote_weight
                elif paternal_hap == 2:
                    h2_votes += vote_weight

            if not overlaps:
                status = "random"
                paternal_haplotype = "random"
                h1_parent = "random"
                h2_parent = "random"
                basis = "no_directional_imprinting_region_overlap"
            elif bam is None:
                status = "random"
                paternal_haplotype = "random"
                h1_parent = "random"
                h2_parent = "random"
                basis = "directional_regions_found_but_no_bam_evidence"
            elif decisive_region_count == 0:
                status = "random"
                paternal_haplotype = "random"
                h1_parent = "random"
                h2_parent = "random"
                basis = "insufficient_haplotype_methylation_evidence"
            elif h1_votes == h2_votes:
                status = "random"
                paternal_haplotype = "random"
                h1_parent = "random"
                h2_parent = "random"
                basis = "conflicting_directional_votes"
            elif h1_votes > h2_votes:
                status = "assigned"
                paternal_haplotype = "H1"
                h1_parent = "father"
                h2_parent = "mother"
                basis = "directional_imprinting_plus_haplotype_methylation"
            else:
                status = "assigned"
                paternal_haplotype = "H2"
                h1_parent = "mother"
                h2_parent = "father"
                basis = "directional_imprinting_plus_haplotype_methylation"

            rows.append(
                {
                    "chrom": str(block.chrom),
                    "phase_set": int(block.phase_set),
                    "phaseblock_start": int(block.start),
                    "phaseblock_end": int(block.end),
                    "variant_count": int(block.variant_count),
                    "parent_assignment_status": status,
                    "paternal_haplotype": paternal_haplotype,
                    "h1_parent": h1_parent,
                    "h2_parent": h2_parent,
                    "assignment_basis": basis,
                    "directional_region_count": len(overlaps),
                    "decisive_region_count": decisive_region_count,
                    "h1_vote_weight": round(h1_votes, 4),
                    "h2_vote_weight": round(h2_votes, 4),
                    "vote_margin": round(abs(h1_votes - h2_votes), 4),
                    "overlap_bp": int(overlap_bp),
                    "overlap_source_tracks": ";".join(sorted(source_tracks)),
                    "overlap_source_databases": ";".join(sorted(source_databases)),
                    "overlap_region_ids": ";".join(sorted(overlap_region_ids)),
                }
            )
    finally:
        if bam is not None:
            bam.close()

    return pd.DataFrame(rows, columns=ASSIGNMENT_COLUMNS)


def run_parent_assignment(config: ParentAssignmentConfig) -> Dict[str, str]:
    regions = load_imprinting_regions(config.imprinting_bed)
    directional_regions = extract_directional_regions(regions)
    phaseblocks = extract_phaseblocks_from_vcf(config.output_vcf)
    assignments = annotate_phaseblocks(
        phaseblocks=phaseblocks,
        directional_regions=directional_regions,
        input_bam=config.input_bam,
        min_sites_per_haplotype=max(1, int(config.min_sites_per_haplotype)),
        min_methylation_diff=max(0.0, float(config.min_methylation_diff)),
    )

    output_path = Path(config.phaseblock_assignment_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(output_path, sep="\t", index=False)

    bed_path = Path(config.imprinting_bed).resolve() if config.imprinting_bed else DEFAULT_IMPRINTING_BED
    return {
        "imprinting_database_used": str(bed_path),
        "directional_imprinting_region_count": str(len(directional_regions)),
        "phaseblock_parent_assignment": str(output_path.resolve()),
    }
