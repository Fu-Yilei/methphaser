from __future__ import annotations

from dataclasses import dataclass
import math
import os
from multiprocessing import Pool
from typing import List, Sequence, Tuple

import pandas as pd

from . import step1_core

GTF_COLUMNS = [
    "chr",
    "phasing",
    "ex/intron",
    "start",
    "end",
    "1",
    "strand",
    "2",
    "info",
]
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


@dataclass(frozen=True)
class Step1Config:
    bam_file: str
    reference: str
    gtf: str
    output_dir: str
    threads: int = 1
    max_gap_len: int = -1
    minimum_vote: float = 0.65
    assignment_min: int = 2
    max_iterations: int = 10
    region_strategy: str = "full-span"
    boundary_window_bp: int = 1000
    chromosomes: Sequence[str] | None = None


def _resolve_target_chromosomes(
    phased_block_df: pd.DataFrame, chromosomes: Sequence[str] | None
) -> List[str]:
    available = list(dict.fromkeys(phased_block_df["chr"]))
    if chromosomes:
        requested = [chrom for chrom in chromosomes if chrom in available]
        return requested
    return [chrom for chrom in AUTOSOMES if chrom in available]


def compute_skipping_pairs(phased_chr_df: pd.DataFrame, max_gap_len: int) -> List[int]:
    skipping_pair_start_list = [-1]
    if len(phased_chr_df) <= 1:
        return skipping_pair_start_list

    distances = []
    rows = list(phased_chr_df.itertuples(index=False))
    for index in range(len(rows) - 1):
        distances.append((index, rows[index + 1].start - rows[index].end))
    if not distances:
        return skipping_pair_start_list

    if max_gap_len == -1:
        largest_gap_index = max(distances, key=lambda item: item[1])[0]
        skipping_pair_start_list.append(largest_gap_index)
        return skipping_pair_start_list

    if max_gap_len == -2:
        return skipping_pair_start_list

    for index, gap in distances:
        if gap >= max_gap_len:
            skipping_pair_start_list.append(index)
    return skipping_pair_start_list


def build_intervals(start_block: int, end_block: int, workers: int) -> List[Tuple[int, int]]:
    total = end_block - start_block
    if total <= 0:
        return []
    workers = max(1, min(workers, total))
    block_size = math.ceil(total / workers)
    intervals: List[Tuple[int, int]] = []
    for worker_index in range(workers):
        start_pos = start_block + (worker_index * block_size)
        end_pos = min(start_block + ((worker_index + 1) * block_size), end_block)
        if start_pos >= end_pos:
            continue
        intervals.append((start_pos, end_pos))
    return intervals


def _run_block_range(
    config: Step1Config,
    start_pos: int,
    end_pos: int,
    chromosome: str,
    skipping_pair: str,
) -> None:
    os.makedirs(os.path.join(config.output_dir, chromosome), exist_ok=True)
    read_assignment_dir = os.path.join(config.output_dir, f"{chromosome}_read_assignment")
    os.makedirs(read_assignment_dir, exist_ok=True)

    relationship_csv = os.path.join(config.output_dir, chromosome, f"{start_pos}_{end_pos}.csv")
    step1_core.main(
        [
            "-b",
            config.bam_file,
            "-r",
            config.reference,
            "-p",
            config.gtf,
            "-n",
            f"{start_pos},{end_pos}",
            "-m",
            chromosome,
            "-o",
            relationship_csv,
            "-k",
            str(config.max_iterations),
            "-s",
            skipping_pair,
            "--region_strategy",
            config.region_strategy,
            "--boundary_window_bp",
            str(config.boundary_window_bp),
            "-a",
            str(config.assignment_min),
            "-c",
            str(config.minimum_vote),
            "-ra",
            read_assignment_dir,
        ]
    )


def run_step1(config: Step1Config) -> List[str]:
    phased_block_df = pd.read_csv(
        config.gtf,
        header=None,
        sep="\t",
        names=GTF_COLUMNS,
    )

    target_chromosomes = _resolve_target_chromosomes(phased_block_df, config.chromosomes)
    if not target_chromosomes:
        raise ValueError("no target chromosomes found in GTF")

    processed: List[str] = []
    for chromosome in target_chromosomes:
        phased_chr_df = phased_block_df[phased_block_df["chr"] == chromosome].reset_index(
            drop=True
        )
        chr_block_num = len(phased_chr_df)
        if chr_block_num < 2:
            continue

        skipping_pair_start_list = compute_skipping_pairs(phased_chr_df, config.max_gap_len)
        intervals = build_intervals(0, chr_block_num, config.threads)
        if not intervals:
            continue

        worker_count = max(1, min(config.threads, len(intervals)))
        tasks = [
            (config, start_pos, end_pos, chromosome, str(skipping_pair_start_list))
            for start_pos, end_pos in intervals
        ]

        if worker_count == 1:
            for task in tasks:
                _run_block_range(*task)
        else:
            with Pool(worker_count) as pool:
                pool.starmap(_run_block_range, tasks)

        processed.append(chromosome)

    if not processed:
        raise ValueError("step1 did not process any chromosomes")
    return processed
