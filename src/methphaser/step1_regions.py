from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Iterable, List, Sequence, Tuple

Region = Tuple[int, int]


@dataclass(frozen=True)
class PhaseRegionPlan:
    span: Region
    windows: Tuple[Region, ...]


def _merge_windows(windows: Iterable[Region]) -> Tuple[Region, ...]:
    ordered = sorted((int(start), int(end)) for start, end in windows if int(start) < int(end))
    if not ordered:
        return tuple()

    merged: List[Region] = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return tuple(merged)


def parse_skipping_pairs(skipping_pair_str: str) -> tuple[List[int], List[int]]:
    cleaned = skipping_pair_str.strip().replace(" ", "")
    if not (cleaned.startswith("[") and cleaned.endswith("]")):
        raise ValueError(f"invalid skipping pair list: {skipping_pair_str}")

    inner = cleaned[1:-1]
    if not inner:
        return [], []

    starts = [int(item) for item in inner.split(",") if item]
    if starts == [-2]:
        return [], []

    ends = [item + 1 for item in starts]
    return starts, ends


def build_full_span_region_plans(
    phased_region_list: Sequence[Region],
    skipping_pair_start: Sequence[int],
) -> List[PhaseRegionPlan]:
    if not phased_region_list:
        return []
    if len(phased_region_list) == 1:
        region = phased_region_list[0]
        return [PhaseRegionPlan(span=region, windows=(region,))]

    skipping_pair_end = [item + 1 for item in skipping_pair_start]
    expanded: List[Region] = [(phased_region_list[0][0], phased_region_list[1][0])]

    for index, _ in enumerate(phased_region_list):
        if index == 0 or index == len(phased_region_list) - 1:
            continue
        if index in skipping_pair_start:
            left_most = phased_region_list[index - 1][1]
            right_most = phased_region_list[index][1]
        elif index in skipping_pair_end:
            left_most = phased_region_list[index][0]
            right_most = phased_region_list[index + 1][0]
        else:
            left_most = phased_region_list[index - 1][1]
            right_most = phased_region_list[index + 1][0]
        expanded.append((left_most, right_most))

    expanded.append((phased_region_list[-2][1], phased_region_list[-1][1]))
    return [PhaseRegionPlan(span=region, windows=(region,)) for region in expanded]


def build_boundary_window_region_plans(
    phased_region_list: Sequence[Region],
    skipping_pair_start: Sequence[int],
    boundary_window_bp: int,
) -> List[PhaseRegionPlan]:
    if boundary_window_bp <= 0:
        raise ValueError("boundary_window_bp must be > 0")

    skip_pairs = set(skipping_pair_start)
    plans: List[PhaseRegionPlan] = []
    total_blocks = len(phased_region_list)

    for block_index, (block_start, block_end) in enumerate(phased_region_list):
        windows: List[Region] = []

        include_left = block_index > 0 and (block_index - 1) not in skip_pairs
        include_right = block_index < (total_blocks - 1) and block_index not in skip_pairs

        if include_left:
            left_end = min(block_end, block_start + boundary_window_bp)
            if block_start < left_end:
                windows.append((block_start, left_end))

        if include_right:
            right_start = max(block_start, block_end - boundary_window_bp)
            if right_start < block_end:
                windows.append((right_start, block_end))

        merged = _merge_windows(windows)
        if not merged:
            merged = ((block_start, block_end),)

        span = (merged[0][0], merged[-1][1])
        plans.append(PhaseRegionPlan(span=span, windows=merged))

    return plans


def encode_windows_for_filename(windows: Sequence[Region]) -> str:
    return "_".join(f"{start}-{end}" for start, end in windows)


def parse_assignment_filename(file_name: str) -> tuple[int, Region, List[Region]]:
    stem = os.path.splitext(file_name)[0]
    parts = stem.split("_")
    if len(parts) < 3:
        raise ValueError(f"invalid assignment file name: {file_name}")

    block_index = int(parts[0])
    span = (int(parts[1]), int(parts[2]))
    windows: List[Region] = []

    for token in parts[3:]:
        if "-" not in token:
            continue
        start_s, end_s = token.split("-", 1)
        start = int(start_s)
        end = int(end_s)
        if start < end:
            windows.append((start, end))

    if not windows:
        windows = [span]

    return block_index, span, windows
