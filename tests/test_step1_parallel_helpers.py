import pandas as pd

from methphaser.step1_parallel import (
    _resolve_target_chromosomes,
    build_intervals,
    compute_skipping_pairs,
)


def test_resolve_target_chromosomes_defaults_to_autosomes() -> None:
    df = pd.DataFrame({"chr": ["chrX", "chr1", "chr2", "chr21", "chrM"]})
    assert _resolve_target_chromosomes(df, None) == ["chr1", "chr2", "chr21"]


def test_build_intervals_evenly() -> None:
    assert build_intervals(0, 10, 3) == [(0, 4), (4, 8), (8, 10)]


def test_compute_skipping_pairs_largest_gap() -> None:
    df = pd.DataFrame(
        [
            {"start": 100, "end": 200},
            {"start": 210, "end": 300},
            {"start": 1000, "end": 1100},
        ]
    )
    # Largest gap is between block 1 and 2 => index 1
    assert compute_skipping_pairs(df, -1) == [-1, 1]
