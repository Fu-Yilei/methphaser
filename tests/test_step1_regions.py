from methphaser.step1_regions import (
    build_boundary_window_region_plans,
    build_full_span_region_plans,
    encode_windows_for_filename,
    parse_assignment_filename,
)


def test_build_full_span_region_plans_matches_legacy_layout() -> None:
    phased_regions = [
        (5012610, 5014036),
        (5047142, 5053990),
        (5068630, 5154647),
        (5216322, 5390428),
        (6068657, 6074639),
        (6080717, 6132733),
    ]

    plans = build_full_span_region_plans(phased_regions, skipping_pair_start=[-1])

    assert [plan.span for plan in plans] == [
        (5012610, 5047142),
        (5014036, 5068630),
        (5053990, 5216322),
        (5154647, 6068657),
        (5390428, 6080717),
        (6074639, 6132733),
    ]
    assert all(plan.windows == (plan.span,) for plan in plans)


def test_build_boundary_window_region_plans_uses_block_edges_only() -> None:
    phased_regions = [
        (1000, 2000),
        (3000, 7000),
        (9000, 10000),
    ]

    plans = build_boundary_window_region_plans(
        phased_regions,
        skipping_pair_start=[],
        boundary_window_bp=500,
    )

    assert [plan.windows for plan in plans] == [
        ((1500, 2000),),
        ((3000, 3500), (6500, 7000)),
        ((9000, 9500),),
    ]


def test_assignment_filename_round_trip() -> None:
    windows = [(100, 200), (800, 900)]
    suffix = encode_windows_for_filename(windows)
    file_name = f"3_100_900_{suffix}.csv"

    index, span, parsed_windows = parse_assignment_filename(file_name)

    assert index == 3
    assert span == (100, 900)
    assert parsed_windows == windows
