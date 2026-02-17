#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Sequence
import urllib.request


HUMAN_ICR_BASE = "https://jb2.humanicr.org/data/json/humanicr"
TRACKS_ROOT = f"{HUMAN_ICR_BASE}/tracks"
TRACKS = [
    "known_ICR_25",
    "putative_ICRs_35_65_final_N1488",
    "ICRs_2_or_more_method_val_N1348",
    "ICRs_2_or_more_method_N4239",
]
OOCYTE_TRACK = "Oocyte_100"
SPERM_TRACK = "Sperm_with_sue_100"
PRIMARY_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
METHYLATION_HIGH_THRESHOLD = 0.75
METHYLATION_LOW_THRESHOLD = 0.25
METHYLATION_MIN_DIFF = 0.35

BED_HEADER = (
    "# chrom\tstart\tend\tregion_id\tsource_track\tsource_database\t"
    "oocyte_mean\tsperm_mean\tmethylated_parent\tdirection_confidence\t"
    "oocyte_cpg_sites\tsperm_cpg_sites\tpaternal_haplotype_rule"
)


def fetch_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def format_score(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.4f}"


def infer_methylated_parent(oocyte_mean: float | None, sperm_mean: float | None) -> str:
    if oocyte_mean is None or sperm_mean is None:
        return "unknown"
    diff = abs(oocyte_mean - sperm_mean)
    if diff < METHYLATION_MIN_DIFF:
        return "unknown"
    if (
        oocyte_mean >= METHYLATION_HIGH_THRESHOLD
        and sperm_mean <= METHYLATION_LOW_THRESHOLD
    ):
        return "maternal"
    if (
        sperm_mean >= METHYLATION_HIGH_THRESHOLD
        and oocyte_mean <= METHYLATION_LOW_THRESHOLD
    ):
        return "paternal"
    return "unknown"


def paternal_rule_for_parent(methylated_parent: str) -> str:
    if methylated_parent == "maternal":
        return "low_methylation_haplotype"
    if methylated_parent == "paternal":
        return "high_methylation_haplotype"
    return "unknown"


class NCListScoreReader:
    def __init__(self, track_name: str):
        self.track_name = track_name
        self.track_meta_cache: Dict[str, dict] = {}
        self.chunk_cache: Dict[tuple[str, int], list] = {}

    def _track_url(self, chrom: str) -> str:
        return f"{TRACKS_ROOT}/{self.track_name}/{chrom}/trackData.json"

    def _chunk_url(self, chrom: str, chunk_id: int) -> str:
        return f"{TRACKS_ROOT}/{self.track_name}/{chrom}/lf-{chunk_id}.json"

    def _load_track_meta(self, chrom: str) -> dict | None:
        if chrom in self.track_meta_cache:
            return self.track_meta_cache[chrom]

        try:
            track_data = fetch_json(self._track_url(chrom))
        except Exception:
            self.track_meta_cache[chrom] = None
            return None

        intervals = track_data.get("intervals", {})
        classes = intervals.get("classes", [])
        lazy_class = int(intervals.get("lazyClass", 1))
        root_entries = intervals.get("nclist", [])

        class_attr_index: Dict[int, Dict[str, int]] = {}
        for class_id, class_def in enumerate(classes):
            attributes = [str(item) for item in class_def.get("attributes", [])]
            class_attr_index[class_id] = {name: idx for idx, name in enumerate(attributes)}

        meta = {
            "root_entries": root_entries,
            "lazy_class": lazy_class,
            "class_attr_index": class_attr_index,
        }
        self.track_meta_cache[chrom] = meta
        return meta

    def _get_chunk_entries(self, chrom: str, chunk_id: int) -> list:
        cache_key = (chrom, int(chunk_id))
        if cache_key in self.chunk_cache:
            return self.chunk_cache[cache_key]
        entries = fetch_json(self._chunk_url(chrom, int(chunk_id)))
        self.chunk_cache[cache_key] = entries
        return entries

    def _entry_span(self, entry: list, attr_index: Dict[str, int]) -> tuple[int, int]:
        # NCList arrays always store class ID at index 0, attributes start at index 1.
        start_i = attr_index["Start"] + 1
        end_i = attr_index["End"] + 1
        return int(entry[start_i]), int(entry[end_i])

    def _overlap_region_indexes(
        self,
        start: int,
        end: int,
        regions: Sequence[tuple[int, int]],
    ) -> list[int]:
        overlap_indexes: list[int] = []
        for idx, (region_start, region_end) in enumerate(regions):
            if region_end <= start:
                continue
            if region_start >= end:
                break
            overlap_indexes.append(idx)
        return overlap_indexes

    def _aggregate_scores_from_entries(
        self,
        chrom: str,
        entries: list,
        meta: dict,
        regions: Sequence[tuple[int, int]],
        sums: list[float],
        counts: list[int],
    ) -> None:
        class_attr_index: Dict[int, Dict[str, int]] = meta["class_attr_index"]
        lazy_class = int(meta["lazy_class"])

        for entry in entries:
            if not entry:
                continue
            class_id = int(entry[0])
            attr_index = class_attr_index.get(class_id, {})
            if "Start" not in attr_index or "End" not in attr_index:
                continue
            start, end = self._entry_span(entry, attr_index)
            overlap_indexes = self._overlap_region_indexes(start, end, regions)
            if not overlap_indexes:
                continue

            if class_id == lazy_class:
                if "Chunk" not in attr_index:
                    continue
                chunk_i = attr_index["Chunk"] + 1
                child_chunk = int(entry[chunk_i])
                child_entries = self._get_chunk_entries(chrom, child_chunk)
                self._aggregate_scores_from_entries(
                    chrom,
                    child_entries,
                    meta,
                    regions,
                    sums,
                    counts,
                )
                continue

            if "Score" not in attr_index:
                continue
            score_i = attr_index["Score"] + 1
            try:
                score = float(entry[score_i])
            except Exception:
                continue
            for idx in overlap_indexes:
                sums[idx] += score
                counts[idx] += 1

    def aggregate_scores_for_regions(
        self,
        chrom: str,
        regions: Sequence[tuple[int, int]],
    ) -> tuple[list[float], list[int]]:
        if not regions:
            return [], []

        sums = [0.0] * len(regions)
        counts = [0] * len(regions)
        meta = self._load_track_meta(chrom)
        if meta is None:
            return sums, counts
        self._aggregate_scores_from_entries(
            chrom,
            meta["root_entries"],
            meta,
            regions,
            sums,
            counts,
        )

        return sums, counts


def _iter_track_rows(track: str, chromosomes: Sequence[str]) -> Iterable[tuple[str, int, int, str]]:
    for chrom in chromosomes:
        url = f"{TRACKS_ROOT}/{track}/{chrom}/trackData.json"
        try:
            data = fetch_json(url)
        except Exception:
            continue

        for entry in data.get("intervals", {}).get("nclist", []):
            if not entry or entry[0] != 0:
                continue
            start = int(entry[1])
            end = int(entry[2])
            if start >= end:
                continue
            name = str(entry[4]) if len(entry) > 4 else f"{track}:{chrom}:{start}-{end}"
            yield chrom, start, end, name


def build_database(
    output_path: Path,
    tracks: Sequence[str],
    chromosomes: Sequence[str],
) -> int:
    sperm_reader = NCListScoreReader(SPERM_TRACK)
    oocyte_reader = NCListScoreReader(OOCYTE_TRACK)

    region_rows = []
    for track in tracks:
        for chrom, start, end, name in _iter_track_rows(track, chromosomes):
            region_rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "region_id": name,
                    "source_track": track,
                    "source_database": "HumanICR",
                }
            )

    rows = []
    for chrom in chromosomes:
        chrom_rows = [row for row in region_rows if row["chrom"] == chrom]
        if not chrom_rows:
            continue
        print(f"processing {chrom} ({len(chrom_rows)} regions)")
        chrom_rows.sort(key=lambda item: (item["start"], item["end"], item["source_track"], item["region_id"]))
        spans = [(int(item["start"]), int(item["end"])) for item in chrom_rows]
        sperm_sums, sperm_counts = sperm_reader.aggregate_scores_for_regions(chrom, spans)
        oocyte_sums, oocyte_counts = oocyte_reader.aggregate_scores_for_regions(chrom, spans)

        for idx, row in enumerate(chrom_rows):
            sperm_mean = (sperm_sums[idx] / sperm_counts[idx]) if sperm_counts[idx] > 0 else None
            oocyte_mean = (oocyte_sums[idx] / oocyte_counts[idx]) if oocyte_counts[idx] > 0 else None
            methylated_parent = infer_methylated_parent(oocyte_mean, sperm_mean)
            direction_confidence = (
                abs(float(oocyte_mean) - float(sperm_mean))
                if oocyte_mean is not None and sperm_mean is not None
                else None
            )
            rows.append(
                (
                    row["chrom"],
                    row["start"],
                    row["end"],
                    row["region_id"],
                    row["source_track"],
                    row["source_database"],
                    format_score(oocyte_mean),
                    format_score(sperm_mean),
                    methylated_parent,
                    format_score(direction_confidence),
                    oocyte_counts[idx],
                    sperm_counts[idx],
                    paternal_rule_for_parent(methylated_parent),
                )
            )

    rows.sort(key=lambda row: (row[0], row[1], row[2], row[4], row[3], row[8]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(BED_HEADER)
        handle.write("\n")
        for row in rows:
            handle.write("\t".join(map(str, row)))
            handle.write("\n")
    return len(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build merged imprinting BED from HumanICR track JSON files",
    )
    parser.add_argument(
        "--output",
        default="src/methphaser/data/imprinting_regions.hg38.bed",
        help="output BED path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = Path(args.output).resolve()
    count = build_database(output, TRACKS, PRIMARY_CHROMS)
    print(f"wrote {count} regions to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
