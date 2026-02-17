from __future__ import annotations

import ast
import os
from dataclasses import dataclass
from itertools import groupby
from multiprocessing import Pool
from operator import itemgetter
from typing import Dict, List, Sequence, Tuple

import pandas as pd
import pysam
from pysam import VariantFile
from tqdm.auto import tqdm

from .step1_regions import parse_assignment_filename


def _parse_block(block_repr: str) -> Tuple[int, int]:
    start, end = ast.literal_eval(block_repr)
    return int(start), int(end)


def _iter_reads_from_regions(
    bam_file: pysam.AlignmentFile,
    chrom: str,
    regions: Sequence[Tuple[int, int]],
):
    seen_alignment_keys = set()
    for start, end in regions:
        for reads in bam_file.fetch(chrom, int(start), int(end), multiple_iterators=True):
            key = _alignment_key(reads)
            if key in seen_alignment_keys:
                continue
            seen_alignment_keys.add(key)
            yield reads


def _alignment_key(reads: pysam.AlignedSegment) -> Tuple[str, int, int, bool, bool]:
    return (
        reads.query_name,
        reads.reference_start,
        reads.reference_end,
        reads.is_secondary,
        reads.is_supplementary,
    )


@dataclass(frozen=True)
class Step2Config:
    input_bam_file: str
    meth_phasing_input_folder: str
    output_vcf: str
    output_bam_prefix: str
    vcf_called: str
    threads: int = 1
    minimum_coverage: int = 0
    voting_difference: float = 0.5
    include_unprocessed_contigs: bool = True
    chromosomes: Sequence[str] | None = None


def apply_flip_counter(relationship_to_block: str, flip_counter: int) -> int:
    if relationship_to_block == "same":
        return flip_counter
    return -1 * flip_counter


def result_filtering(comparison_df: pd.DataFrame, min_reads_num: int, min_variance: float) -> pd.DataFrame:
    low_support = (
        comparison_df["same_hap_num"] + comparison_df["diff_hap_num"] <= min_reads_num
    )
    comparison_df[low_support] = comparison_df[low_support].assign(
        myth_phasing_relationship="cannot decide"
    )

    vote_diff = (
        abs(comparison_df["same_hap_num"] - comparison_df["diff_hap_num"])
        / (comparison_df["same_hap_num"] + comparison_df["diff_hap_num"])
    )
    low_variance = vote_diff <= min_variance
    comparison_df[low_variance] = comparison_df[low_variance].assign(
        myth_phasing_relationship="cannot decide"
    )
    return comparison_df


def get_altered_block_start_loc(final_blocks: Sequence[Tuple[int, int]], original_block_start_loc: int) -> str:
    for start, end in final_blocks:
        if start <= original_block_start_loc <= end:
            return str(start)
    return "-1"


def get_all_final_blocks_dict(block_relationship_df: Dict[str, pd.DataFrame]) -> Dict[str, List[Tuple[int, int]]]:
    final_block_dict: Dict[str, List[Tuple[int, int]]] = {}
    for chrom in block_relationship_df:
        current_block_relationship_df = block_relationship_df[chrom]
        current_block_relationship_df_all = current_block_relationship_df[
            current_block_relationship_df.myth_phasing_relationship != "cannot decide"
        ]

        data = list(current_block_relationship_df_all.index)
        current_chr_index = []
        for _, group in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
            current_chr_index.append(list(map(itemgetter(1), group)))

        current_block_list: List[Tuple[int, int]] = []
        for group_index in current_chr_index:
            start = _parse_block(
                current_block_relationship_df_all.loc[group_index[0]].snp_phased_block_1
            )[0]
            end = _parse_block(
                current_block_relationship_df_all.loc[group_index[-1]].snp_phased_block_2
            )[1]
            current_block_list.append((start, end))

        for row_index in list(current_block_relationship_df.index):
            if row_index in data:
                continue
            if row_index == 0:
                start, end = _parse_block(
                    current_block_relationship_df.loc[row_index].snp_phased_block_1
                )
                current_block_list.append((start, end))
            start, end = _parse_block(
                current_block_relationship_df.loc[row_index].snp_phased_block_2
            )
            current_block_list.append((start, end))

        final_block_dict[chrom] = current_block_list

    return final_block_dict


def get_altered_vcf(
    original_vcf: str,
    output_vcf: str,
    block_relationship_dfs: Dict[str, pd.DataFrame],
    include_unprocessed_contigs: bool = True,
) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, List[str]], Dict[str, List[Tuple[str, int]]]]:
    called_vcf_file = VariantFile(original_vcf)
    final_block_dict: Dict[str, List[Tuple[int, int]]] = {}
    flipping_dict: Dict[str, List[Tuple[str, int]]] = {}
    remaining_dict: Dict[str, List[str]] = {}

    output_dir = os.path.dirname(output_vcf)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    updated_chroms = set()
    all_final_blocks = get_all_final_blocks_dict(block_relationship_dfs)

    with open(output_vcf, "w", encoding="utf-8") as altered_vcf_file:
        altered_vcf_file.writelines(str(called_vcf_file.header))
        for chrom in tqdm(block_relationship_dfs.keys()):
            updated_chroms.add(chrom)
            current_chrom_final_block = all_final_blocks[chrom]
            current_chrom_block_relationship_df = block_relationship_dfs[chrom]

            final_block_list: List[Tuple[int, int]] = []
            flipping_list: List[Tuple[str, int]] = []
            block_num = 0
            while block_num < len(current_chrom_block_relationship_df):
                current_block_row = current_chrom_block_relationship_df.iloc[block_num]
                if current_block_row.myth_phasing_relationship == "cannot decide":
                    block_num += 1
                    continue

                current_block_start = current_block_row.snp_phased_block_1
                block_start = _parse_block(current_block_start)[0]
                flipping_list.append((current_block_start, 1))
                flip_counter = apply_flip_counter(
                    current_block_row.myth_phasing_relationship,
                    1,
                )
                current_block_num = block_num
                block_num += 1

                while (
                    block_num < len(current_chrom_block_relationship_df)
                    and current_chrom_block_relationship_df.iloc[
                        block_num
                    ].myth_phasing_relationship
                    != "cannot decide"
                ):
                    flipping_list.append(
                        (
                            current_chrom_block_relationship_df.iloc[
                                block_num
                            ].snp_phased_block_1,
                            flip_counter,
                        )
                    )
                    flip_counter = apply_flip_counter(
                        current_chrom_block_relationship_df.iloc[
                            block_num
                        ].myth_phasing_relationship,
                        flip_counter,
                    )
                    block_num += 1

                flipping_list.append(
                    (
                        current_chrom_block_relationship_df.iloc[
                            block_num - 1
                        ].snp_phased_block_2,
                        flip_counter,
                    )
                )
                block_end = _parse_block(
                    current_chrom_block_relationship_df.iloc[
                        block_num - 1
                    ].snp_phased_block_2
                )[1]
                final_block_list.append((block_start, block_end))

            final_block_dict[chrom] = final_block_list
            flipping_dict[chrom] = flipping_list

            unphased_list = []
            for index, block_repr in enumerate(
                current_chrom_block_relationship_df.snp_phased_block_1[:-1]
            ):
                left_end = _parse_block(block_repr)[1]
                right_start = _parse_block(
                    current_chrom_block_relationship_df.snp_phased_block_1[index + 1]
                )[0]
                unphased_list.append((left_end, right_start))

            for start, end in unphased_list:
                for rec in called_vcf_file.fetch(chrom, start, end - 1):
                    altered_vcf_file.writelines(str(rec))

            for block_repr, flip_flag in flipping_list:
                block_start, block_end = _parse_block(block_repr)
                for rec_obj in called_vcf_file.fetch(chrom, block_start - 1, block_end):
                    rec = str(rec_obj)
                    if "1|0" in rec or "0|1" in rec:
                        if "PS" in rec.split("\t")[-2]:
                            split_rec = rec.split("\t")
                            ps_tag_location = split_rec[-2].split(":").index("PS")
                            start_loc = split_rec[-1].split(":")[ps_tag_location]
                            split_rec[-1] = split_rec[-1].replace(
                                start_loc,
                                get_altered_block_start_loc(
                                    current_chrom_final_block,
                                    int(start_loc),
                                ),
                            )
                            rec = "\t".join(split_rec)

                        if flip_flag == -1:
                            rec = rec.replace("1|0", "tmp").replace("0|1", "1|0").replace("tmp", "0|1")
                    altered_vcf_file.writelines(rec)

            flipped_blocks = {block_repr for block_repr, _ in flipping_list}
            remaining_list = [
                block_repr
                for block_repr in list(current_chrom_block_relationship_df.snp_phased_block_1)
                if block_repr not in flipped_blocks
            ]
            last_block = current_chrom_block_relationship_df.iloc[-1].snp_phased_block_2
            if last_block not in flipped_blocks:
                remaining_list.append(last_block)
            remaining_dict[chrom] = remaining_list

            for block_repr in remaining_list:
                block_start, block_end = _parse_block(block_repr)
                for rec in called_vcf_file.fetch(chrom, block_start - 1, block_end):
                    altered_vcf_file.writelines(str(rec))

        if include_unprocessed_contigs:
            for chrom in called_vcf_file.header.contigs:
                if chrom in updated_chroms:
                    continue
                for record in called_vcf_file.fetch(chrom):
                    altered_vcf_file.write(str(record))

    return final_block_dict, remaining_dict, flipping_dict


def get_altered_bam(
    flipping_dict: Dict[str, List[Tuple[str, int]]],
    remaining_dict: Dict[str, List[str]],
    block_relationship_dfs: Dict[str, pd.DataFrame],
    input_bam_file: str,
    original_bam_file: str,
    modified_bam_file: str,
    assignment_path: str,
    chrom: str,
) -> None:
    input_bam = pysam.AlignmentFile(input_bam_file, "rb")
    original_bam = pysam.AlignmentFile(original_bam_file, "rb")
    modified_bam = pysam.AlignmentFile(modified_bam_file, "wb", template=original_bam)
    written_alignment_keys = set()

    def write_read(reads: pysam.AlignedSegment) -> None:
        key = _alignment_key(reads)
        if key in written_alignment_keys:
            return
        modified_bam.write(reads)
        written_alignment_keys.add(key)

    sorted_read_assignment_files = sorted(
        os.listdir(assignment_path),
        key=lambda file_name: parse_assignment_filename(file_name)[0],
    )

    flipping_dict_chrom = dict(flipping_dict[chrom])
    snp_block_flipping_chrom = dict(zip(remaining_dict[chrom], [1] * len(remaining_dict[chrom])))
    snp_block_flipping_chrom.update(flipping_dict_chrom)

    current_all_unphased_reads_list_new: List[str] = []
    current_all_phased_reads_list_new: List[str] = []
    for index, file_name in tqdm(enumerate(sorted_read_assignment_files)):
        _, _, assignment_windows = parse_assignment_filename(file_name)
        current_block_reads_f_path = os.path.join(assignment_path, file_name)

        if index == len(sorted_read_assignment_files) - 1:
            current_snp_block = block_relationship_dfs[chrom].snp_phased_block_2[index - 1]
        else:
            current_snp_block = block_relationship_dfs[chrom].snp_phased_block_1[index]

        current_block_reads_df = pd.read_csv(current_block_reads_f_path)
        read_to_hp_dict = dict(
            zip(current_block_reads_df["read_id"], current_block_reads_df["haplotype"])
        )
        current_flip_flag = snp_block_flipping_chrom[current_snp_block]

        snp_start, snp_end = _parse_block(current_snp_block)
        current_snp_block_reads = input_bam.fetch(chrom, snp_start, snp_end)
        current_all_phased_reads_list: List[str] = []
        if current_flip_flag == 1:
            for reads in current_snp_block_reads:
                if reads.query_name in current_all_phased_reads_list_new or not reads.has_tag("HP"):
                    continue
                current_all_phased_reads_list.append(reads.query_name)
                write_read(reads)
        elif current_flip_flag == -1:
            for reads in current_snp_block_reads:
                if reads.query_name in current_all_phased_reads_list_new or not reads.has_tag("HP"):
                    continue
                current_all_phased_reads_list.append(reads.query_name)
                if reads.get_tag("HP") == 1:
                    reads.set_tag(tag="HP", value=2, value_type="i")
                elif reads.get_tag("HP") == 2:
                    reads.set_tag(tag="HP", value=1, value_type="i")
                write_read(reads)
        current_all_phased_reads_list_new = current_all_phased_reads_list

        current_all_unphased_reads_list: List[str] = []

        if current_flip_flag == 1:
            for reads in _iter_reads_from_regions(input_bam, chrom, assignment_windows):
                if (
                    reads.query_name in current_all_unphased_reads_list_new
                    or reads.query_name in current_all_unphased_reads_list
                    or reads.has_tag("HP")
                ):
                    continue
                current_all_unphased_reads_list.append(reads.query_name)
                if reads.query_name in read_to_hp_dict:
                    if read_to_hp_dict[reads.query_name] == 1:
                        reads.set_tag(tag="HP", value=1, value_type="i")
                    elif read_to_hp_dict[reads.query_name] == 2:
                        reads.set_tag(tag="HP", value=2, value_type="i")
                write_read(reads)

        elif current_flip_flag == -1:
            for reads in _iter_reads_from_regions(input_bam, chrom, assignment_windows):
                if (
                    reads.query_name in current_all_unphased_reads_list_new
                    or reads.query_name in current_all_unphased_reads_list
                    or reads.has_tag("HP")
                ):
                    continue
                current_all_unphased_reads_list.append(reads.query_name)
                if reads.query_name in read_to_hp_dict:
                    if read_to_hp_dict[reads.query_name] == 1:
                        reads.set_tag(tag="HP", value=2, value_type="i")
                    elif read_to_hp_dict[reads.query_name] == 2:
                        reads.set_tag(tag="HP", value=1, value_type="i")
                write_read(reads)

        current_all_unphased_reads_list_new = current_all_unphased_reads_list

    for reads in input_bam.fetch(chrom):
        write_read(reads)

    modified_bam.close()
    input_bam.close()
    original_bam.close()


def get_block_relationships(
    input_folder: str,
    min_required_read: int = 0,
    min_diff_perc: float = 0.5,
    chromosomes: Sequence[str] | None = None,
) -> Dict[str, pd.DataFrame]:
    relationship_df_by_chr: Dict[str, pd.DataFrame] = {}
    selected = set(chromosomes) if chromosomes else None

    for chr_name in [x for x in os.listdir(input_folder) if "_" not in x]:
        if selected and chr_name not in selected:
            continue

        relationship_df = pd.DataFrame()
        chr_folder_path = os.path.join(input_folder, chr_name)
        if not os.path.isdir(chr_folder_path):
            continue

        for csv_file in sorted(os.listdir(chr_folder_path), key=lambda x: int(x.split("_")[0])):
            csv_file_path = os.path.join(chr_folder_path, csv_file)
            rel_df = pd.read_csv(csv_file_path, index_col=0)[1:].reset_index()
            rel_df = result_filtering(rel_df, min_required_read, min_diff_perc)
            relationship_df = pd.concat([relationship_df, rel_df], ignore_index=True)

        if not relationship_df.empty:
            relationship_df_by_chr[chr_name] = relationship_df

    return relationship_df_by_chr


def run_step2(config: Step2Config) -> List[str]:
    block_relationship_dfs = get_block_relationships(
        config.meth_phasing_input_folder,
        min_required_read=config.minimum_coverage,
        min_diff_perc=config.voting_difference,
        chromosomes=config.chromosomes,
    )
    if not block_relationship_dfs:
        raise ValueError("no relationship CSVs found for step2 post-processing")

    _, remaining_dict, flipping_dict = get_altered_vcf(
        config.vcf_called,
        config.output_vcf,
        block_relationship_dfs,
        include_unprocessed_contigs=config.include_unprocessed_contigs,
    )

    interval_list = []
    generated_bams: List[str] = []
    for chrom in block_relationship_dfs.keys():
        chrom_assignment_path = os.path.join(
            config.meth_phasing_input_folder,
            f"{chrom}_read_assignment",
        )
        chrom_output_bam = f"{config.output_bam_prefix}.{chrom}.methtagged.bam"
        generated_bams.append(chrom_output_bam)
        interval_list.append(
            (
                flipping_dict,
                remaining_dict,
                block_relationship_dfs,
                config.input_bam_file,
                config.input_bam_file,
                chrom_output_bam,
                chrom_assignment_path,
                chrom,
            )
        )

    if interval_list:
        with Pool(min(config.threads, len(interval_list))) as pool:
            pool.starmap(get_altered_bam, interval_list)

    return generated_bams
