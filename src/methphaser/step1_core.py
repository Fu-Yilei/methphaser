import pysam
import re
import warnings
import argparse
import sys
import os
import pandas as pd
from typing import Iterable, Sequence, Tuple

# from tqdm.auto import tqdm
from scipy import stats

from .step1_regions import (
    build_boundary_window_region_plans,
    build_full_span_region_plans,
    encode_windows_for_filename,
    parse_skipping_pairs,
)

warnings.filterwarnings("ignore")


"""
methphaser: phase reads based on methlytion informaiton
@author: Yilei Fu
@Email: yilei.fu@nanoporetech.com, yf20@rice.edu
"""

METHYLATION_IDENTIFIER_0 = ("C", 0, "m")
METHYLATION_IDENTIFIER_1 = ("C", 1, "m")
ReadWindow = Tuple[int, int]
ASSIGNMENT_DF_COLUMNS = [
    "read_id",
    "read_len",
    "phase_block",
    "ref_start",
    "ref_end",
    "haplotype",
    "hp_supportring_cgs",
    "total_cgs",
    "voting",
]
PROBABILITY_DF_COLUMNS = [
    "hp_1_prob",
    "hp_2_prob",
    "p-value",
    "hp_0_prob",
    "assignment",
]


def parse_arg(argv):
    """
    Function for pasing arguments
    """

    parser = argparse.ArgumentParser(
        description="methphaser: phase reads based on methlytion informaiton"
    )
    required_args = parser.add_argument_group("Required arguments")
    # input set
    required_args.add_argument(
        "-b",
        "--bam_file",
        type=str,
        help="input methylation annotated bam file",
        required=True,
        metavar="",
    )
    required_args.add_argument(
        "-r",
        "--reference",
        type=str,
        help="reference genome",
        required=True,
        metavar="",
    )
    required_args.add_argument(
        "-p",
        "--phased_blocks",
        type=str,
        help="gtf file from whatshap visualization",
        required=True,
        metavar="",
    )


    parser.add_argument(
        "-t", "--threads", type=int, help="threads, default 1", default=1, metavar=""
    )
    parser.add_argument(
        "-c",
        "--cut_off",
        type=float,
        help="the minimum percentage of vote to determine a read's haplotype, default 0.65",
        default=0.65,
        metavar="",
    )
    parser.add_argument(
        "-a",
        "--assignment_min",
        type=int,
        help="minimum assigned read number for ranksum test, default 2",
        default=2,
        metavar="",
    )
    parser.add_argument(
        "-m",
        "--chromosome",
        type=str,
        help="the chromosome for read phasing, default chr1",
        default="chr1",
        metavar="",
    )
    parser.add_argument(
        "-n",
        "--targeting_blocks",
        type=str,
        help="only process the blocks from <m, n>, for testing,, default: all",
        default="all",
        metavar="",
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        type=str,
        help="output block-relationship CSV file",
        default="methphaser_relationships.csv",
        metavar="",
    )
    parser.add_argument(
        "-s",
        "--skipping_pair",
        type=str,
        help="a list of number [a, b, ...]. The program will skip the block pair of <a, a+1>, <b, b+1>, ... mainly for avoiding centromere region.",
        default="[0]",
        metavar="",
    )
    parser.add_argument(
        "-k",
        "--k_iterations",
        type=int,
        help="use at most k iterations, default: 10, use -1 for unlimited iterations",
        default=10,
        metavar="",
    )
    parser.add_argument(
        "-ra",
        "--read_assignment",
        type=str,
        help="output read assignment csv folder. The output csv will be folder/phase-block.csv",
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--region_strategy",
        choices=["full-span", "boundary-window"],
        default="full-span",
        help="full-span keeps legacy expanded regions; boundary-window uses SNP-block edge windows",
    )
    parser.add_argument(
        "--boundary_window_bp",
        type=int,
        default=1000,
        help="window size for boundary-window strategy",
        metavar="",
    )

    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)

    return args


def _resolve_phase_windows(
    phase_region: ReadWindow,
    phase_windows: Sequence[ReadWindow] | None = None,
) -> list[ReadWindow]:
    if phase_windows:
        return [(int(start), int(end)) for start, end in phase_windows if int(start) < int(end)]
    start, end = int(phase_region[0]), int(phase_region[1])
    if start >= end:
        return []
    return [(start, end)]


def _phase_region_bounds(
    phase_region: ReadWindow,
    phase_windows: Sequence[ReadWindow] | None = None,
) -> ReadWindow:
    windows = _resolve_phase_windows(phase_region, phase_windows)
    if not windows:
        return phase_region
    return windows[0][0], windows[-1][1]


def _iter_reads_in_windows(
    bam_file: pysam.AlignmentFile,
    chromosome: str,
    phase_region: ReadWindow,
    phase_windows: Sequence[ReadWindow] | None = None,
    multiple_iterators: bool = False,
) -> Iterable[pysam.AlignedSegment]:
    windows = _resolve_phase_windows(phase_region, phase_windows)
    seen_alignment_keys = set()
    for start, end in windows:
        for reads in bam_file.fetch(
            chromosome,
            start,
            end,
            multiple_iterators=multiple_iterators,
        ):
            key = (
                reads.query_name,
                reads.reference_start,
                reads.reference_end,
                reads.is_secondary,
                reads.is_supplementary,
            )
            if key in seen_alignment_keys:
                continue
            seen_alignment_keys.add(key)
            yield reads


def _collect_cpg_positions(
    ref_seq: pysam.FastaFile,
    chromosome: str,
    phase_region: ReadWindow,
    phase_windows: Sequence[ReadWindow] | None = None,
) -> list[int]:
    windows = _resolve_phase_windows(phase_region, phase_windows)
    positions = set()
    for start, end in windows:
        phased_block_ref = ref_seq.fetch(chromosome, start, end)
        cg_loc = [match.start(0) for match in re.finditer("CG", str(phased_block_ref))]
        positions.update(item + start + 1 for item in cg_loc)
    return sorted(positions)


def _pick_methylation_identifier(mm):
    if METHYLATION_IDENTIFIER_0 in mm:
        return METHYLATION_IDENTIFIER_0
    if METHYLATION_IDENTIFIER_1 in mm:
        return METHYLATION_IDENTIFIER_1
    return None


def get_base_modification_dictionary(
    bam_file,
    ref_seq,
    chromosome,
    phase_region,
    snp_phase_region,
    phase_windows=None,
):
    """
    This is for the first time phasing within SNP phased block
    Return value: a dictionary that contains cpg location and its haplotype related base modification score
    """
    cg_loc = _collect_cpg_positions(
        ref_seq,
        chromosome,
        phase_region,
        phase_windows=phase_windows,
    )
    hp_myth_dict = {}
    """
        Data structure:
        {i:[[], [], [], [0, 0, 0]]} 
            i: CpG locations
            []: per haplotype base modification score, max 255
            0: per haplotype base coverage
        Use dictionary so that when querying CpG locations the time complexity is O(1)
    """
    for i in cg_loc:
        # build the dictionary
        hp_myth_dict.update({i: [[], [], [], [0, 0, 0]]})
    phased_block_alignment = bam_file.fetch(
        chromosome, snp_phase_region[0], snp_phase_region[1], multiple_iterators=True
    )

    for reads in phased_block_alignment:
        read_base_ref_loc = reads.get_reference_positions(
            full_length=True
        )  # use full_length=True or the positions won't match
        mm = (
            reads.modified_bases
        )  # mm is a dictionary that contains {score type: [(location, score)]}. score is 255-based
        HP = 0
        if reads.has_tag("HP"):  # update read number list
            if reads.get_tag("HP") == 1:
                HP = 1
                for i in read_base_ref_loc:
                    if i in hp_myth_dict:  # O(1) search
                        hp_myth_dict[i][3][
                            1
                        ] += 1  # increase the per haplotype base coverage on each CpG locaitons
            else:
                HP = 2
                for i in read_base_ref_loc:
                    if i in hp_myth_dict:
                        hp_myth_dict[i][3][2] += 1
        else:
            HP = 0
            for i in read_base_ref_loc:
                if i in hp_myth_dict:
                    hp_myth_dict[i][3][0] += 1
        if (mm != -1) and (mm != {}):  # update base modification scores list
            methylation_identifier = _pick_methylation_identifier(mm)
            if not methylation_identifier:
                continue
            for i in mm[methylation_identifier]:  # Remora only output one type of score: c 1 m/c 0 m, but this part can be improved for other methlyation callers
                if read_base_ref_loc[i[0]]:  # i format: (loc, score)
                    if reads.is_forward:  # cg/gc on forward and reverse reads
                        mm_ref_loc = read_base_ref_loc[i[0]] + 1
                    else:
                        mm_ref_loc = read_base_ref_loc[i[0]]
                    if mm_ref_loc in hp_myth_dict:
                        modification_chance = i[1]  # 0 - 255 based
                        if HP == 1:
                            hp_myth_dict[mm_ref_loc][1].append(
                                modification_chance
                            )  # add the score to the list
                        elif HP == 2:
                            hp_myth_dict[mm_ref_loc][2].append(
                                modification_chance)
                        else:
                            hp_myth_dict[mm_ref_loc][0].append(
                                modification_chance)
    return hp_myth_dict


def get_modified_list(
    assignment_df,
    phase_region,
    sam_file,
    chromosome,
    ref_seq,
    phase_windows=None,
):
    """
    Build the same dictionary as SNP phased region one,
    and update this one during the iterations
    difference between get_base_modification_dictionary:
        different regions, different HP0 filtering
    """
    cg_loc = _collect_cpg_positions(
        ref_seq,
        chromosome,
        phase_region,
        phase_windows=phase_windows,
    )
    hp_myth_dict = {}
    for i in cg_loc:
        hp_myth_dict.update({i: [[], [], [], [0, 0, 0]]})
    assigned_haplotype_by_read = dict(
        zip(assignment_df["read_id"], assignment_df["haplotype"])
    )

    for reads in _iter_reads_in_windows(
        sam_file,
        chromosome,
        phase_region,
        phase_windows=phase_windows,
    ):
        read_base_ref_loc = reads.get_reference_positions(full_length=True)
        mm = reads.modified_bases
        HP = 0

        if reads.query_name in assigned_haplotype_by_read:
            # which means the read was marked as HP 0 by Whatshap
            read_reassign_haplotype = assigned_haplotype_by_read[reads.query_name]
            if read_reassign_haplotype == 1:
                HP = 1
                for i in read_base_ref_loc:
                    if i in hp_myth_dict:
                        hp_myth_dict[i][3][1] += 1
            elif read_reassign_haplotype == 2:
                HP = 2
                for i in read_base_ref_loc:
                    if i in hp_myth_dict:
                        hp_myth_dict[i][3][2] += 1
            else:
                HP = 0
                for i in read_base_ref_loc:
                    if i in hp_myth_dict:
                        hp_myth_dict[i][3][0] += 1

            if mm != -1 and mm != {}:
                methylation_identifier = _pick_methylation_identifier(mm)
                if not methylation_identifier:
                    continue
                for i in mm[methylation_identifier]:
                    if read_base_ref_loc[i[0]]:
                        if reads.is_forward:
                            mm_ref_loc = read_base_ref_loc[i[0]] + 1
                        else:
                            mm_ref_loc = read_base_ref_loc[i[0]]
                        if mm_ref_loc in hp_myth_dict:
                            modification_chance = i[1]
                            if HP == 1:
                                hp_myth_dict[mm_ref_loc][1].append(
                                    modification_chance)
                            elif HP == 2:
                                hp_myth_dict[mm_ref_loc][2].append(
                                    modification_chance)
                            else:
                                hp_myth_dict[mm_ref_loc][0].append(
                                    modification_chance)
    return hp_myth_dict


def build_df(previous_assignment_df, hp_list):
    if not hp_list:
        return previous_assignment_df

    incoming_df = pd.DataFrame(hp_list, columns=list(previous_assignment_df.columns))
    existing_ids = set(previous_assignment_df["read_id"])

    append_df = incoming_df[~incoming_df["read_id"].isin(existing_ids)]
    if not append_df.empty:
        if previous_assignment_df.empty:
            previous_assignment_df = append_df.reset_index(drop=True)
        else:
            previous_assignment_df = pd.concat(
                [previous_assignment_df, append_df], ignore_index=True
            )

    update_df = incoming_df[incoming_df["read_id"].isin(existing_ids)]
    if not update_df.empty:
        indexed_existing = previous_assignment_df.set_index("read_id")
        indexed_update = update_df.set_index("read_id")
        indexed_existing.update(indexed_update)
        previous_assignment_df = indexed_existing.reset_index()

    return previous_assignment_df


def _build_read_probability_df(reads, hp_base_modification_probablity):
    read_base_ref_loc = reads.get_reference_positions(full_length=True)
    read_base_ref_loc_aligned = reads.get_reference_positions(full_length=False)
    mm = reads.modified_bases
    if mm == -1 or mm == {}:
        return None, read_base_ref_loc_aligned

    methylation_identifier = _pick_methylation_identifier(mm)
    if not methylation_identifier:
        return None, read_base_ref_loc_aligned

    hp_0_probablity = {}
    for i in mm[methylation_identifier]:
        if read_base_ref_loc[i[0]]:
            if reads.is_forward:
                mm_ref_loc = read_base_ref_loc[i[0]] + 1
            else:
                mm_ref_loc = read_base_ref_loc[i[0]]

            if mm_ref_loc in hp_base_modification_probablity:
                modification_chance = i[1] / 255
                hp_1_prob = hp_base_modification_probablity[mm_ref_loc][0]
                hp_2_prob = hp_base_modification_probablity[mm_ref_loc][1]
                assignment = 0
                if hp_1_prob is not None and hp_2_prob is not None:
                    if abs(hp_1_prob - modification_chance) > abs(
                        hp_2_prob - modification_chance
                    ):
                        assignment = 2
                    elif abs(hp_1_prob - modification_chance) < abs(
                        hp_2_prob - modification_chance
                    ):
                        assignment = 1
                probablity_result = hp_base_modification_probablity[mm_ref_loc] + [
                    modification_chance,
                    assignment,
                ]
                hp_0_probablity.update({mm_ref_loc: probablity_result})

    assignment_df = pd.DataFrame.from_dict(
        hp_0_probablity,
        orient="index",
        columns=PROBABILITY_DF_COLUMNS,
    )
    return assignment_df, read_base_ref_loc_aligned


def _append_assignment_record(
    assignment_list,
    reads,
    phase_region_start,
    assignment_df,
    assignmet_threshold,
    assignment_min_number,
    forced_haplotype=None,
):
    read_base_ref_loc_aligned = reads.get_reference_positions(full_length=False)
    ref_start = read_base_ref_loc_aligned[0] if read_base_ref_loc_aligned else None
    ref_end = read_base_ref_loc_aligned[-1] if read_base_ref_loc_aligned else None
    read_length = reads.query_length

    if len(assignment_df) == 0:
        assignment_list.append(
            [
                reads.query_name,
                read_length,
                phase_region_start,
                ref_start,
                ref_end,
                forced_haplotype if forced_haplotype is not None else 0,
                None,
                0,
                None,
            ]
        )
        return

    hp1_num = len(assignment_df[assignment_df.assignment == 1])
    hp2_num = len(assignment_df[assignment_df.assignment == 2])
    winner_num = hp1_num if hp1_num >= hp2_num else hp2_num
    winner_hap = 1 if hp1_num >= hp2_num else 2
    winner_ratio = winner_num / len(assignment_df)

    if forced_haplotype is not None:
        assignment_list.append(
            [
                reads.query_name,
                read_length,
                phase_region_start,
                ref_start,
                ref_end,
                forced_haplotype,
                winner_num,
                len(assignment_df),
                winner_ratio,
            ]
        )
        return

    if len(assignment_df) <= assignment_min_number:
        assignment_list.append(
            [
                reads.query_name,
                read_length,
                phase_region_start,
                ref_start,
                ref_end,
                0,
                winner_num,
                len(assignment_df),
                winner_ratio,
            ]
        )
        return

    hp1_ratio = hp1_num / len(assignment_df)
    hp2_ratio = hp2_num / len(assignment_df)
    if hp1_ratio >= assignmet_threshold:
        assignment_list.append(
            [
                reads.query_name,
                read_length,
                phase_region_start,
                ref_start,
                ref_end,
                1,
                hp1_num,
                len(assignment_df),
                hp1_ratio,
            ]
        )
        return
    if hp2_ratio >= assignmet_threshold:
        assignment_list.append(
            [
                reads.query_name,
                read_length,
                phase_region_start,
                ref_start,
                ref_end,
                2,
                hp2_num,
                len(assignment_df),
                hp2_ratio,
            ]
        )
        return

    assignment_list.append(
        [
            reads.query_name,
            read_length,
            phase_region_start,
            ref_start,
            ref_end,
            0,
            winner_num,
            len(assignment_df),
            winner_ratio,
        ]
    )


def _assign_reads_by_methylation(
    bam_file,
    chromosome,
    phase_region,
    hp_base_modification_probablity,
    previous_assignment_df,
    assignmet_threshold,
    assignment_min_number,
    read_filter,
    phase_windows=None,
    force_haplotype_zero_for_low_support=False,
):
    phase_region_start = _phase_region_bounds(phase_region, phase_windows=phase_windows)[0]
    assignment_list = []
    for reads in _iter_reads_in_windows(
        bam_file,
        chromosome,
        phase_region,
        phase_windows=phase_windows,
        multiple_iterators=True,
    ):
        if not read_filter(reads):
            continue
        assignment_df, _ = _build_read_probability_df(
            reads,
            hp_base_modification_probablity,
        )
        if assignment_df is None:
            continue
        _append_assignment_record(
            assignment_list,
            reads,
            phase_region_start,
            assignment_df,
            assignmet_threshold,
            assignment_min_number,
            forced_haplotype=0 if force_haplotype_zero_for_low_support else None,
        )
    return build_df(previous_assignment_df, assignment_list)


def get_base_modification_list(
    bam_file,
    ref_seq,
    chromosome,
    phase_region,
    hp_base_modification_probablity,
    previous_assignment_df,
    assignmet_threshold,
    assignment_min_number,
    phase_windows=None,
):
    """
    Assign previously-unresolved reads during iterative refinement.
    """
    del ref_seq
    unresolved_read_ids = set(
        previous_assignment_df[previous_assignment_df.haplotype == 0]["read_id"]
    )
    return _assign_reads_by_methylation(
        bam_file=bam_file,
        chromosome=chromosome,
        phase_region=phase_region,
        hp_base_modification_probablity=hp_base_modification_probablity,
        previous_assignment_df=previous_assignment_df,
        assignmet_threshold=assignmet_threshold,
        assignment_min_number=assignment_min_number,
        read_filter=lambda reads: reads.query_name in unresolved_read_ids,
        phase_windows=phase_windows,
        force_haplotype_zero_for_low_support=False,
    )


def get_base_modification_list_snp_block(
    bam_file,
    ref_seq,
    chromosome,
    phase_region,
    hp_base_modification_probablity,
    previous_assignment_df,
    assignmet_threshold,
    assignment_min_number,
    phase_windows=None,
):
    """
    Assign initial SNP-unphased reads before iterative refinement.
    """
    del ref_seq
    return _assign_reads_by_methylation(
        bam_file=bam_file,
        chromosome=chromosome,
        phase_region=phase_region,
        hp_base_modification_probablity=hp_base_modification_probablity,
        previous_assignment_df=previous_assignment_df,
        assignmet_threshold=assignmet_threshold,
        assignment_min_number=assignment_min_number,
        read_filter=lambda reads: not reads.has_tag("HP"),
        phase_windows=phase_windows,
        force_haplotype_zero_for_low_support=False,
    )


def get_siginificant_probablity_dict(hp_myth_dict, hp_myth_list_new=None, hp_min_num=3):
    """
    This function takes base modificaiton score list as input, output a dict that
    also contains base modification probablity with SNP unphased reads
    In this function, we perform ranksum test

    data structure: {i:[hp_1_prob, hp_2_prob, p-value]}
        i: CpG locations
        hp_1/2_prob = sum(base modification score)/(per-base per-haplotype coverage)
        p-value: calculated with scipy package


    args:
        hp_myth_list_new: the base modificaiton list from last iteration
        hp_myth_dict: format
            {
                location: ([[hp0_scores], [hp1_scores], [hp2_scores], [hp0_read_n, hp1_read_n, hp2_read_n]] )
            }
    outputs:
        hp_base_modification_probablity: a dictionary that contians the locations where the 2 hps have siginificant
        different behaves

    """
    if hp_myth_list_new is None:
        hp_myth_list_new = {}

    hp_base_modification_probablity = dict()
    for index, (location, modifications) in enumerate(hp_myth_dict.items()):  # see args
        # the aggregated score of hp1
        hp_1_probablity_sum = sum(modifications[1])
        hp_2_probablity_sum = sum(modifications[2])
        hps = [[], [], []]
        if (
            location in hp_myth_list_new.keys()
        ):  # new and old should have the same index (CpG locations)
            hps = hp_myth_list_new[location]
            # the aggregated score of hp1
            hp_1_probablity_sum_new = sum(hps[1])
            hp_1_probablity_sum += hp_1_probablity_sum_new  # old + new
            hp_2_probablity_sum_new = sum(hps[2])
            hp_2_probablity_sum += hp_2_probablity_sum_new
            # hp_1_new_num = hps[3][1] # new read number
            # hp_2_new_num = hps[3][2]
        # hp_1_num = modifications[3][1] + hp_1_new_num
        # hp_2_num = modifications[3][2] + hp_2_new_num
        h1_full_list = modifications[1] + hps[1]
        h2_full_list = modifications[2] + hps[2]
        if (len(h1_full_list) >= hp_min_num) and (
            len(h2_full_list) >= hp_min_num
        ):  # enough hp info in this location
            hp_1_probablity = hp_1_probablity_sum / (255 * len(h1_full_list))
            hp_2_probablity = hp_2_probablity_sum / (255 * len(h2_full_list))
            ttest_value = stats.ranksums(
                h1_full_list, h2_full_list
            )  # Do a ranksums test of two hp score lists
            if (
                ttest_value.pvalue < 0.05
            ):  # if two score sets are significantly different
                # if abs(hp_1_probablity-hp_2_probablity) >= 0.5:
                hp_base_modification_probablity.update(
                    {
                        location: [
                            hp_1_probablity,
                            hp_2_probablity,
                            ttest_value.pvalue,
                            # abs(hp_1_probablity-hp_2_probablity)
                        ]
                    }
                )

    return hp_base_modification_probablity


def get_assignment_max(
    chromosome,
    tagged_bam,
    ref_seq,
    phased_region_l,
    snp_phased_region_l,
    hp_threshold,
    assignment_threshold,
    k_iterations,
    phase_windows_l=None,
):
    '''
        The main function of the program: assign reads and assign relationships between phased blocks.
    '''
    all_not_assigned_reads_increasing_dict = {}
    hp0_assignment_df_dict = {}
    modif_d = {}
    hp_prob_d = {}
    assignment_d = {}
    if k_iterations == -1:
        k_iterations = float("inf")
    for index, phased_regions in enumerate(phased_region_l):
        current_phase_windows = None
        if phase_windows_l is not None:
            current_phase_windows = phase_windows_l[index]

        # test
        modif_l = []
        hp_prob_l = []
        assignment_list = []

        not_assigned_reads_list = []
        increasing_assigned_num = float("inf")
        cnt = 0
        hp0_cnt = float("inf")
        hp0_cnt_new = 0
        modified_list = {}
        assignment_df = pd.DataFrame(columns=ASSIGNMENT_DF_COLUMNS)
        # only add snp phased region first

        cnt += 1
        snp_phased_region = snp_phased_region_l[index]  # get snp phased block
        # print(
        #     f"processing phased region:{snp_phased_region[0]}-{snp_phased_region[1]}, iteration {cnt}",
        #     end="\r",
        # )
        base_modification_list = get_base_modification_dictionary(  # build the dictionary with snp phased reads
            tagged_bam,
            ref_seq,
            chromosome,
            phased_regions,
            snp_phased_region,
            phase_windows=current_phase_windows,
        )
        hp_base_modification_prob = get_siginificant_probablity_dict(  # calculate the statistically significantly different CpGs
            base_modification_list, modified_list, hp_min_num=hp_threshold
        )
        assignment_df = get_base_modification_list_snp_block(
            tagged_bam,
            ref_seq,
            chromosome,
            phased_regions,
            hp_base_modification_prob,
            assignment_df,
            assignment_threshold,
            hp_threshold,
            phase_windows=current_phase_windows,
        )

        modified_list = get_modified_list(
            assignment_df,
            phased_regions,
            tagged_bam,
            chromosome,
            ref_seq,
            phase_windows=current_phase_windows,
        )
        hp0_cnt_new = len(assignment_df[assignment_df.haplotype == 0])
        increasing_assigned_num = hp0_cnt_new - hp0_cnt
        hp0_cnt = hp0_cnt_new
        not_assigned_reads_list.append(hp0_cnt)

        while (increasing_assigned_num < 0) and (cnt <= k_iterations):
            '''
                start the iteration
            '''
            modif_l.append(modified_list)
            hp_prob_l.append(hp_base_modification_prob)
            assignment_list.append(assignment_df)

            cnt += 1
            # print(
            #     f"processing unphased region:{phased_regions[0]}-{phased_regions[1]}, iteration {cnt}",
            #     end="\r",
            # )

            hp_base_modification_prob = get_siginificant_probablity_dict(
                base_modification_list, modified_list, hp_min_num=hp_threshold
            )
            assignment_df = get_base_modification_list(
                tagged_bam,
                ref_seq,
                chromosome,
                phased_regions,
                hp_base_modification_prob,
                assignment_df,
                assignment_threshold,
                hp_threshold,
                phase_windows=current_phase_windows,
            )
            modified_list = get_modified_list(
                assignment_df,
                phased_regions,
                tagged_bam,
                chromosome,
                ref_seq,
                phase_windows=current_phase_windows,
            )
            hp0_cnt_new = len(assignment_df[assignment_df.haplotype == 0])
            increasing_assigned_num = hp0_cnt_new - hp0_cnt
            # print(
            #     f"increasing hp assignment: {hp0_cnt_new} - {hp0_cnt} = {increasing_assigned_num}",
            #     end="\r",
            # )
            hp0_cnt = hp0_cnt_new
            not_assigned_reads_list.append(hp0_cnt)

        assignment_list.append(assignment_df)
        modif_l.append(modified_list)
        hp_prob_l.append(hp_base_modification_prob)
        modif_d.update({index: modif_l})
        hp_prob_d.update({index: hp_prob_l})
        assignment_d.update({index: assignment_list})

        all_not_assigned_reads_increasing_dict.update(
            {phased_regions: not_assigned_reads_list}
        )
        hp0_assignment_df_dict.update({phased_regions: assignment_df})

    return (
        all_not_assigned_reads_increasing_dict,
        hp0_assignment_df_dict,
        modif_d,
        hp_prob_d,
        assignment_d,
    )


'''
Benchmarking functions
'''


def get_overlap_reads_df(region_1, region_2, assignment_df):
    df1 = assignment_df[region_1]
    df2 = assignment_df[region_2]
    df2_overlap = df2[df2.read_id.isin(df1.read_id)].copy()
    another_hap = []
    for i in df2_overlap.read_id:
        another_hap.append(df1[df1.read_id == i].iloc[0].haplotype)
    df2_overlap.loc[:, "haplotype_in_connected_block"] = another_hap
    overlap_df = df2_overlap[
        (df2_overlap.haplotype != 0) & (
            df2_overlap.haplotype_in_connected_block != 0)
    ]
    same_assignment_num = len(
        overlap_df[overlap_df.haplotype ==
                   overlap_df.haplotype_in_connected_block]
    )
    diff_assignment_num = len(
        overlap_df[overlap_df.haplotype !=
                   overlap_df.haplotype_in_connected_block]
    )
    return overlap_df, same_assignment_num, diff_assignment_num


def output_block_csv(chromosome, region_start_index, snp_phased_region_list, read_assignmentg_df, output_csv):
    result_comparison_vcf = pd.DataFrame(
        columns=[
            "snp_phased_block_1",
            "snp_phased_block_2",
            "extended_phased_block_1",
            "extended_phased_block_2",
            # "vcf_file_relationship",
            "myth_phasing_relationship",
            "same_hap_num",
            "diff_hap_num",
            "olp_hp1_total_read_len",
            "olp_hp2_total_read_len",
            "olp_hp1_CpG_num",
            "olp_hp2_CpG_num",
            "olp_same_read_len",
            "olp_not_same_read_len",
            "olp_same_CpG_num",
            "olp_not_same_CpG_num",
        ]
    )
    for index, phased_region_item in enumerate(list(read_assignmentg_df.keys())):
        result_comparison_l = []
        # result_hp_df = read_assignmentg_df[phased_region_item]
        # previous_result_hp_df_index_snp = phased_region_list[index-1]
        previous_result_hp_df_index = list(
            read_assignmentg_df.keys())[index - 1]
        phased_snp_region_item = list(read_assignmentg_df.keys())[index]
        overlaps_reads_df = get_overlap_reads_df(
            previous_result_hp_df_index, phased_region_item, read_assignmentg_df
        )[0]
        phased_region_1 = snp_phased_region_list[region_start_index + index - 1]
        phased_region_2 = snp_phased_region_list[region_start_index + index]
        result_comparison_l += [  # display the regions
            phased_region_1,
            phased_region_2,
            previous_result_hp_df_index,
            phased_snp_region_item,
        ]

        # if (vcf_relationship_list[region_start_index + index - 1][1] == None) or (
        #     vcf_relationship_list[region_start_index + index][1] == None
        # ):
        #     result_comparison_l += ["cannot decide"]
        # elif (
        #     vcf_relationship_list[region_start_index + index - 1][1]
        #     == vcf_relationship_list[region_start_index + index][1]
        # ):
        #     result_comparison_l += ["same"]
        # else:
        #     result_comparison_l += ["not same"]

        if len(overlaps_reads_df) > 0:
            same_hap_num = len(
                overlaps_reads_df[
                    overlaps_reads_df.haplotype
                    == overlaps_reads_df.haplotype_in_connected_block
                ]
            )

            diff_hap_num = len(
                overlaps_reads_df[
                    overlaps_reads_df.haplotype
                    != overlaps_reads_df.haplotype_in_connected_block
                ]
            )

            if same_hap_num > diff_hap_num:
                # 'meth_phasing_relationship' 
                result_comparison_l += ["same", same_hap_num, diff_hap_num]
            elif same_hap_num < diff_hap_num:
                result_comparison_l += ["not same", same_hap_num, diff_hap_num]
            else:
                result_comparison_l += ["cannot decide",
                                        same_hap_num, diff_hap_num]
        else:
            result_comparison_l += ["cannot decide", None, None]
        total_hp_1_read_len = sum(
            overlaps_reads_df[overlaps_reads_df.haplotype == 1].read_len
        )
        total_hp_2_read_len = sum(
            overlaps_reads_df[overlaps_reads_df.haplotype == 2].read_len
        )
        total_same_read_len = sum(
            overlaps_reads_df[
                overlaps_reads_df.haplotype
                == overlaps_reads_df.haplotype_in_connected_block
            ].read_len
        )
        total_not_same_read_len = sum(
            overlaps_reads_df[
                overlaps_reads_df.haplotype
                != overlaps_reads_df.haplotype_in_connected_block
            ].read_len
        )
        total_same_CpG_num = sum(
            overlaps_reads_df[
                overlaps_reads_df.haplotype
                == overlaps_reads_df.haplotype_in_connected_block
            ].total_cgs
        )
        total_not_same_CpG_num = sum(
            overlaps_reads_df[
                overlaps_reads_df.haplotype
                != overlaps_reads_df.haplotype_in_connected_block
            ].total_cgs
        )
        total_hp_1_CpG_num = sum(
            overlaps_reads_df[overlaps_reads_df.haplotype == 1].total_cgs
        )
        total_hp_2_CpG_num = sum(
            overlaps_reads_df[overlaps_reads_df.haplotype == 2].total_cgs
        )
        result_comparison_l += [
            total_hp_1_read_len,
            total_hp_2_read_len,
            total_hp_1_CpG_num,
            total_hp_2_CpG_num,
            total_same_read_len,
            total_not_same_read_len,
            total_same_CpG_num,
            total_not_same_CpG_num,
        ]
        result_comparison_vcf.loc[len(
            result_comparison_vcf.index)] = result_comparison_l  # type: ignore
    result_comparison_vcf.to_csv(output_csv)

    return

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # parse args
    args = parse_arg(argv)
    bam_file = args.bam_file
    ref_seq = pysam.FastaFile(args.reference)
    phased_block_df = pd.read_csv(
        args.phased_blocks,
        header=None,
        sep="\t",
        names=[
            "chr",
            "phasing",
            "ex/intron",
            "start",
            "end",
            "1",
            "strand",
            "2",
            "info",
        ],
    )
    skipping_block_pair_str = args.skipping_pair
    threads = args.threads
    tagged_bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)
    chromosome = args.chromosome
    phased_block_df = phased_block_df[phased_block_df.chr == chromosome]
    n_blocks = args.targeting_blocks
    hp_min_number = args.assignment_min
    minimum_voting = args.cut_off
    read_assignment_output = args.read_assignment
    region_strategy = args.region_strategy
    boundary_window_bp = max(1, args.boundary_window_bp)
    targetd_blocks = []
    k_iterations = args.k_iterations
    if n_blocks == "all":
        targetd_blocks = [0, len(phased_block_df)]
    else:
        targetd_blocks = [int(x) for x in n_blocks.split(",")]
    # print(targetd_blocks)
    # if targetd_blocks[1]-targetd_blocks[0] == 1:
    #     print(targetd_blocks)

    output_csv = args.output_csv
    # initial phasing regions
    phased_region_list = list(
        zip(phased_block_df.start, phased_block_df.end)
    )  # snp phased blocks
    if len(phased_region_list) < 2:
        pd.DataFrame(
            columns=[
                "snp_phased_block_1",
                "snp_phased_block_2",
                "extended_phased_block_1",
                "extended_phased_block_2",
                "myth_phasing_relationship",
                "same_hap_num",
                "diff_hap_num",
                "olp_hp1_total_read_len",
                "olp_hp2_total_read_len",
                "olp_hp1_CpG_num",
                "olp_hp2_CpG_num",
                "olp_same_read_len",
                "olp_not_same_read_len",
                "olp_same_CpG_num",
                "olp_not_same_CpG_num",
            ]
        ).to_csv(output_csv)
        return
    skipping_block_pair_start, _ = parse_skipping_pairs(skipping_block_pair_str)
    if region_strategy == "full-span":
        region_plans = build_full_span_region_plans(
            phased_region_list,
            skipping_block_pair_start,
        )
    else:
        region_plans = build_boundary_window_region_plans(
            phased_region_list,
            skipping_block_pair_start,
            boundary_window_bp=boundary_window_bp,
        )

    selected_region_plans = region_plans[targetd_blocks[0]: targetd_blocks[1]]
    selected_phase_regions = [plan.span for plan in selected_region_plans]
    selected_phase_windows = [list(plan.windows) for plan in selected_region_plans]

    (
        assigned_reads_increasing_list,
        unphased_reads_assignment_dataframes,
        modification_dict,
        hp_prob_dict,
        assignment_dict,
    ) = get_assignment_max(
        chromosome,
        tagged_bam,
        ref_seq,
        selected_phase_regions,
        phased_region_list[targetd_blocks[0]: targetd_blocks[1]],
        hp_min_number,
        minimum_voting,
        k_iterations,
        phase_windows_l=selected_phase_windows,
    )
    output_block_csv(
        chromosome,
        targetd_blocks[0],
        phased_region_list,
        unphased_reads_assignment_dataframes,
        output_csv,
    )
    if read_assignment_output:
        os.makedirs(read_assignment_output, exist_ok=True)
        for index, phase_block_name in enumerate(unphased_reads_assignment_dataframes.keys()):
            current_block_num = index + targetd_blocks[0]
            window_suffix = encode_windows_for_filename(selected_phase_windows[index])
            output_name = f"{current_block_num}_{phase_block_name[0]}_{phase_block_name[1]}"
            if window_suffix:
                output_name = f"{output_name}_{window_suffix}"
            unphased_reads_assignment_dataframes[phase_block_name].to_csv(
                os.path.join(read_assignment_output, f"{output_name}.csv")
            )


if __name__ == "__main__":
    main(sys.argv[1:])
