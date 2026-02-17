from pathlib import Path

import pandas as pd

from methphaser.imprinting import (
    ParentAssignmentConfig,
    annotate_phaseblocks,
    default_imprinting_bed_path,
    extract_directional_regions,
    extract_phaseblocks_from_vcf,
    load_imprinting_regions,
    run_parent_assignment,
)


def test_extract_directional_regions_filters_unknown() -> None:
    regions = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "start": 100,
                "end": 200,
                "region_id": "A",
                "source_track": "known_ICR_25",
                "source_database": "HumanICR",
                "oocyte_mean": 0.90,
                "sperm_mean": 0.10,
                "methylated_parent": "maternal",
                "direction_confidence": 0.80,
                "oocyte_cpg_sites": 10,
                "sperm_cpg_sites": 10,
                "paternal_haplotype_rule": "low_methylation_haplotype",
            },
            {
                "chrom": "chr1",
                "start": 300,
                "end": 400,
                "region_id": "B",
                "source_track": "known_ICR_25",
                "source_database": "HumanICR",
                "oocyte_mean": 0.5,
                "sperm_mean": 0.5,
                "methylated_parent": "unknown",
                "direction_confidence": 0.0,
                "oocyte_cpg_sites": 2,
                "sperm_cpg_sites": 2,
                "paternal_haplotype_rule": "unknown",
            },
        ]
    )

    directional = extract_directional_regions(regions)
    assert len(directional) == 1
    assert directional.iloc[0]["region_id"] == "A"


def test_extract_phaseblocks_from_vcf(tmp_path: Path) -> None:
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=1000000>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:100",
                "chr1\t150\t.\tC\tT\t.\tPASS\t.\tGT:PS\t1|0:100",
                "chr1\t500\t.\tG\tA\t.\tPASS\t.\tGT:PS\t0|1:500",
            ]
        ),
        encoding="utf-8",
    )

    phaseblocks = extract_phaseblocks_from_vcf(str(vcf_path))
    assert len(phaseblocks) == 2
    assert list(phaseblocks["phase_set"]) == [100, 500]
    assert list(phaseblocks["variant_count"]) == [2, 1]


def test_run_parent_assignment_without_bam_is_random(tmp_path: Path) -> None:
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=1000000>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:100",
                "chr1\t120\t.\tC\tT\t.\tPASS\t.\tGT:PS\t1|0:100",
                "chr1\t900\t.\tG\tA\t.\tPASS\t.\tGT:PS\t0|1:900",
            ]
        ),
        encoding="utf-8",
    )

    bed_path = tmp_path / "imprinting.bed"
    bed_path.write_text(
        "\n".join(
            [
                "# chrom\tstart\tend\tregion_id\tsource_track\tsource_database\toocyte_mean\tsperm_mean\tmethylated_parent\tdirection_confidence\toocyte_cpg_sites\tsperm_cpg_sites\tpaternal_haplotype_rule",
                "chr1\t80\t150\tTEST_1\tknown_ICR_25\tHumanICR\t0.9000\t0.1000\tmaternal\t0.8000\t20\t20\tlow_methylation_haplotype",
            ]
        ),
        encoding="utf-8",
    )

    out_path = tmp_path / "phaseblock_parent_assignment.tsv"
    outputs = run_parent_assignment(
        ParentAssignmentConfig(
            output_vcf=str(vcf_path),
            phaseblock_assignment_output=str(out_path),
            imprinting_bed=str(bed_path),
            input_bam=None,
        )
    )

    assert Path(outputs["phaseblock_parent_assignment"]).exists()
    out_df = pd.read_csv(out_path, sep="\t")
    assert list(out_df["parent_assignment_status"]) == ["random", "random"]
    assert out_df.iloc[0]["assignment_basis"] == "directional_regions_found_but_no_bam_evidence"
    assert out_df.iloc[1]["assignment_basis"] == "no_directional_imprinting_region_overlap"


def test_annotate_phaseblocks_empty_phaseblocks() -> None:
    empty_phaseblocks = pd.DataFrame(
        columns=["chrom", "phase_set", "start", "end", "variant_count"]
    )
    directional_regions = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "start": 100,
                "end": 200,
                "region_id": "R1",
                "source_track": "x",
                "source_database": "y",
                "methylated_parent": "maternal",
                "direction_confidence": 0.8,
            }
        ]
    )
    output = annotate_phaseblocks(
        phaseblocks=empty_phaseblocks,
        directional_regions=directional_regions,
        input_bam=None,
        min_sites_per_haplotype=3,
        min_methylation_diff=0.1,
    )
    assert output.empty


def test_default_imprinting_database_exists_and_loads() -> None:
    db_path = default_imprinting_bed_path()
    assert db_path.exists()
    regions = load_imprinting_regions()
    assert not regions.empty
    assert "methylated_parent" in regions.columns
    assert regions["methylated_parent"].isin(["maternal", "paternal", "unknown"]).all()
