from pathlib import Path

import pytest

from methphaser.pipeline import PipelineConfig, run_pipeline


TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data"
REFERENCE = DATA_DIR / "reference" / "GRCh38.chr21.fa"
BAM = DATA_DIR / "HG002_2025q1_PAW70337_30x.chr21_4.9M_6.2M.subset.bam"
GTF = DATA_DIR / "PAW70337_output" / "SAMPLE.wf_snp.haploblocks.gtf"
VCF = DATA_DIR / "PAW70337_output" / "SAMPLE.wf_snp.vcf.gz"

REQUIRED_FILES = [
    REFERENCE,
    BAM,
    BAM.with_suffix(BAM.suffix + ".bai"),
    GTF,
    VCF,
    Path(str(VCF) + ".tbi"),
]


@pytest.mark.integration
def test_hg002_smoke_pipeline(tmp_path: Path) -> None:
    missing = [str(path) for path in REQUIRED_FILES if not path.exists()]
    assert not missing, f"missing integration test files: {missing}"

    output_dir = tmp_path / "smoke_run"
    config = PipelineConfig(
        bam=str(BAM),
        reference=str(REFERENCE),
        output_dir=str(output_dir),
        threads=1,
    )

    outputs = run_pipeline(config)

    assert Path(outputs["work_dir"]).exists()
    assert Path(outputs["output_vcf"]).exists()
    assert Path(outputs["output_bam"]).exists()
    assert Path(outputs["output_bam_index"]).exists()

    rel_csvs = list((output_dir / "work" / "chr21").glob("*.csv"))
    assert rel_csvs, "expected relationship CSV output for chr21"
