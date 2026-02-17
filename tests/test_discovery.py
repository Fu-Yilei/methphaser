from pathlib import Path

from methphaser.discovery import discover_phase_inputs, infer_run_token


def test_infer_run_token_from_bam_name() -> None:
    bam = Path("/tmp/HG002_2025q1_PAW70337_30x.bam")
    assert infer_run_token(bam) == "PAW70337"


def test_discover_phase_inputs_auto(tmp_path: Path) -> None:
    bam = tmp_path / "HG002_2025q1_PAW70337_30x.bam"
    bam.write_text("", encoding="utf-8")

    out_dir = tmp_path / "PAW70337_output"
    out_dir.mkdir()
    gtf = out_dir / "SAMPLE.wf_snp.haploblocks.gtf"
    vcf = out_dir / "SAMPLE.wf_snp.vcf.gz"
    gtf.write_text("gtf", encoding="utf-8")
    vcf.write_text("vcf", encoding="utf-8")

    discovered = discover_phase_inputs(str(bam))
    assert discovered.gtf == gtf.resolve()
    assert discovered.vcf == vcf.resolve()


def test_discover_phase_inputs_explicit(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("", encoding="utf-8")
    gtf = tmp_path / "blocks.gtf"
    vcf = tmp_path / "calls.vcf.gz"
    gtf.write_text("gtf", encoding="utf-8")
    vcf.write_text("vcf", encoding="utf-8")

    discovered = discover_phase_inputs(str(bam), gtf=str(gtf), vcf=str(vcf))
    assert discovered.gtf == gtf.resolve()
    assert discovered.vcf == vcf.resolve()
