import os
import uuid
import subprocess
import pandas as pd


def depth(
    input_bam: str,
    output_path: str,
    region_str: str = None,
    only_forward: bool = False,
    only_reverse: bool = False,
):
    """
    Run `samtools depth` on a given `input_bam`;

    """
    if only_forward and only_reverse:
        raise ValueError("`only_forward` and `only_reverse` are mutually exclusive.")

    cmd = "samtools depth"
    cmd += " -aa"
    if region_str is not None:
        cmd += f" -r {region_str}"
    if only_forward:
        cmd += " --excl-flags 0x10"
    if only_reverse:
        cmd += " --incl-flags 0x10"
    cmd += f" -o {output_path}"
    cmd += f" {input_bam}"

    subprocess.run(cmd, shell=True, check=True)


def bedcov(
    bam_path: str, bed_path: str, output_csv: str, cov_threshold: int = 100
) -> None:
    """
    Run `samtools bedcov`, and munge results into a CSV

    """

    # We choose to first write to temporary BED file
    # which is reloaded, columns are named, then saved as CSV
    temp_bed = f"{bam_path[:-4]}.temp.{str(uuid.uuid4())[:8]}.bed"

    cmd = f"samtools bedcov -d {cov_threshold} -c"
    cmd += f" {bed_path} {bam_path} > {temp_bed}"
    subprocess.run(cmd, shell=True, check=True)

    df = pd.read_csv(temp_bed, sep="\t", header=None)
    breadth_col = f"breadth_{cov_threshold}X_bp"
    df.columns = ["chrom", "start", "end", "name", "total_cov", breadth_col, "n_reads"]

    df.insert(3, "length_bp", df["end"] - df["start"])
    df.insert(5, "mean_cov", df["total_cov"] / df["length_bp"])
    df.insert(
        7, breadth_col.replace("_bp", "_per"), 100 * df[breadth_col] / df["length_bp"]
    )
    df.to_csv(output_csv, index=False)

    os.remove(temp_bed)


def mpileup(
    input_bam: str,
    ref_fasta: str,
    output_pileup: str,
    bed_path: str = None,
    region_str: str = None,
):
    """
    Run `samtools mpileup`

    Note:
    - Setting '-Q 0' as common use is evaluating sequencing
    performance

    """
    if bed_path is not None and region_str is not None:
        raise ValueError("`bed_path` and `region_str` are mutually exclusive.")

    cmd = "samtools mpileup"
    cmd += " -a"
    cmd += " -Q 0"  # Min. Quality
    cmd += f" -f {ref_fasta}"
    if bed_path is not None:
        cmd += f" -l {bed_path}"
    if region_str is not None:
        cmd += f" -r {region_str}"
    cmd += f" -o {output_pileup}"
    cmd += f" {input_bam}"

    subprocess.run(cmd, shell=True, check=True)
