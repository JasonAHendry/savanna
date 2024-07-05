import os
import click
from pathlib import Path


BARCODING_KITS = ["SQK-NBD114-96", "SQK-RBK114-96"]
PIPELINES = ["plasmo", "ento"]


slurm_basecall = Path("slurm/templates/savanna-basecall.slurm")
slurm_analyse = Path("slurm/templates/savanna-analyse.slurm")
slurm_summary = Path("slurm/templates/savanna-summary.slurm")
slurm_submit = Path("slurm/templates/savanna-submit.sh")
RUNS_DIR = Path("slurm/runs")


def load_format_write(input_file: Path, output_file: Path, **kwargs) -> None:
    """Load a file that has named formating fields, e.g. {job_name}, format it, and write"""

    input_str = "".join(open(input_file, "r").readlines())
    output_str = input_str.format(**kwargs)
    with open(output_file, "w") as output:
        output.write(output_str)


@click.command(short_help="Create a submission script for Slurm.")
@click.option(
    "-p",
    "--pod5_dir",
    type=str,
    required=False,
    help="Directory containing POD5 files generated by MinKNOW (e.g. '<path>/<to>/pod5_pass'). Searched recursively.",
)
@click.option(
    "-k",
    "--barcode_kit",
    type=click.Choice(BARCODING_KITS),
    default="SQK-NBD114-96",
    show_default=True,
    help="Barcoding Kit.",
)
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'.",
)
@click.option(
    "-m",
    "--metadata_csv",
    type=str,
    required=True,
    help="Path to metadata CSV file containing barcode and sample information.",
)
@click.option(
    "-r",
    "--region_bed",
    type=str,
    required=True,
    help="Path to BED file specifying genomic regions of interest.",
)
@click.option(
    "-P",
    "--pipeline",
    type=click.Choice(PIPELINES),
    default=PIPELINES[0],
    show_default=True,
    help="Set the pipeline to be run.",
)
def main(
    pod5_dir: str,
    barcode_kit: str,
    expt_name: str,
    metadata_csv: str,
    region_bed: str,
    pipeline: str,
):

    # Imports here to reduce CLI latency
    from savanna.util.dirs import produce_dir, ExperimentDirectories
    from savanna.util.metadata import MetadataTableParser

    # Fixed
    include_unclassified = False
    basecall_model = "sup"

    print("Generating Slurm submission scripts...")
    print(f"  POD5 directory: {pod5_dir}")
    print(f"  Basecalling model: {basecall_model}")
    print(f"  Barcoding Kit: {barcode_kit}")
    print(f"  Experiment name: {expt_name}")
    print(f"  Metadata CSV: {metadata_csv}")
    print(f"  Regions BED: {region_bed}")
    print(f"  Pipeline: {pipeline}")
    print(f"  Include unclassified? {include_unclassified}")

    # Parse inputs
    _ = ExperimentDirectories(
        expt_name
    )  # builds the output  directories; can help with races
    metadata = MetadataTableParser(metadata_csv, include_unclassified)
    array_str = ",".join([str(int(b[-2:])) for b in metadata.barcodes])

    print(f"Found {metadata.n_barcodes} barcodes to process: {array_str}")

    # Prepare output directory
    run_count = sum([1 for d in os.listdir(RUNS_DIR) if d.startswith(expt_name)])
    output_dir = Path(produce_dir(RUNS_DIR, f"{expt_name}-r{run_count}"))

    load_format_write(
        slurm_basecall,
        output_dir / slurm_basecall.name,
        expt_name=expt_name,
        pod5_dir=pod5_dir,
        model=basecall_model,
        kit=barcode_kit,
        min_qscore=10,
    )

    load_format_write(
        slurm_analyse,
        output_dir / slurm_analyse.name,
        expt_name=expt_name,
        metadata=metadata_csv,
        regions=region_bed,
        pipeline=pipeline,
        array_str=array_str,
    )

    load_format_write(
        slurm_summary,
        output_dir / slurm_summary.name,
        expt_name=expt_name,
        metadata=metadata_csv,
        regions=region_bed,
        pipeline=pipeline,
    )

    submit_path = output_dir / slurm_submit.name
    load_format_write(slurm_submit, submit_path, run_dir=output_dir)
    submit_path.chmod(0o777)

    print(f"Final submission script written to: {submit_path}")
    print("Done.")


if __name__ == "__main__":
    main()
