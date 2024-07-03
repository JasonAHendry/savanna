#import os
import click


# Duplicated because slows API considerably
BARCODING_KIT_MAPPING = {"native96": "SQK-NBD114-96", "rapid96": "SQK-RBK114-96"}

@click.command(short_help="Run demultiplexing for an experiment.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, e.g. '2024-01-01_pilot'. Used to create output folder.",
)
@click.option(
    "-k",
    "--barcode_kit",
    type=click.Choice(BARCODING_KIT_MAPPING),
    default="native96",
    show_default=True,
    help="Barcoding Kit.",
)
@click.option(
    "-i",
    "--input_fastq_dir",
    type=str,
    default=None,
    help="Input FASTQ directory for demultiplexing. If omitted, will assume location based on experiment name.",
)
@click.option(
    "--strict",
    type=bool,
    default=True,
    is_flag=True,
    help="Increase strictness of barcode allocations.",
)  # Can add other optionals here
def demux(expt_name: str, barcode_kit: str, input_fastq_dir: str, strict: bool) -> None:
    """
    Download target reference genome by specifying a `reference_name`, or download
    all available genomes with the `--all` flag

    """

    from savanna.util.dirs import ExperimentDirectories
    from .demultiplexing import GuppyBarcoder

    expt_dirs = ExperimentDirectories(expt_name)

    # Set the input directory, if it doesn't exist
    if input_fastq_dir is None:
        input_fastq_dir = expt_dirs.basecall_dir
        # Improvement is to check for a FASTQ file...
        n_fastqs = sum(
            [
                1
                for f in os.listdir(input_fastq_dir)
                if f.endswith(".fastq") or f.endswith(".fastq.gz")
            ]
        )
        if n_fastqs == 0:
            raise FileNotFoundError(
                f"No input FASTQ files found in {input_fastq_dir}. Please explicitly set `-i` flag."
            )

    # Define output directory
    output_dir = expt_dirs.demux_dir

    print("Building basecaller...")
    demuxer = GuppyBarcoder(fastq_dir=input_fastq_dir, kit=barcode_kit)
    demuxer.run(output_dir=output_dir, strict=strict)
