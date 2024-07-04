import os
import click


# TODO:
# - Would be great to include some quality control here
# -> using either barcoding_summary.txt
# -> or scanning directories
# - Former is probably better


DEMUX_TOOLS = ["dorado", "guppy"]
BARCODING_KITS = ["SQK-NBD114-96", "SQK-RBK114-96"]


@click.command(short_help="FASTQ to per-sample FASTQ.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, e.g. '2024-01-01_pilot'. Used to create output folder.",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=str,
    default=None,
    required=False,
    help="Directory containing FASTQ file(s) for demultiplexing. May be from MinKNOW, savanna basecall, or other source. When omitted, location is assumed from '--expt_name'."
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
    "-t",
    "--tool",
    type=click.Choice(DEMUX_TOOLS),
    default="dorado",
    show_default=True,
    help="Tool to use for demultiplexing. Note that dorado does not support '--strict' at present."
)
@click.option(
    "--strict",
    type=bool,
    default=True,
    is_flag=True,
    help="Increase strictness of barcode allocations.",
)
def demultiplex(expt_name: str, fastq_dir: str, barcode_kit: str, tool: str, strict: bool) -> None:
    """
    Performing demultiplexing: split FASTQ data based on sample barcodes detected
    in reads

    """

    from savanna.util.dirs import ExperimentDirectories
    from .demultiplexing import DEMUXER_COLLECTION

    expt_dirs = ExperimentDirectories(expt_name)

    # Assume FASTQ directory location, if not provided
    if fastq_dir is None:
        fastq_dir = expt_dirs.basecall_dir
        n_fastqs = sum(
            [
                1
                for f in os.listdir(fastq_dir)
                if f.endswith(".fastq") or f.endswith(".fastq.gz")
            ]
        )
        if n_fastqs == 0:
            raise FileNotFoundError(
                f"No input FASTQ files found in {fastq_dir}. Please explicitly set '-f' / '--fastq_dir' flag."
            )
        
    print(fastq_dir)

    # Define output directory
    output_dir = expt_dirs.demux_dir

    print("Building basecaller...")
    demuxer = DEMUXER_COLLECTION[tool](fastq_dir=fastq_dir, kit=barcode_kit)
    print(demuxer.fastq_dir)
    demuxer.run(output_dir=output_dir, strict=strict, dry_run=True)


