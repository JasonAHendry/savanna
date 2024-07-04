import click


PIPELINES = ["plasmo", "ento"]


@click.command(short_help="Per-sample FASTQ to results.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'.",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=str,
    required=True,
    help="Path to directory containing demultiplexed FASTQ files (e.g. '<path>/<to>/fastq_pass'). Typically produced by MinKNOW, dorado, guppy or with savanna demultiplex.",
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
    "-p",
    "--pipeline",
    type=click.Choice(PIPELINES),
    default=PIPELINES[0],
    show_default=True,
    help="Set the pipeline to be run.",
)
@click.option(
    "-b",
    "--barcode",
    type=int,
    required=False,
    help="Optionally run the pipeline for only one barcode (e.g. --barcode 3).",
)
@click.option(
    "-s",
    "--summary_only",
    type=bool,
    is_flag=True,
    default=False,
    help="Optionally only summarise existing results across barcodes.",
)
def analyse(
    expt_name, fastq_dir, metadata_csv, region_bed, pipeline, barcode, summary_only
):

    from .main import main

    main(
        expt_name, fastq_dir, metadata_csv, region_bed, pipeline, barcode, summary_only
    )
