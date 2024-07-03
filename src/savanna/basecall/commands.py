import click


MODEL_CHOICES = ["fast", "hac", "sup"]


@click.command(short_help="Run basecalling for an experiment.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, e.g. '2024-01-01_pilot'. Used to create output folder.",
)
@click.option(
    "-p",
    "--pod5_dir",
    type=str,
    required=True,
    help="Parent directory inside which POD5 files exist. This directory is searched recursively.",
)
@click.option(
    "-m",
    "--model",
    type=click.Choice(MODEL_CHOICES),
    default="sup",
    show_default=True,
    help="Basecalling model.",
)
@click.option(
    "-q",
    "--min_qscore",
    type=int,
    default=10,
    show_default=True,
    help="Minimum quality score threshold.",
)
def basecall(expt_name: str, pod5_dir: str, model: str, min_qscore: int) -> None:
    """
    Run basecalling for an experiment using `dorado basecaller`

    """

    from .basecalling import Dorado
    from savanna.util.dirs import ExperimentDirectories

    print(f"Loaded experiment: {expt_name}")
    expt_dirs = ExperimentDirectories(expt_name)
    output_fastq = f"{expt_dirs.basecall_dir}/dorado.fastq.gz"
    print(f"Output FASTQ: {output_fastq}")

    print("Building basecaller...")
    basecaller = Dorado(pod5_dir=pod5_dir, model=model)
    print("Running...")
    basecaller.run(
        output_fastq=output_fastq,
        min_qscore=min_qscore,
    )

