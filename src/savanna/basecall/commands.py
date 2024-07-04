import click


# Redefined here because import slows the CLI
BASECALL_MODELS = ["fast", "hac", "sup"]


# TODO:
# - In future, model download should be automatically part of the CLI


@click.command(short_help="Run basecalling for an experiment.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=False,
    help="Name of the experiment, e.g. '2024-01-01_pilot'. Used to create output folder.",
)
@click.option(
    "-p",
    "--pod5_dir",
    type=str,
    required=False,
    help="Parent directory inside which POD5 files exist. This directory is searched recursively.",
)
@click.option(
    "-m",
    "--model",
    type=click.Choice(BASECALL_MODELS),
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
@click.option(
    "--download_models",
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help="Download all basecalling models. Must be run before using this command: `savanna basecall --download_models`.",
)
def basecall(
    expt_name: str, pod5_dir: str, model: str, min_qscore: int, download_models: bool
) -> None:
    """
    Run basecalling for an experiment using `dorado basecaller`

    """

    from .basecalling import Dorado
    from savanna.util.dirs import ExperimentDirectories

    if download_models:
        basecaller = Dorado()
        for model in BASECALL_MODELS:
            basecaller.download_model(model)
        return

    if (expt_name is None) or (pod5_dir is None):
        raise click.UsageError(
            "If not using '--download_models', must provide options '-e' / '--expt_name' and '-p' / '--pod5_dir'."
        )

    print(f"Loaded experiment: {expt_name}")
    expt_dirs = ExperimentDirectories(expt_name)
    output_fastq = f"{expt_dirs.basecall_dir}/dorado.fastq.gz"
    print(f"Output FASTQ: {output_fastq}")

    print("Runninng basecaller...")
    basecaller = Dorado()
    basecaller.run(
        model=model,
        pod5_dir=pod5_dir,
        output_fastq=output_fastq,
        min_qscore=min_qscore,
    )
    print("Done.")
