import click
from savanna.gather.commands import gather
from savanna.basecall.commands import basecall
from savanna.demux.commands import demux
from savanna.run.commands import run


@click.group()
def cli():
    """
    Analyse targeted nanopore sequencing data for genomic
    surveillance

    """
    pass


cli.add_command(gather)
cli.add_command(basecall)
cli.add_command(demux)
cli.add_command(run)
