import click
from collections import OrderedDict
from savanna.download.command import download
from savanna.basecall.command import basecall
from savanna.demultiplex.command import demultiplex
from savanna.analyse.command import analyse


# From: https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
class OrderedGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super(OrderedGroup, self).__init__(name, commands, **attrs)
        #: the registered subcommands by their exported names.
        self.commands = commands or OrderedDict()

    def list_commands(self, ctx):
        return self.commands


@click.group(cls=OrderedGroup)
@click.version_option(message="%(prog)s-v%(version)s")
def cli():
    """
    Analyse targeted nanopore sequencing data for genomic
    surveillance

    """
    pass


cli.add_command(download)
cli.add_command(basecall)
cli.add_command(demultiplex)
cli.add_command(analyse)


# This will be for the specific steps
# @click.group()
# def only():
#     """ Only run a specific command"""
#     pass


# cli.add_command(only)
# only.add_command(map)