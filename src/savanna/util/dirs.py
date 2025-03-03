import os
import shutil
import time
from pathlib import Path

from .metadata import MetadataTableParser
from .regions import RegionBEDParser


ROOT_DIR = Path(__file__).absolute().parent.parent.parent.parent


def produce_dir(*args):
    """
    Produce a new directory by concatenating `args`,
    if it does not already exist

    params
        *args: str1, str2, str3 ...
            Comma-separated strings which will
            be combined to produce the directory,
            e.g. str1/str2/str3

    returns
        dir_name: str
            Directory name created from *args.

    """

    # Define directory path
    dir_name = os.path.join(*args)

    # Create directory; handle possible file-system latency
    n_tries = 5
    for _ in range(n_tries):
        try:
            os.makedirs(dir_name, exist_ok=True)
            return dir_name
        except FileExistsError:
            time.sleep(0.1)
    raise FileExistsError(f"Failed to create directory {dir_name} despite {n_tries} attempts.")


class ExperimentDirectories:
    """
    Put all the information about experimental
    directory structure in a single place

    Advantage is any future changes to how the directories
    are organised can be managed here

    Best practice is that objects/functions *do not* directly depend on
    this object, better to use its members as arguments to other
    objects or functions

    """

    results_dir = produce_dir(ROOT_DIR, "results")

    def __init__(
        self,
        expt_name: str,
        metadata: MetadataTableParser = None,
        regions: RegionBEDParser = None,
        fastq_dir: str = None,
        approach_name: str = "",
    ):
        """
        Create directory heirarchy for an experiment

        """

        self.expt_name = expt_name
        self.expt_dir = produce_dir(self.results_dir, expt_name)

        # This enables different Guppy versions, barcoding strategies, &c
        self.approach_name = approach_name
        self.approach_dir = produce_dir(self.expt_dir, approach_name)

        # Output for basecalling & demultiplexing
        self.basecall_dir = produce_dir(self.approach_dir, "basecalling")
        if fastq_dir is not None:
            self.demux_dir = fastq_dir
        else:
            self.demux_dir = produce_dir(self.approach_dir, "demultiplexing")

        # Setup metadata
        self.metadata_dir = produce_dir(self.expt_dir, "metadata")
        self._populate_metadata_dir(metadata, regions)

        # For logs
        self.log_dir = produce_dir(self.expt_dir, "logs")

        # This is for experimental summaries
        self.summary_dir = produce_dir(self.approach_dir, "summary")

        # Output for individual barcode analyses
        self.barcodes_dir = produce_dir(self.approach_dir, "barcodes")
        self._barcode_dirs = self._create_barcode_dirs_dict(metadata)

    def get_barcode_dir(self, barcode_name: str):
        return self._barcode_dirs[barcode_name]

    def _populate_metadata_dir(
        self, metadata: MetadataTableParser, regions: RegionBEDParser
    ) -> None:
        """
        Move metadata CSV and regions BED into the metadata directory,
        and store their paths as attributes
        """
        if metadata is not None:
            self.metadata_csv = f"{self.metadata_dir}/{os.path.basename(metadata.csv)}"
            if not os.path.exists(self.metadata_csv):
                metadata.df.to_csv(self.metadata_csv, index=False)

        if regions is not None:
            self.regions_bed = f"{self.metadata_dir}/{os.path.basename(regions.path)}"
            if not os.path.exists(self.regions_bed):
                shutil.copy(regions.path, self.regions_bed)

    def _create_barcode_dirs_dict(self, metadata: MetadataTableParser) -> None:
        """
        This can lead to race conditios on a cluster;
        do I need it?
        """
        if metadata is None:
            return None
        return {b: produce_dir(self.barcodes_dir, b) for b in metadata.barcodes}
