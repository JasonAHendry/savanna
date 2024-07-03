from abc import ABC, abstractmethod

from savanna.util.dirs import ExperimentDirectories
from savanna.util.metadata import MetadataTableParser
from savanna.util.regions import RegionBEDParser

from savanna.basecall.basecalling import Dorado
from savanna.demux.demultiplexing import GuppyBarcoder
from savanna.run.experiment.fastq import ExperimentFastQCount
from savanna.run.experiment.map import (
    ExperimentMapToPfalciparum,
    ExperimentMapUnmappedToHSapiens,
)
from savanna.run.experiment.bedcov import ExperimentBedCoverage


class ExperimentPipeline(ABC):
    """
    Encapsulate an experiment pipeline

    Should I have some slurm functioality here?
    Or at least Slurm functionality could use this

    """

    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        reference_name: str = "Pf3D7",
        only_barcode: str = None,
        only_summary: bool = False,
        make_plot: bool = True,
    ):
        # Storage
        self.expt_dirs = expt_dirs
        self.metadata = metadata
        self.regions = regions
        self.ref_name = reference_name

        # Behaviour
        self.only_barcode = only_barcode
        self.only_summary = only_summary
        self.make_plot = make_plot

    @abstractmethod
    def run(self):
        pass


class PfNativeBarcodingPipeline(ExperimentPipeline):
    def __init__(
        self,
        fastq_dir: str,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        reference_name: str = "Pf3D7",
        only_barcode: str = None,
        only_summary: bool = False,
        make_plot: bool = True,
    ):
        super().__init__(
            expt_dirs,
            metadata,
            regions,
            reference_name,
            only_barcode,
            only_summary,
            make_plot,
        )
        self.fastq_dir = fastq_dir

    def run(self):
        """
        Run demultiplexing and analysis
        """

        # Need to add basecalling with dorado here

        # Demultiplex
        demuxer = GuppyBarcoder(self.fastq_dir, kit="rapid96")
        demuxer.run(
            output_dir=self.expt_dirs.demux_dir,
            use_gpu=False,
            dry_run=False,
            strict=True,
            trim_barcodes=True,
        )

        # These are common behavioural keyword arguments
        kwargs = {self.only_barcode, self.only_summary, self.make_plot}

        # Map to P.f.
        map_pf = ExperimentMapToPfalciparum(self.expt_dirs, self.metadata, **kwargs)
        map_pf.run()

        # Remap to H.s.
        map_hs = ExperimentMapUnmappedToHSapiens(
            self.expt_dirs, self.metadata, **kwargs
        )
        map_hs.run()

        # Calculate BED Coverage
        bedcov = ExperimentBedCoverage(
            self.expt_dirs, self.metadata, self.regions, **kwargs
        )
        bedcov.run()
