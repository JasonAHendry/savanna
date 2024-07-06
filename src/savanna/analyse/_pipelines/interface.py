import logging
import yaml
from abc import ABC, abstractmethod

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference
from savanna.analyse.fastq.experiment import ExperimentFASTQCount
from savanna.analyse.map.experiment import ExperimentMapToReference
from savanna.analyse.bedcov.experiment import ExperimentBedCoverage
from savanna.analyse.bamstats.experiment import ExperimentBamStats
from savanna.analyse.bamfilt.experiment import ExperimentFilterBAM
from savanna.analyse.call.experiment import ExperimentCallWithBcftools


log = logging.getLogger()


class Pipeline(ABC):
    """
    Intereface for a pipline

    Notes:
    * Must define self.run()
    * Some common attributes assigned (self.expt_dirs)
    * Some common steps are packaged already (self.map_to_reference())

    """

    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        barcode: int = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):
        # Storage
        self.expt_dirs = expt_dirs
        self.metadata = metadata
        self.regions = regions
        self.params = None

        # Behaviour
        self.barcode = barcode
        self.summary_only = summary_only
        self.make_plot = make_plot
        self.kwargs = {
            "barcode": barcode,
            "summary_only": summary_only,
            "make_plot": make_plot,
        }

    @abstractmethod
    def run(self) -> None:
        pass

    def _load_parameters(self, parameter_path: str):
        with open(parameter_path, "r") as file:
            self.params = yaml.safe_load(file)

    def _count_fastqs(self) -> None:
        log.info("Counting FASTQ files")
        analyse_fastq = ExperimentFASTQCount(
            self.expt_dirs, self.metadata, **self.kwargs
        )
        analyse_fastq.run()
        log.info("Done.")

    def _map_to_reference(self, reference: Reference) -> None:
        log.info("Mapping reads")
        log.info(f" to: {reference.name}")
        mapper = ExperimentMapToReference(
            self.expt_dirs, self.metadata, reference, **self.kwargs
        )
        mapper.run()
        log.info("Done.")

    def _filter_bam(self, reference: Reference) -> None:
        log.info("Filtering BAM file")
        log.info(f"  Reference: {reference.name}")
        filter = ExperimentFilterBAM(
            self.expt_dirs, self.metadata, self.reference, **self.kwargs
        ) # so they would get passed here
        filter.run()
        log.info("Done.")

    def _calc_bamstat(self, reference: Reference) -> None:
        log.info("Calculating mapping statistics")
        log.info(f" to: {reference.name}")
        bamstats = ExperimentBamStats(
            self.expt_dirs, self.metadata, reference, **self.kwargs
        )
        bamstats.run()
        log.info("Done.")

    def _calc_bedcov(self, reference: Reference) -> None:
        log.info("Calculating coverage over amplicons")
        log.info(f"  Reference: {reference.name}")
        log.info(f"  Amplicons (BED): {self.regions.path}")
        bedcov = ExperimentBedCoverage(
            self.expt_dirs, self.metadata, self.regions, reference, **self.kwargs
        )
        bedcov.run()

    def _call_with_bcftools(self, reference: Reference) -> None:
        log.info("Calling variants with bcftools")
        log.info(f"  Reference: {reference.name}")
        log.info(f"  Amplicons (BED): {self.regions.path}")
        call = ExperimentCallWithBcftools(
            self.expt_dirs,
            self.metadata,
            self.regions,
            reference,
            **self.kwargs,
        )
        call.run()
        log.info("Done.")
