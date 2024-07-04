from abc import ABC, abstractmethod

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import (
    PlasmodiumFalciparum3D7,
    AnophelesGambiaePEST,
)

from savanna.analyse.fastq.experiment import ExperimentFASTQCount
from savanna.analyse.map.experiment import (
    ExperimentMapToReference,
)
from savanna.analyse.bedcov.experiment import ExperimentBedCoverage


class Pipeline(ABC):
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


class PlasmoPipeline(Pipeline):
    """
    Basic pipeline for plasmodium falciparum

    """

    reference = PlasmodiumFalciparum3D7()

    def run(self):
        print(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()

        print("FASTQ...")
        analyse_fastq = ExperimentFASTQCount(
            self.expt_dirs, self.metadata, **self.kwargs
        )
        analyse_fastq.run()

        print("Map...")
        mapper = ExperimentMapToReference(
            self.expt_dirs, self.metadata, self.reference, **self.kwargs
        )
        mapper.run()

        print("Coverage...")
        bedcov = ExperimentBedCoverage(
            self.expt_dirs,
            self.metadata,
            self.regions,
            self.reference,
            **self.kwargs,
        )
        bedcov.run()


class EntoPipeline(Pipeline):
    """
    Basic pipeline for vectors falciparum

    """

    reference = AnophelesGambiaePEST()

    def run(self):
        print(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()

        print("FASTQ...")
        analyse_fastq = ExperimentFASTQCount(self.expt_dirs, self.metadata, **self.kwargs)
        analyse_fastq.run()

        print("Map...")
        mapper = ExperimentMapToReference(self.expt_dirs, self.metadata, self.reference, **self.kwargs)
        mapper.run()

        print("Coverage...")
        bedcov = ExperimentBedCoverage(
            self.expt_dirs, self.metadata, self.regions, self.reference, **self.kwargs
        )
        bedcov.run()


PIPELINE_COLLECTION = {
    "plasmo": PlasmoPipeline,
    "ento": EntoPipeline,
}
