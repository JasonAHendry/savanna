from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.metadata import MetadataTableParser
from savanna.download.references import Reference
from .barcode import (
    BarcodeMapToReference,
    BarcodeMapToPfalciparum,
    BarcodeMapUnmappedToHSapiens,
)


class ExperimentMapToReference(ExperimentAnalysis):
    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        reference: Reference,
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):
        self.reference = reference
        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapToReference(barcode_name, self.expt_dirs, self.reference)

    def _summarise(self):
        pass

    def _plot(self):
        pass


class ExperimentMapToPfalciparum(ExperimentAnalysis):
    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapToPfalciparum(barcode_name, self.expt_dirs)

    def _summarise(self):
        pass

    def _plot(self):
        pass


class ExperimentMapUnmappedToHSapiens(ExperimentAnalysis):
    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapUnmappedToHSapiens(barcode_name, self.expt_dirs)

    def _summarise(self):
        pass

    def _plot(self):
        pass
