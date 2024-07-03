from savanna.run.experiment._interface import ExperimentAnalysis
from savanna.run.barcode._interface import BarcodeAnalysis
from savanna.run.barcode.map import (
    BarcodeMapToPfalciparum,
    BarcodeMapUnmappedToHSapiens,
)


class ExperimentMapToPfalciparum(ExperimentAnalysis):
    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapToPfalciparum(barcode_name, self.expt_dirs)

    def _summarise(self):
        pass

    def _plot(self):
        pass


class ExperimentMapUnmappedToHSapiens(ExperimentAnalysis):
    """
    In this case, we really are just wrapping the mapping steps
    We don't produce any summary files or plots of the BAM files themselves

    """

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapUnmappedToHSapiens(barcode_name, self.expt_dirs)

    def _summarise(self):
        pass

    def _plot(self):
        pass

