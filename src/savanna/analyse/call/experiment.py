from typing import Dict

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from .callers import VariantCaller
from .barcode import BarcodeVariantCaller


class ExperimentVariantCaller(ExperimentAnalysis):
    name = "call"

    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        caller: VariantCaller,
        reference: Reference = PlasmodiumFalciparum3D7(),
        make_plot: bool = True,
        barcode: str = None,
        summary_only: bool = False,
    ):
        self.reference = reference
        self.caller = caller
        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeVariantCaller(
            barcode_name=barcode_name,
            expt_dirs=self.expt_dirs,
            reference=self.reference,
            caller=self.caller,
        )

    def _summarise(self):
        """Merge the VCF files and convert to a TSV"""
        pass

    def _plot(self):
        pass
