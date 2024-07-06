from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.metadata import MetadataTableParser
from savanna.download.references import Reference
from .barcode import FilterBAM, MIN_MAPQ, EXCLUDE_FLAGS


class ExperimentFilterBAM(ExperimentAnalysis):
    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        reference: Reference,
        min_mapq: int = MIN_MAPQ,
        exclude_flags: str = EXCLUDE_FLAGS,
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):
        self.reference = reference
        self.min_mapq = min_mapq
        self.exclude_flags = exclude_flags
        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return FilterBAM(
            barcode_name, 
            self.expt_dirs, 
            self.reference,
            min_mapq=self.min_mapq,
            exclude_flags=self.exclude_flags
        )

    def _summarise(self):
        pass

    def _plot(self):
        pass
