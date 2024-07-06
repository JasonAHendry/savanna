from typing import List
from abc import ABC, abstractmethod

from .barcode import BarcodeAnalysis
from .result import BarcodeAnalysisResults
from savanna.util.metadata import MetadataTableParser, check_barcode_format
from savanna.util.dirs import produce_dir, ExperimentDirectories


class ExperimentAnalysis(ABC):
    """
    Perform a given per-barcode analysis across an entire experiment (e.g.
    all of the barcodes in that exeperiment); and the optionally
    summarise the results and produce a plot
    """

    name = ""

    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):

        # Storage
        self.expt_dirs = expt_dirs
        self.metadata = metadata

        # Behaviour modifiers
        self.barcode = self._set_barcode(barcode)
        self.summary_only = summary_only
        self.make_plot = make_plot

        # Object pointing to expected results
        self.results = None

        # Outputs of _summarise() and _plot() go here
        self.summary_dir = produce_dir(self.expt_dirs.summary_dir, self.name)

    def _set_barcode(self, barcode: str):
        """
        Check that the `barcode` argument conforms to a valid
        barcode within the metadata

        """
        if barcode in [None, "unclassified"]:
            return barcode

        barcode = check_barcode_format(barcode, try_to_fix=True)
        if not barcode in self.metadata.barcodes:
            raise ValueError(
                f"Barcode {barcode} cannot be analysed, it is not in {self.metadata.csv}"
            )
        return barcode

    @abstractmethod
    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        pass

    def _collect_outputs(self, barcodes: List[str]) -> BarcodeAnalysisResults:
        results = BarcodeAnalysisResults(barcodes)
        for b in barcodes:
            barcode_analysis = self._get_barcode_analysis(b)
            results[b] = (barcode_analysis.outputs_exist, barcode_analysis.outputs)
        return results

    def _run(self, barcodes: List[str]) -> BarcodeAnalysisResults:
        """
        Run a specific `BarcodeAnalysis` for a list of barcodes

        Note on usage:
        - In some subclasses, this will get overridden with 'pass' or simply
        a check for the correct input files.
        - Ultimately, self._summarise() will often expect a success = List[bool],
        as it will only perform summaries for barcodes where the analysis has run
        successfully or the expected input files were found; otherwise it
        should log a warning.

        """

        results = BarcodeAnalysisResults(barcodes)
        for b in barcodes:
            barcode_analysis = self._get_barcode_analysis(b)
            s = barcode_analysis.run()
            results[b] = (s, barcode_analysis.outputs)
        return results

    @abstractmethod
    def _summarise(self):
        pass

    @abstractmethod
    def _plot(self):
        pass

    def run(self):
        """
        Run the experiment analysis, optionally only running for a single
        barcode, performing just the summary, or not plotting

        """

        if self.summary_only:
            self.results = self._collect_outputs(self.metadata.barcodes)
        else:
            if self.barcode is not None:
                self._run([self.barcode])
                return # we don't summarise if only one barcode was run.
            self.results = self._run(self.metadata.barcodes)

        self._summarise()
        if self.make_plot:
            self._plot()
