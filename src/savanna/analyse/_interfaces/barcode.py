import os
from typing import List
from abc import ABC, abstractmethod
from savanna.util.dirs import ExperimentDirectories, produce_dir
from savanna.util.exceptions import InputDoesNotExit


class BarcodeAnalysis(ABC):
    """
    Run an analysis for a single barcode

    """

    name = ""

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        make_plot: bool = True,
    ):

        # Define the specific barcode
        self.barcode_name = barcode_name
        self.barcode_dir = expt_dirs.get_barcode_dir(barcode_name)

        # Define output directory
        self.output_dir = produce_dir(self.barcode_dir, self.name)

        self.expt_dir = expt_dirs
        self.make_plot = make_plot

        # This must be defined
        self.inputs = self._define_inputs()
        self.outputs = self._define_outputs()

    @abstractmethod
    def _define_inputs(self) -> List[str]:
        """
        Define the expected input files in a list

        """
        pass

    @abstractmethod
    def _define_outputs(self) -> List[str]:
        """
        Define the generated output files in a list;
        this helps downstream methods to collect the output files
        without needing specific knowledge of the output file paths

        """
        pass

    @property
    def inputs_exist(self):
        return all([os.path.exists(input) for input in self.inputs])

    @property
    def outputs_exist(self):
        return all([os.path.exists(output) for output in self.outputs])

    @abstractmethod
    def _run(self) -> None:
        """
        Run analysis on the input files

        """
        pass

    @abstractmethod
    def _plot(self) -> None:
        """
        Plot results of the analysis

        """
        pass

    def _check_inputs(self) -> None:
        """
        Check that all the necessary input files exist for
        a given barcode analysis

        """
        if self.inputs is None:
            raise RuntimeError("`self.inputs` must be defined.")

        if not self.inputs_exist:
            for (
                input
            ) in self.inputs:  # verbose but important to have good feedback here
                if not os.path.exists(input):
                    raise InputDoesNotExit(
                        f"Input '{input}' does not exist, but is needed for {self.__class__.__name__} analysis."
                    )

    def run(self) -> bool:
        self._check_inputs()
        self._run()
        if self.make_plot:
            self._plot()
        return self.outputs_exist


class NoBarcodeAnalysisNeeded(BarcodeAnalysis):
    """
    This class can be used when there is only experiment-level analysisi
    performed as part of
    """

    def _define_inputs(self) -> List[str]:
        """No inputs"""
        return [self.barcode_dir]

    def _define_outputs(self) -> List[str]:
        return [self.barcode_dir]

    def _run(self):
        """Nothing is computed"""
        return True

    def _plot(self):
        pass
