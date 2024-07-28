import os
from typing import List
from abc import ABC, abstractmethod
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories, produce_dir
from savanna.util.exceptions import InputDoesNotExit


class SummaryAnalysis(ABC):
    """
    Perform an analysis on experiment summary files

    NB:
    - Bad code duplication with BarcodeAnalysis
    - But for now just keeping, will consider reorganising
      these abstract classes later

    """
    
    name = ""

    def __init__(
            self,
            expt_dirs: ExperimentDirectories,
            metadata: MetadataTableParser,
            make_plot: bool = True,
        ):


            self.expt_dirs = expt_dirs
            self.metadata = metadata
            self.make_plot = make_plot

            # Define output directory
            self.output_dir = produce_dir(self.expt_dirs.summary_dir, self.name)

            # This must be defined in child class
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

