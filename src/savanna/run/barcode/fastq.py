import os
import json
from typing import List

from savanna.run.barcode._interface import BarcodeAnalysis


class BarcodeFASTQCounter(BarcodeAnalysis):
    """
    Count the number or FASTQ files generated for a particular barcode

    """

    name = "fastqs"

    def _define_inputs(self) -> List[str]:
        """
        Check that the expected input directory is present

        """

        self.fastq_dir = f"{self.expt_dir.demux_dir}/{self.barcode_name}"

        return [self.fastq_dir]

    def _define_outputs(self) -> List[str]:

        self.output_json = f"{self.output_dir}/{self.barcode_name}.fastq.json"

        return [self.output_json]

    @staticmethod
    def _is_fastq(input_file: str) -> bool:
        return input_file.endswith(".fastq") or input_file.endswith(".fastq.gz")

    def _run(self) -> None:
        """
        Count the number of fastq_files

        """

        fastq_files = [f for f in os.listdir(self.fastq_dir) if self._is_fastq(f)]
        n_fastq = len(fastq_files)

        json.dump(
            {
                "barcode": self.barcode_name,
                "fastq_dir": self.fastq_dir,
                "n_fastq": n_fastq,
            },
            open(self.output_json, "w"),
        )

    def _plot(self):
        pass
