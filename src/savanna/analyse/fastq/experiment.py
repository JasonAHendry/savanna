import json
import pandas as pd
import matplotlib.pyplot as plt

from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from .barcode import BarcodeFASTQCounter


class ExperimentFASTQCount(ExperimentAnalysis):
    name = "fastq"
    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeFASTQCounter(barcode_name, self.expt_dirs)

    def _summarise(self):
        """
        Combine all of the JSON files which contain counts of the
        number of FASTQ files observed for each barcode

        """
        if self.results is None:
            raise RuntimeError("Did not find any results.")

        # Interate and load all JSON as Dict
        fastq_dts = []
        for b, (success, outputs) in self.results.items():
            if not success:
                # should be a 'FailedAnalysisError'
                raise ValueError(f"Analysis {self.__class__.__name__} failed for {b}!")

            output_json = outputs[0]
            fastq_dts.append(json.load(open(output_json, "r")))

        # Concat and write
        self.fastq_df = pd.DataFrame(fastq_dts)
        df_path = f"{self.summary_dir}/summary.fastq.csv"
        self.fastq_df.to_csv(df_path, index=False)

    def _plot(self):
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

        # Could merge with sample id
        self.fastq_df.index = self.fastq_df.barcode
        self.fastq_df["n_fastq"].plot(kind="bar", ax=ax)
        fig.savefig(
            f"{self.summary_dir}/plot.fastq.pdf",
            dpi=300,
            pad_inchs=0.5,
            bbox_inches="tight",
        )