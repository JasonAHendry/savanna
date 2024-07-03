import json
import pandas as pd

import matplotlib.pyplot as plt

from savanna.run.barcode._interface import BarcodeAnalysis

from ._interface import ExperimentAnalysis
from savanna.run.barcode.fastq import BarcodeFASTQCounter

# Update to work with new interface
class ExperimentFASTQCount(ExperimentAnalysis):
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
        df_path = f"{self.expt_dirs.approach_dir}/summary.fastq.csv"
        self.fastq_df.to_csv(df_path, index=False)

    def _plot(self):
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

        # Could merge with sample id
        self.fastq_df.index = self.fastq_df.barcode
        self.fastq_df["n_fastq"].plot(kind="bar", ax=ax)
        fig.savefig(
            f"{self.expt_dirs.approach_dir}/plot.fastq.pdf",
            dpi=300,
            pad_inchs=0.5,
            bbox_inches="tight",
        )


# # Update to work with new interface
# class ExperimentFASTQCount(ExperimentAnalysis):
#     def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
#         return BarcodeFASTQCounter(barcode_name, self.expt_dirs)

#     def _get_output_json(self, barcode_name: str):
#         """ This is coupled to BarcodeFASTQCounter """
#         barcode_dir = self.expt_dirs.get_barcode_dir(barcode_name)
#         analysis_dir = f"{barcode_dir}/{BarcodeFASTQCounter.name}"
#         return f"{analysis_dir}/{barcode_name}.fastq.json"

#     def _summarise(self):
#         """
#         Combine all of the JSON files which contain counts of the
#         number of FASTQ files observed for each barcode

#         """

#         # Interate and load all JSON as Dict
#         fastq_dts = []
#         for s, b in zip(self.success, self.metadata.barcodes):
#             if not s:
#                 # should be a 'FailedAnalysisError'
#                 raise ValueError(f"Analysis {self.__class__.__name__} failed for {b}!")
#             output_json = self._get_output_json(b)
#             fastq_dts.append(json.load(open(output_json, "r")))

#         # Concat and write
#         self.fastq_df = pd.DataFrame(fastq_dts)
#         df_path = f"{self.expt_dirs.approach_dir}/summary.fastq.csv"
#         self.fastq_df.to_csv(df_path, index=False)

#     def _plot(self):
#         fig, ax = plt.subplots(1, 1, figsize=(6, 6))

#         # Could merge with sample id
#         self.fastq_df.index = self.fastq_df.barcode
#         self.fastq_df["n_fastq"].plot(kind="bar", ax=ax)
#         fig.savefig(
#             f"{self.expt_dirs.approach_dir}/plot.fastq.pdf",
#             dpi=300,
#             pad_inchs=0.5,
#             bbox_inches="tight",
#         )


# class ExperimentFastQCount(ExperimentAnalysis):
#     def _run(self):
#         self.outputs = []
#         for b in self.metadata.barcodes:
#             analysis = BarcodeFASTQCounter(b, self.expt_dirs)
#             success = analysis.run()
#             if success:
#                 self.outputs.append(analysis.output_json)

#     def _summarise(self):
#         """Combine into a single CSV"""
#         fastq_dts = []
#         for output in self.outputs:
#             try:
#                 dt = json.load(open(output, "r"))
#                 fastq_dts.append(dt)
#             except FileNotFoundError:
#                 continue
#         df = pd.DataFrame(fastq_dts)
#         df_path = f"{self.expt_dirs.approach_dir}/summary.fastq.csv"
#         df.to_csv(df_path, index=False)

#         self.fastq_df = df
