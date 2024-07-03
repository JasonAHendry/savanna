import json
import pandas as pd

import matplotlib.pyplot as plt

from .interface import ExperimentAnalysis
from savanna.run.barcode.interface import BarcodeAnalysis
from savanna.run.barcode.bedcov import BarcodeBEDCoverage
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser


class ExperimentBedCoverage(ExperimentAnalysis):
    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        only_barcode: str = None,
        only_summary: bool = False,
        make_plot: bool = True,
    ):
        self.regions = regions
        super().__init__(expt_dirs, metadata, only_barcode, only_summary, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeBEDCoverage(barcode_name, self.expt_dirs, self.regions)

    def _summarise(self):
        """Combine into a single CSV"""

        bedcov_dfs = []
        for b, (success, outputs) in self.results.items():
            output_csv = outputs[0]
            bedcov_dfs.append(pd.read_csv(output_csv))
        self.bedcov_df = pd.concat(bedcov_dfs)

        df_path = f"{self.expt_dirs.approach_dir}/summary.bedcov.csv"
        self.bedcov_df.to_csv(df_path, index=False)

    def _plot(self):
        fig, ax = plt.subplots(1, 1, figsize=(4, 8))

        self.bedcov_df.index = self.bedcov_df["barcode"]
        self.bedcov_df[["mean_cov", "n_reads"]].plot(lw=0, marker="o", ax=ax)
        fig.savefig(f"{self.expt_dirs.approach_dir}/plot.bedcov.pdf")


# class ExperimentBedCoverage(ExperimentAnalysis):
#     def __init__(
#         self,
#         expt_dirs: ExperimentDirectories,
#         metadata: MetadataTableParser,
#         regions: RegionBEDParser,
#         make_plot: bool = True,
#     ):
#         self.regions = regions
#         super().__init__(expt_dirs, metadata, make_plot)

#     def _run(self):
#         self.outputs = []
#         for b in self.metadata.barcodes:
#             analysis = BarcodeBEDCoverage(b, self.expt_dirs, regions=self.regions)
#             success = analysis.run()
#             if success:
#                 self.outputs.append(analysis.output_csv)

#     def _summarise(self):
#         """Combine into a single CSV"""

#         bedcov_dfs = []
#         for output_csv in self.outputs:
#             bedcov_dfs.append(pd.read_csv(output_csv))
#         bedcov_df = pd.concat(bedcov_dfs)
#         df_path = f"{self.expt_dirs.approach_dir}/summary.bedcov.csv"
#         bedcov_df.to_csv(df_path, index=False)

#         # barcode_dir = self.expt_dirs.get_barcode_dir(b)
#         # try:
#         #     df = pd.read_csv(f"{barcode_dir}/bedcov/{b}.{self.ref_name}.bedcov.csv")
#         #     bedcov_dfs.append(df)
#         # except FileNotFoundError:
#         #     continue

#     def _plot(self):
#         fig, ax = plt.subplots(1, 1, figsize=(4, 8))


#         #fig, ax = plt.subplots(1, 1, figsize=(6, 6))

#         # Could merge with sample id
