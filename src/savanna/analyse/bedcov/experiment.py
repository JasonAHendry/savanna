import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser

from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from .barcode import BarcodeBEDCoverage


class ExperimentBedCoverage(ExperimentAnalysis):

    name = "bedcov"
    
    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        reference: Reference = PlasmodiumFalciparum3D7(),
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):
        self.regions = regions
        self.reference = reference
        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeBEDCoverage(
            barcode_name, self.expt_dirs, self.regions, self.reference
        )

    def _summarise(self):
        """Combine into a single CSV"""

        bedcov_dfs = []
        for b, (success, outputs) in self.results.items():
            output_csv = outputs[0]
            bedcov_dfs.append(pd.read_csv(output_csv))
        self.bedcov_df = pd.concat(bedcov_dfs)

        df_path = f"{self.summary_dir}/summary.bedcov.csv"
        self.bedcov_df.to_csv(df_path, index=False)

    def _plot(self):
        fig, ax = plt.subplots(1, 1, figsize=(4, 6))

        sns.stripplot(
            y="barcode", x="mean_cov", hue="name", jitter=True, data=self.bedcov_df
        )
        ax.legend(bbox_to_anchor=(1, 1))

        fig.savefig(
            f"{self.summary_dir}/plot.bedcov.pdf",
            bbox_inches="tight",
            pad_inches=0.5,
        )

