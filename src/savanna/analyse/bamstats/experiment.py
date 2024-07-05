import json
import pandas as pd
import matplotlib.pyplot as plt

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories

from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from .barcode import BamStats


class ExperimentBamStats(ExperimentAnalysis):
    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        reference: Reference = PlasmodiumFalciparum3D7(),
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):
        self.reference = reference
        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BamStats(barcode_name, self.expt_dirs, self.reference)

    def _summarise(self):
        """
        Combine all of the JSON files which contain counts of the
        number of FASTQ files observed for each barcode
        """
        if self.results is None:
            raise RuntimeError("Did not find any results.")

        # Load
        bamstat_dts = []
        for b, (success, outputs) in self.results.items():
            if not success:
                # should be a 'FailedAnalysisError'
                # or could skip?
                raise ValueError(f"Analysis {self.__class__.__name__} failed for {b}!")

            output_json = outputs[0]
            bamstat_dts.append(json.load(open(output_json, "r")))

        # Concat and write
        self.bamstats_df = pd.DataFrame(bamstat_dts)
        df_path = (
            f"{self.expt_dirs.approach_dir}/summary.bamstats.{self.reference.name}.csv"
        )
        self.bamstats_df.to_csv(df_path, index=False)

    def _plot(self):
        self.bamstats_df.index = self.bamstats_df.barcode
        PLOT_STATS = ["n_primary", "n_secondary", "n_chimera", "n_unmapped"]

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

        self.bamstats_df[PLOT_STATS].plot(
            kind="bar", stacked=True, ec="black", lw=0.5, width=0.8, ax=ax
        )
        ax.set_xlabel("No. Alignments")
        ax.set_ylabel("")
        ax.legend(bbox_to_anchor=(1, 1))

        fig.savefig(
            f"{self.expt_dirs.approach_dir}/plot.bamstats.{self.reference.name}.pdf",
            dpi=300,
            pad_inchs=0.5,
            bbox_inches="tight",
        )
