import pandas as pd
from dataclasses import dataclass
from itertools import product
from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.analyse._interfaces.barcode import NoBarcodeAnalysisNeeded
from .genotypes import SampleGenotypes
from .drugs import ANTIMALARIALS


class ExperimentResistance(ExperimentAnalysis):

    name = "resistance"

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        """No barcode-level analysis needed"""
        return NoBarcodeAnalysisNeeded(barcode_name, self.expt_dirs)

    def _summarise(self):
        """Combine into a single CSV"""

        # Load and clean
        df = pd.read_csv(
            f"{self.expt_dirs.approach_dir}/summary/call/bcftools.filtered.annotated.tsv",
            sep="\t",
        )
        df.query("mut_type == 'missense'", inplace=True)
        df.insert(11, "gene", [amplicon.split("-")[0] for amplicon in df["amplicon"]])

        # Extract sample genotypes
        sample_genotypes = [
            SampleGenotypes.from_sample_dataframe(sample, sample_df)
            for sample, sample_df in df.groupby("sample")
        ]

        # Create an output resistance table
        @dataclass
        class TableRow:
            sample_name: str
            drug_name: str
            prediction: str

        results = [
            TableRow(
                sample.name,
                drug.name,
                drug.predict(sample.genotypes).value,
            )
            for sample, drug in product(sample_genotypes, ANTIMALARIALS)
        ]
        result_df = pd.DataFrame(results)
        result_df.to_csv(
            f"{self.expt_dirs.summary_dir}/{self.name}/bcftools.filtered.resistance_status.csv",
            index=False,
        )

    def _plot(self):
        pass
