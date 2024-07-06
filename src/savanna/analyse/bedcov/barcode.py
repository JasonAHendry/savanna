import pandas as pd
from typing import List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from savanna.wrappers import samtools


class BarcodeBEDCoverage(BarcodeAnalysis):
    """
    Compute coverage summary statistics from a BAM file over
    regions defined in a BED file

    """

    name = "bedcov"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        reference: Reference = PlasmodiumFalciparum3D7(),
        make_plot: bool = True,
    ):
        self.regions = regions
        self.reference = reference
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self):
        self.bam_path = (
            f"{self.barcode_dir}/bams/{self.barcode_name}.{self.reference.name}.filtered.bam"
        )
        return [self.bam_path, self.regions.path]

    def _define_outputs(self) -> List[str]:
        self.output_csv = f"{self.output_dir}/table.bedcov.csv"  # may want more here
        return [self.output_csv]

    def _run(self):
        samtools.bedcov(
            bam_path=self.bam_path,
            bed_path=self.regions.path,
            output_csv=self.output_csv,
        )

        # Load the dataframe to add a barcode column indicator
        df = pd.read_csv(self.output_csv)
        df.insert(0, "barcode", self.barcode_name)
        df.to_csv(self.output_csv, index=False)

    def _plot(self):
        pass
