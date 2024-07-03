import pandas as pd
from typing import List
from savanna.run.barcode._interface import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.wrappers import samtools


class BarcodeBEDCoverage(BarcodeAnalysis):
    """
    Count the number or FASTQ files generated for a particular barcode

    """

    name = "bedcov"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        make_plot: bool = True,
    ):
        self.regions = regions
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self):
        """
        Check that the expected input directory is present

        """

        self.bam_path = f"{self.barcode_dir}/bams/{self.barcode_name}.Pf3D7.bam"

        return [self.bam_path, self.regions.path]

    def _define_outputs(self) -> List[str]:

        self.output_csv = f"{self.output_dir}/table.bedcov.csv"

        return [self.output_csv]

    def _run(self):
        """
        Count the number of fastq_files

        """

        # Run `samtools bedcov`

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
