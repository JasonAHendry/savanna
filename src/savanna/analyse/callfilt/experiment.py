from typing import Dict

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from .barcode import BarcodeVcfFilter


# MERGE (without filtering)
# FILTER --> ANNOTATE

# OUTPUTS
# complete.vcf.gz
# filtered.vcf.gz
# filtered.tsv
#


class ExperimentVcfFilter(ExperimentAnalysis):
    name = "call"

    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        caller_name: str,  # for now, only to get VCF path... need a better strategy
        to_snps: bool = True,
        to_biallelic: bool = True,
        min_depth: int = 50,
        min_qual: int = 15,
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
    ):
        # Inputs
        self.regions = regions
        self.caller_name = caller_name

        # Parameters
        self.filter_params = {
            "to_snps": to_snps,
            "to_biallelic": to_biallelic,
            "min_depth": min_depth,
            "min_qual": min_qual,
        }

        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeVcfFilter(
            barcode_name,
            expt_dirs=self.expt_dirs,
            regions=self.regions,
            caller_name=self.caller_name,
            **self.filter_params
        )

    def _summarise(self):
        pass

    def _plot(self):
        pass
