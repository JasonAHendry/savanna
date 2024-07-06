import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser

from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from .barcode import CallWithBcftools
from .merger import VariantMerger
from .annotator import VariantAnnotator


# Goals:
# - Carefully review filtering
# - Make sure all of the parameters used are made explicit somehow
# - Would be good to improve logging
# Merge all of the VCF files
# Annotate
# TSV
# Make a basic plot
# - Heatmap
# - Dot plot
# - Interactive plot (like a plotly)


class ExperimentCallWithBcftools(ExperimentAnalysis):
    name = "call"
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
        return CallWithBcftools(
            barcode_name, self.expt_dirs, self.regions, self.reference
        )

    def _summarise(self):
        """Merge the VCF files and convert to a TSV"""

        # Collect VCF paths
        vcfs = []
        for b, (success, outputs) in self.results.items():
            filtered_vcf = outputs[0]
            vcfs.append(filtered_vcf)

        # Merge VCFs
        output_vcf = f"{self.summary_dir}/bcftools.filtered.biallelic.vcf.gz"
        merger = VariantMerger(vcfs)
        merger.run(output_vcf)

        # Annotate
        annotator = VariantAnnotator(
            output_vcf,  # Note that we annotated only filtered VCF
            self.regions.path,
            self.reference,
            output_dir=self.summary_dir,
        )
        annotator.run()
        annotator.convert_to_tsv()

    def _plot(self):
        pass
