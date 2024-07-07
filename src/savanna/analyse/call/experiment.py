from typing import Dict

from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from savanna.download.references import Reference, PlasmodiumFalciparum3D7

from .callers import MIN_DEPTH, MIN_QUAL, TO_BIALLELIC_SNPS
from .barcode import CallWithMethod
from .merger import VariantMerger
from .annotator import VariantAnnotator



class ExperimentCallWithMethod(ExperimentAnalysis):
    name = "call"
    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        caller_name: str,
        min_depth: int = MIN_DEPTH,
        min_qual: int = MIN_QUAL,
        to_biallelic: bool = TO_BIALLELIC_SNPS,
        reference: Reference = PlasmodiumFalciparum3D7(),
        barcode: str = None,
        summary_only: bool = False,
        make_plot: bool = True,
        **caller_params: Dict
    ):
        self.regions = regions
        self.reference = reference

        self.min_depth = min_depth
        self.min_qual = min_qual
        self.to_biallelic = to_biallelic

        self.caller_name = caller_name
        self.caller_params = caller_params

        super().__init__(expt_dirs, metadata, barcode, summary_only, make_plot)

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return CallWithMethod(
            barcode_name=barcode_name, 
            expt_dirs=self.expt_dirs, 
            regions=self.regions, 
            reference=self.reference,
            caller_name=self.caller_name,
            min_depth=self.min_depth,
            min_qual=self.min_qual,
            to_biallelic=self.to_biallelic,
            **self.caller_params
        )
    
    def _summarise(self):
        """Merge the VCF files and convert to a TSV"""

        # Collect VCF paths
        vcfs = []
        for b, (success, outputs) in self.results.items():
            filtered_vcf = outputs[1] 
            vcfs.append(filtered_vcf)

        # Merge VCFs
        output_vcf = f"{self.summary_dir}/{self.caller_name}.filtered.biallelic.vcf.gz"
        merger = VariantMerger(vcfs)
        merger.run(output_vcf,
                   min_depth=self.min_depth,
                   min_qual=self.min_qual,
                   to_biallelic=self.to_biallelic
        )

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

