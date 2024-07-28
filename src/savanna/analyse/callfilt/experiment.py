import os
import subprocess
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference
from savanna.wrappers import bcftools
from savanna.analyse._interfaces import BarcodeAnalysis, ExperimentAnalysis
from .annotator import VariantAnnotator
from .barcode import BarcodeVcfFilter


class ExperimentVcfFilter(ExperimentAnalysis):
    name = "call"

    def __init__(
        self,
        expt_dirs: ExperimentDirectories,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        reference: Reference,
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
        self.reference = reference
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
            reference=self.reference,
            caller_name=self.caller_name,
            **self.filter_params
        )

    def _summarise(self):
        """
        Merge filtered VCFs and annotate

        Outputs:
        a) A merged VCF, not filtered
        b) A merged VCF
            - retaining variants wth at least one PASS and two alleles
        c) An annotated VCF and TSV of b)
        """
        # Collect VCF paths
        vcfs = []
        for _, (_, outputs) in self.results.items():
            filtered_vcf = outputs[0] 
            vcfs.append(filtered_vcf)

        # Merge all the VCFs
        output_vcf = f"{self.summary_dir}/{self.caller_name}.vcf.gz"
        bcftools.merge(vcfs, output_vcf)
        bcftools.index(output_vcf)

        # Filter to variant passing sites before anntoating VCF
        filtered_vcf = output_vcf.replace(".vcf.gz", ".filtered.vcf.gz")
        cmd = (
            "bcftools view"
            " --apply-filters PASS"
            " --min-alleles 2"
            f" -Oz -o {filtered_vcf}"
            f" {output_vcf}"
        )
        subprocess.run(cmd, shell=True, check=True)

        # Annotate
        annotator = VariantAnnotator(
            filtered_vcf,  # Note that we annotated only filtered VCF
            self.regions.path,
            self.reference,
            output_dir=self.summary_dir,
        )
        annotator.run()
        annotator.convert_to_tsv()

        # Once annotation completed, can remove unnanotated version
        os.remove(filtered_vcf)

    def _plot(self):
        pass
