import pandas as pd
from typing import List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from savanna.wrappers import samtools
from .callers import BcfTools
from .annotator import VariantAnnotator

# Run variant calling for the barcode
# make sure to return a gVCF (so we retain homozygous positios within the target regions)
class CallWithBcftools(BarcodeAnalysis):
    name = "call"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        reference: Reference = PlasmodiumFalciparum3D7(),
        make_plot: bool = True,
    ):
        self.Caller = BcfTools
        self.regions = regions
        self.reference = reference
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self):
        self.bam_path = (
            f"{self.barcode_dir}/bams/{self.barcode_name}.{self.reference.name}.bam"
        )
        self.fasta_path = self.reference.fasta_path
        self.gff_path = self.reference.gff_standard_path
        return [self.bam_path, self.fasta_path, self.gff_path, self.regions.path]

    def _define_outputs(self):
        self.vcf = (
            f"{self.output_dir}/bcftools.vcf.gz"  # not included in output list, removed
        )
        self.filtered_vcf = self.vcf.replace(".vcf.gz", ".filtered.vcf.gz")
        self.filtered_biallelic_vcf = self.filtered_vcf.replace(
            ".vcf.gz", ".biallelic.vcf.gz"
        )
        self.output_tsv = self.filtered_biallelic_vcf.replace(".vcf.gz", ".tsv")
        return [self.filtered_vcf, self.filtered_biallelic_vcf]

    def _run(self):
        caller = self.Caller(fasta_path=self.fasta_path)

        print("Calling variants...")
        caller.run(self.bam_path, self.vcf, sample_name=self.barcode_name)

        print("Filtering...")
        caller.filter(output_vcf=self.filtered_vcf, bed_path=self.regions.path)
        caller.filter(
            output_vcf=self.filtered_biallelic_vcf,
            bed_path=self.regions.path,
            to_biallelic=True,
        )

        annotator = VariantAnnotator(
            vcf_path=self.filtered_biallelic_vcf,
            bed_path=self.regions.path,
            reference=self.reference,
            output_dir=self.output_dir,
        )
        annotator.run()
        annotator.convert_to_tsv()

    def _plot(self):
        pass
