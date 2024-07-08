from typing import Dict, List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from .callers import VariantCaller
from .annotator import VariantAnnotator


class BarcodeVariantCaller(BarcodeAnalysis):
    name = "call"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        caller: VariantCaller,
        reference: Reference = PlasmodiumFalciparum3D7(),
        make_plot: bool = True,
    ):
        self.caller = caller
        self.reference = reference
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self) -> List[str]:
        """
        Define input files needed for amplicon variant calling
        """
        # BAM file
        bam_dir = f"{self.barcode_dir}/bams"
        bam_name = f"{self.barcode_name}.{self.reference.name}.filtered.bam"
        self.bam_path = f"{bam_dir}/{bam_name}"

        # FASTA file
        self.fasta_path = self.reference.fasta_path

        return [self.bam_path, self.fasta_path]

    def _define_outputs(self) -> List[str]:
        """
        Define outputs of amplicon variant calling
        """
        self.vcf = f"{self.output_dir}/{self.caller.name}.vcf.gz"
        return [self.vcf]

    def _run(self) -> None:
        """
        Run the variant calling method
        """
        self.caller.run(self.bam_path, self.vcf, sample_name=self.barcode_name)

    def _plot(self):
        pass
