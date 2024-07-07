from typing import Dict, List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from .callers import CALLER_COLLECTION, MIN_DEPTH, MIN_QUAL, TO_BIALLELIC_SNPS
from .annotator import VariantAnnotator


class CallWithMethod(BarcodeAnalysis):
    name = "call"
    def __init__(self,
                 barcode_name: str,
                 expt_dirs: ExperimentDirectories,
                 regions: RegionBEDParser,
                 caller_name: str,
                 min_depth: int = MIN_DEPTH,
                 min_qual: int = MIN_QUAL,
                 to_biallelic: bool = TO_BIALLELIC_SNPS,
                 reference: Reference = PlasmodiumFalciparum3D7(),
                 make_plot: bool = True,
                 **caller_params: Dict):
        """
        Select the appropriate VariantCaller and store the parameters
        """
        self.caller_name = caller_name
        self.Caller = CALLER_COLLECTION[caller_name]
        self.caller_params = caller_params

        self.regions = regions
        self.reference = reference

        self.min_depth = min_depth
        self.min_qual = min_qual
        self.to_biallelic = to_biallelic

        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self) -> List[str]:
        """
        Define input files needed for amplicon variant calling
        """
        # BAM file
        self.bam_dir = f"{self.barcode_dir}/bams"
        self.bam_path = f"{self.bam_dir}/{self.barcode_name}.{self.reference.name}.filtered.bam"

        # FASTA file
        self.fasta_path = self.reference.fasta_path

        # GFF file (used for filtering)
        self.gff_path = self.reference.gff_standard_path

        return [self.bam_path, self.fasta_path, self.gff_path, self.regions.path]
    
    def _define_outputs(self) -> List[str]:
        """
        Define outputs of amplicon variant calling
        """
        core = f"{self.output_dir}/{self.caller_name}"
        self.vcf = f"{core}.vcf.gz"
        self.filtered_vcf = f"{core}.filtered.vcf.gz"
        self.output_tsv = f"{core}.filtered.tsv"
        return [self.vcf, self.filtered_vcf, self.output_tsv]
    
    def _run(self) -> None:
        """
        Run the variant calling method
        """
        # Initialise caller with parameters and run
        caller = self.Caller(
            fasta_path=self.fasta_path,
            **self.caller_params
        )
        caller.run(
            self.bam_path, 
            self.vcf, 
            sample_name=self.barcode_name
        )

        # Filter to amplicos only (before merging)
        caller.filter_to_amplicons(
            output_vcf=self.filtered_vcf,
            bed_path=self.regions.path,
            min_depth=self.min_depth,
            min_qual=self.min_qual,
            to_biallelic=self.to_biallelic
        )

        # Annotate here as well, in case
        # want to look at a barcode alone
        annotator = VariantAnnotator(
            vcf_path=self.filtered_vcf,
            bed_path=self.regions.path,
            reference=self.reference,
            output_dir=self.output_dir,
        )
        annotator.run()
        annotator.convert_to_tsv()

    def _plot(self):
        pass
