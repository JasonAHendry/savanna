import os
import subprocess
from typing import List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import Reference
from savanna.wrappers import bcftools


class BarcodeVcfFilter(BarcodeAnalysis):
    name = "call"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        regions: RegionBEDParser,
        reference: Reference,
        caller_name: str,  # for now, only to get VCF path... need a better strategy
        to_snps: bool = True,
        to_biallelic: bool = True,
        min_depth: int = 50,
        min_qual: int = 15,
        make_plot: bool = True,
    ):
        # Inputs
        self.regions = regions
        self.reference = reference
        self.caller_name = caller_name

        # Parameters
        self.to_snps = to_snps
        self.to_biallelic = to_biallelic
        self.min_depth = min_depth
        self.min_qual = min_qual

        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self) -> List[str]:
        self.vcf_path = f"{self.barcode_dir}/call/{self.caller_name}.vcf.gz"
        self.bed_path = self.regions.path
        return [self.vcf_path]

    def _define_outputs(self) -> List[str]:
        self.filtered_vcf = self.vcf_path.replace(".vcf.gz", ".filtered.vcf.gz")
        return [self.filtered_vcf]
    
    def _get_lowcomplexity_filter_command(self) -> str:
        """
        Construct the command for low-complexity filtering
        of a VCF

        TODO:
        * Optional flag
        * Allow for input / outputs

        """

        bed_mask_path = self.expt_dir.regions_bed.replace(".bed", ".lowcomplexity_mask.bed")
        if not os.path.exists(bed_mask_path):
            print("Creating low-complexity mask for amplicons...")
            cmd = (
                "bedtools intersect"
                f" -a {self.reference.fasta_mask_path}"
                f" -b {self.expt_dir.regions_bed}"
                f" -wa > {bed_mask_path}"
            )
            subprocess.run(cmd, shell=True, check=True)
            print("Done.")

        # Masking command
        cmd = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowComplexity"
            " --set-GTs ."
            f" --mask-file {bed_mask_path}"
            " -Ou -"
        )

        return cmd

    def _run(self):
        """
        Reduce to sites within amplicons, soft-filter based on
        depth and quality
        """
        # Reduce to amplicons and particular variant types
        cmd_view = "bcftools view"
        cmd_view += f" -R {self.bed_path}"
        if self.to_snps:
            cmd_view += " --types='snps'"
        if self.to_biallelic:
            cmd_view += " --min-alleles 2"
            cmd_view += " --max-alleles 2"
        cmd_view += f" {self.vcf_path}"

        cmd_depth_filter = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowDepth"
            " --set-GTs ."
            f" --exclude 'FORMAT/DP < {self.min_depth}'"
            " -Ou - "
        )

        cmd_lowcomplexity_filter = self._get_lowcomplexity_filter_command()

        cmd_qual_filter = (
            "bcftools filter"
            " --mode +"
            " --soft-filter LowQual"
            " --set-GTs ."
            f" --exclude 'QUAL < {self.min_qual}'"
            f" -Oz -o {self.filtered_vcf} - "
        )

        cmd = f"{cmd_view} | {cmd_depth_filter} | {cmd_lowcomplexity_filter} |{cmd_qual_filter}"
        subprocess.run(cmd, shell=True, check=True)

        # Index
        bcftools.index(self.filtered_vcf)

    def _plot(self):
        pass
