import subprocess
import warnings
from typing import List


class VariantMerger:

    # For now, we only have two filtering criterion
    DEPTH_MIN = 50
    QUAL_MIN = 20

    def __init__(self, vcfs: List[str]):
        self.vcfs = vcfs
        self.vcf_string = " ".join(self.vcfs)

    def run(self, output_vcf: str):
        """
        Merge and filter an input VCF file

        """

        cmd_merge = f"bcftools merge {self.vcf_string}"

        warnings.warn("Output is currently limited to biallelic SNPs!")
        cmd_view = "bcftools view"
        cmd_view += " --min-alleles 2"
        cmd_view += " --max-alleles 2"
        cmd_view += " --types='snps'"
        cmd_view += f" -e 'MAX(FORMAT/DP)<{self.DEPTH_MIN}'"
        cmd_view += " - "  # Pipe

        cmd_filter = "bcftools filter -Oz"
        cmd_filter += " -S ."  # retain failed genotypes as ./.
        cmd_filter += f" -e 'FORMAT/DP<{self.DEPTH_MIN}'"
        cmd_filter += f" -o {output_vcf}"
        cmd_filter += " - "  # Pipe

        cmd = f"{cmd_merge} | {cmd_view} | {cmd_filter}"

        subprocess.run(cmd, check=True, shell=True)
