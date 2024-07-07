import logging
import subprocess
from typing import List
from .callers import MIN_DEPTH, MIN_QUAL, TO_BIALLELIC_SNPS


log = logging.getLogger()


class VariantMerger:
    """
    Merge variants from different barcodes

    When we filter VCFs from each barcode, we return a "."
    for those sites that fail to meet our filtering criterion. We still
    retain, however, all bi-allelic SNP sites across the entire amplicon region.

    Here, after merging, we can remove sites that were homozygous reference for
    all variants (e.g. have only one allele) or sites that did not pass our
    filtering criterion in *any* sample.

    For the later, it would best to use the genontype being "./." as an 
    indicator that the site failed to pass filtering in all samples. But
    need to learn how to do this.

    """

    def __init__(self, vcfs: List[str]):
        self.vcfs = vcfs
        self.vcf_string = " ".join(self.vcfs)

    def run(self, output_vcf: str, 
            min_depth: int = MIN_DEPTH,
            min_qual: int = MIN_QUAL,
            to_biallelic: str = TO_BIALLELIC_SNPS
    ):
        """
        Merge and filter an input VCF file

        """

        cmd_merge = f"bcftools merge {self.vcf_string}"

        cmd_view = "bcftools view"
        if to_biallelic:
            cmd_view += " --types='snps'"
            cmd_view += " --min-alleles 2"
            cmd_view += " --max-alleles 2"
        cmd_view += f" -e 'MAX(FORMAT/DP)<{min_depth} || QUAL<{min_qual}'"
        cmd_view += f" -o {output_vcf}"  # Pipe

        # I don't believe this is needed; I just want to actively
        # remove variants that haven't met filtering criterion
        # for any sample; I don't need to filter again as
        # this has been done per-barcode
        # cmd_filter = "bcftools filter -Oz"
        # cmd_filter += " -S ."  # retain failed genotypes as ./.
        # cmd_filter += f" -e 'FORMAT/DP<{min_depth}'"
        # cmd_filter += f" -o {output_vcf}"
        # cmd_filter += " - "  # Pipe

        cmd = f"{cmd_merge} | {cmd_view}"

        subprocess.run(cmd, check=True, shell=True)
