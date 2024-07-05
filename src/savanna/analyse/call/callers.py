import subprocess
from abc import ABC, abstractmethod
from savanna.wrappers import bcftools


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


class VariantCaller(ABC):
    def __init__(self, fasta_path: str) -> None:
        self.fasta_path = fasta_path
        self.vcf_path = None

    @abstractmethod
    def _run(self, bam_path: str, vcf_path: str) -> None:
        pass

    def run(self, bam_path: str, vcf_path: str, sample_name: str = None):
        """
        Run core variant calling method
        """

        # Store
        self.vcf_path = vcf_path

        # Core method
        self._run(bam_path, vcf_path)

        # Optionally name
        if sample_name is not None:
            bcftools.reheader(self.vcf_path, self.vcf_path, [sample_name])

        # Index
        bcftools.index(self.vcf_path)

    def filter(
        self,
        output_vcf: str,
        bed_path: str,
        min_depth: int = 50,
        # min_qual: int = 20,
        to_biallelic: bool = False,
    ) -> None:
        """
        Filters the output VCF to only regions contained within
        `bed_path`; also optionally excludes sites below a threshold
        """

        if self.vcf_path is None:
            raise ValueError("Must run variant calling before filtering.")

        self.MIN_DEPTH = min_depth
        # self.MIN_QUAL = min_qual

        cmd_view = "bcftools view"
        cmd_view += f" -R {bed_path}"
        if to_biallelic:
            cmd_view += " --types='snps'"
            cmd_view += " --min-alleles 2"
            cmd_view += " --max-alleles 2"
        cmd_view += f" {self.vcf_path}"

        cmd_filter = "bcftools filter"
        cmd_filter += " -S ."
        cmd_filter += f" -e 'FORMAT/DP<{self.MIN_DEPTH}'"  # ||QUAL<{self.MIN_QUAL}'"
        cmd_filter += f" -Oz -o {output_vcf} -"

        cmd = f"{cmd_view} | {cmd_filter}"

        subprocess.run(cmd, check=True, shell=True)

        # Index
        bcftools.index(output_vcf)


# ================================================================
# Concrete implementations
#
# ================================================================


class BcfTools(VariantCaller):
    ANNOTATE_MPILEUP = "FORMAT/DP,FORMAT/AD"
    ANNOTATE_CALL = "FORMAT/GQ"
    MAX_DEPTH = 10_000

    def _run(self, bam_path: str, vcf_path: str) -> None:
        """
        Run variant calling with bcftools

        Note that here I am following recommendations of
        `bcftools` in calling with ONT reads, that is, using
        `-X ont` for `bcftools mpileup`, which sets:

        `-B`  : disable per-base alignment quality
        `-Q5` : Skip bases with base quality < 5
        `--max-BQ 30` : Set the maximum base quality to 30
            - ONT sets homopolymers to 90 for some reason
        `-I`  : Do not call any indels

        Then for `bcftools call`, I use `-P 0.01`,

        `-P` : Prior on mutation rate

        """

        cmd_pileup = "bcftools mpileup -Ou"
        cmd_pileup += f" -X ont"
        cmd_pileup += f" --annotate {self.ANNOTATE_MPILEUP}"
        cmd_pileup += f" --max-depth {self.MAX_DEPTH}"
        cmd_pileup += f" -f {self.fasta_path}"
        cmd_pileup += f" {bam_path}"

        # NB: We are returning *all* variants (not using -v)
        cmd_call = (
            f"bcftools call -m -P 0.01 -a '{self.ANNOTATE_CALL}' -Oz -o {vcf_path} -"
        )

        cmd = f"{cmd_pileup} | {cmd_call}"

        subprocess.run(cmd, check=True, shell=True)


# ================================================================
# Define collection of available callers
#
# ================================================================


# Note, they are already initialised
caller_collection = {"bcftools": BcfTools}
