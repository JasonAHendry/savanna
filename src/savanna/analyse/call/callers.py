import subprocess
from abc import ABC, abstractmethod
from savanna.wrappers import bcftools


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


class VariantCaller(ABC):
    """
    Interface for different variant calling methods

    """

    def __init__(self, fasta_path: str) -> None:
        self.fasta_path = fasta_path

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Give a short name to the variant caller
        """
        pass

    @abstractmethod
    def _run(self, bam_path: str, vcf_path: str) -> None:
        """
        Produce a VCF from a BAM file
        """
        pass

    def run(self, bam_path: str, vcf_path: str, sample_name: str = None):
        """
        Run core variant calling method
        """
        # Core method
        self._run(bam_path, vcf_path)

        # Optionally name
        if sample_name is not None:
            bcftools.reheader(vcf_path, vcf_path, [sample_name])

        # Index
        bcftools.index(vcf_path)


# ================================================================
# Concrete implementations
#
# Guidelines:
# - Parameters we might commonly change are exposed in __init__
# - "Fixed" parameters are made class attributes for visibility
#
# ================================================================


class BcfTools(VariantCaller):
    """
    Call variants using bcftools mpileup | bcftools call

    See:
        https://samtools.github.io/bcftools/howtos/variant-calling.html
    """

    ANNOTATE_MPILEUP = "FORMAT/DP,FORMAT/AD"
    ANNOTATE_CALL = "FORMAT/GQ"

    def __init__(
        self,
        fasta_path: str,
        mpileup_setting: str = "ont",
        max_depth: int = 10_000,
        mutation_rate_prior: float = 0.01,
    ) -> None:
        self.fasta_path = fasta_path
        self.mpileup_setting = mpileup_setting
        self.max_depth = max_depth
        self.mutation_rate_prior = mutation_rate_prior

    @property
    def name(self):
        return "bcftools"

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

        From `bcftools mpileup`
        ont:         -B -Q5 --max-BQ 30 -I [also try eg |bcftools call -P0.01
        """
        # Prepare mpileup command
        cmd_pileup = "bcftools mpileup -Ou"
        cmd_pileup += f" -X {self.mpileup_setting}"
        cmd_pileup += f" --annotate {self.ANNOTATE_MPILEUP}"
        cmd_pileup += f" --max-depth {self.max_depth}"
        cmd_pileup += f" -f {self.fasta_path}"
        cmd_pileup += f" {bam_path}"

        # Prepare call command
        #  - returning *all* variants (not using -v)
        #  - using multiallelic model (-m)
        cmd_call = "bcftools call -m"
        cmd_call += f" -P {self.mutation_rate_prior}"
        cmd_call += f" -a '{self.ANNOTATE_CALL}'"
        cmd_call += f" -Oz -o {vcf_path} -"

        cmd = f"{cmd_pileup} | {cmd_call}"

        subprocess.run(cmd, check=True, shell=True)


# ================================================================
# Define collection of available callers
#
# ================================================================


CALLER_COLLECTION = {"bcftools": BcfTools}
