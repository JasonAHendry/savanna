import subprocess
from abc import ABC, abstractmethod
from savanna.wrappers import bcftools


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


MIN_DEPTH = 50
MIN_QUAL = 20
TO_BIALLELIC_SNPS = True


class VariantCaller(ABC):
    """
    Interface for different variant calling methods

    Define the _run() method, which writes a vcf to
    self.vcf_path, to implement a new variant caller.

    The class includes useful auxiliary steps such as
    adding a sample name to the VCF and indexing it. It
    also has a filtering method useful for reducing
    a VCF to amplicon regions, with optional variant
    type filtering.
    """

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

    def filter_to_amplicons(
        self,
        output_vcf: str,
        bed_path: str,
        min_depth: int = MIN_DEPTH,
        min_qual: int = MIN_QUAL,
        to_biallelic: bool = TO_BIALLELIC_SNPS,
    ) -> None:
        """
        Filters the output VCF to only regions contained within
        `bed_path`; also optionally excludes sites below a threshold
        """
        if self.vcf_path is None:
            raise ValueError("Must run variant calling before filtering.")

        cmd_view = "bcftools view"
        cmd_view += f" -R {bed_path}"
        if to_biallelic:
            cmd_view += " --types='snps'"
            cmd_view += " --min-alleles 2"
            cmd_view += " --max-alleles 2"
        cmd_view += f" {self.vcf_path}"

        cmd_filter = "bcftools filter"
        cmd_filter += " -S ."
        cmd_filter += f" -e 'FORMAT/DP<{min_depth} || QUAL<{min_qual}'"
        cmd_filter += f" -Oz -o {output_vcf} -"

        cmd = f"{cmd_view} | {cmd_filter}"

        subprocess.run(cmd, check=True, shell=True)

        # Index
        bcftools.index(output_vcf)


# ================================================================
# Concrete implementations
#
# Guidelines:
# - Parameters we might commonly change are exposed in __init__
# - "Fixed" parameters are made class attributes for visibility
#
# ================================================================


# Why out here?
MAX_DEPTH = 10_000
MUTATION_RATE_PRIOR = 0.01


class BcfTools(VariantCaller):
    """
    Call variants using bcftools mpileup | bcftools call

    See:
        https://samtools.github.io/bcftools/howtos/variant-calling.html
    """

    ANNOTATE_MPILEUP = "FORMAT/DP,FORMAT/AD"
    ANNOTATE_CALL = "FORMAT/GQ"

    def __init__(self, 
                 fasta_path: str, 
                 max_depth: int = MAX_DEPTH,
                 mutation_rate_prior: float = MUTATION_RATE_PRIOR
    ) -> None:
        self.fasta_path = fasta_path
        self.max_depth = max_depth
        self.mutation_rate_prior = mutation_rate_prior
        self.vcf_path = None

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
        cmd_pileup += f" -X ont"
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


# Note, they are already initialised
CALLER_COLLECTION = {"bcftools": BcfTools}
