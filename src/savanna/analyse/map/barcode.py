import os
import subprocess
from typing import List
from abc import ABC, abstractmethod

from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.download.references import Reference, PlasmodiumFalciparum3D7, HomoSapiens

# This is a reasonably large file.
# At top, I have classes encapsulating mapping algorithms.
# Below, these classes are used in some BarcodeAnalysis
# subclasses.


# --------------------------------------------------------------------------------
# Define abstract base class for different mapping algorithms
#
# --------------------------------------------------------------------------------


class MappingAlgorithm(ABC):
    """
    Abstract base class for implementing different mapping algorithms.

    """

    def __init__(
        self,
        fasta_path: str,
        *,
        fastq_path: str = None,
        fastq_dir: str = None,
        input_bam: str = None,
    ):
        """
        Initialise `fasta_path`, `input_fastqs` and any `initial` processing

        """

        self.fasta_path = fasta_path
        self.inputs = [fastq_path, fastq_dir, input_bam]
        if sum(i is not None for i in self.inputs) != 1:
            raise ValueError(
                "Only one of `fastq_path`, `fastq_dir`, or `input_bam` can be provided."
            )

        self.initial = ""
        self.input_fastqs = "-"

        if input_bam is not None:
            self.initial = self._get_remap_command(input_bam)
        elif fastq_path is not None:
            self.input_fastqs = fastq_path
        else:
            self.input_fastqs = " ".join(self._get_all_fastqs(fastq_dir))

        self.output_bam = None

    def _get_remap_command(self, input_bam: str) -> str:
        """
        Get a command to remap unmapped reads from a given `input_bam`

        """
        return f"samtools view -f 0x004 {input_bam} | samtools fastq |"

    def _get_all_fastqs(self, fastq_dir: str) -> List[str]:
        """
        Get all the FASTQ files in a given `fastq_dir`

        """
        return [
            f"{fastq_dir}/{fastq}"
            for fastq in os.listdir(fastq_dir)
            if fastq.endswith(".fastq") or fastq.endswith(".fastq.gz")
        ]

    @abstractmethod
    def _define_mapping_command(self, output_bam: str, flags: str) -> str:
        """
        Abstract method to define the command line for the mapping algorithm.

        """

    def run(self, output_bam, flags="", dry_run: bool = False):
        """
        Executes the mapping algorithm with the provided parameters.

        Constructs and optionally executes the mapping command using the
        implementation provided by the subclass.

        """

        map_cmd = self._define_mapping_command(output_bam, flags)
        cmd = self.initial + " " + map_cmd

        if dry_run:
            return cmd

        subprocess.run(cmd, shell=True, check=True)
        self.output_bam = output_bam

    def index(self):
        """
        Index the *.bam file generated by mapping

        """
        cmd = f" samtools index {self.output_bam}"
        subprocess.run(cmd, check=True, shell=True)


# --------------------------------------------------------------------------------
# Concrete mapping algorithm implementations
#
# --------------------------------------------------------------------------------


class Minimap2(MappingAlgorithm):
    """
    Map long reads with `minimap2`

    """

    def _define_mapping_command(
        self, output_bam: str, flags: str = "--eqx --MD"
    ) -> str:
        """
        Run minimap2, compress result to .bam file, and sort

        """
        map_cmd = "minimap2"
        map_cmd += f" -ax map-ont {flags} {self.fasta_path} {self.input_fastqs} |"
        map_cmd += " samtools view -S -b - |"
        map_cmd += f" samtools sort -o {output_bam}"

        return map_cmd


class BwaMem(MappingAlgorithm):
    """
    Map short reads with `bwa mem`

    """

    def create_reference_index(self):
        """
        Create an index for bwa

        Produces a set of index files with the same prefix
        as `self.reference.fasta_path`.

        """
        index_cmd = f"bwa index {self.fasta_path}"
        subprocess.run(index_cmd, shell=True, check=True)

    def _define_mapping_command(self, output_bam: str, flags="") -> str:
        """
        Run bwa, compress result to .bam file, and sort

        """
        map_cmd = "bwa mem"
        map_cmd += " -R '@RG\\tID:misc\\tSM:pool'"  # ID and SM tags needed for gatk HaplotypeCaller
        map_cmd += f" {flags} {self.fasta_path} {self.input_fastqs} |"
        map_cmd += " samtools view -S -b - |"
        map_cmd += f" samtools sort -o {output_bam}"

        return map_cmd


# --------------------------------------------------------------------------------
# Create a collection of mapping algorithms
#
# --------------------------------------------------------------------------------


MAPPER_COLLECTION = {"minimap2": Minimap2, "bwa": BwaMem}


# --------------------------------------------------------------------------------
# Create BarcodeAnalysis steps for mapping
#
# --------------------------------------------------------------------------------


# Start of  drafting specific, can add generality later if needed
# Generality would be good actually..


class BarcodeMapToReference(BarcodeAnalysis):
    """
    Map reads to a single reference genome

    """

    name = "bams"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        reference: Reference,
        make_plot: bool = True,
    ):
        self.reference = reference
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self) -> List[str]:
        self.fastq_dir = f"{self.expt_dir.demux_dir}/{self.barcode_name}"
        self.fasta_path = self.reference.fasta_path
        return [self.fastq_dir, self.fasta_path]

    def _define_outputs(self) -> List[str]:
        self.output_bam = (
            f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.bam"
        )
        return [self.output_bam]

    def _run(self):
        """
        Run mapping from FASTQ data using Minimap2

        """

        # Confirm reference genome downloaded
        self.reference.confirm_downloaded()

        # Prepare the mapper
        mapper = Minimap2(fasta_path=self.fasta_path, fastq_dir=self.fastq_dir)

        # Run
        mapper.run(self.output_bam)
        mapper.index()

    def _plot(self):
        pass


class BarcodeMapToPfalciparum(BarcodeAnalysis):
    """
    Map reads from FASTQ data

    """

    name = "bams"
    reference = PlasmodiumFalciparum3D7()

    def _define_inputs(self) -> List[str]:
        """
        Define input files for mapping algorithm

        """

        self.fastq_dir = f"{self.expt_dir.demux_dir}/{self.barcode_name}"

        return [self.fastq_dir]

    def _define_outputs(self) -> List[str]:

        self.output_bam = (
            f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.bam"
        )

        return [self.output_bam]

    def _run(self):
        """
        Run mapping from FASTQ data using Minimap2

        """

        # Confirm reference genome downloaded
        self.reference.confirm_downloaded()

        # Prepare the mapper
        mapper = Minimap2(
            fasta_path=self.reference.fasta_path, fastq_dir=self.fastq_dir
        )

        # Run
        mapper.run(self.output_bam)
        mapper.index()

    def _plot(self):
        pass


class BarcodeMapUnmappedToHSapiens(BarcodeAnalysis):
    """
    Map unmapped reads to the human host

    """

    name = "bams"
    reference = HomoSapiens()

    def _define_inputs(self) -> List[str]:
        """
        Define input files for mapping algorithm

        """

        self.pf_bam = f"{self.output_dir}/{self.barcode_name}.Pf3D7.bam"

        return [self.pf_bam]

    def _define_outputs(self) -> List[str]:

        self.output_bam = (
            f"{self.output_dir}/{self.barcode_name}.{self.reference.name}.bam"
        )

        return [self.output_bam]

    def _run(self):
        """
        Run mapping from FASTQ data using Minimap2

        """

        # Confirm reference genome downloaded
        self.reference.confirm_downloaded()

        # Prepare the mapper
        mapper = Minimap2(fasta_path=self.reference.fasta_path, input_bam=self.pf_bam)

        # Define output bam and run
        mapper.run(self.output_bam)
        mapper.index()

    def _plot(self):
        pass
