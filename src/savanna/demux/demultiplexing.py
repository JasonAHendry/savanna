import subprocess
from savanna.util.dirs import produce_dir
from abc import ABC, abstractmethod


BARCODING_KIT_MAPPING = {"native96": "SQK-NBD114-96", "rapid96": "SQK-RBK114-96"}


class Demultiplexer(ABC):
    """
    Demultiplex FASTQ data from a sequencing experiment;
    split into separate folders for each barcode

    TODO:
        - Input from either single file or directory
        - Can add a simple def _check_tool():
            - There is a self.tool = "guppy_barcoder"
            - check tool just runs subprocess.run(self.tool + " --help")
            - Will give better error message when not installed.
    """

    def __init__(self, fastq_dir: str, kit: str = "native96"):
        self.fastq_dir = fastq_dir
        self.barcode_kit = BARCODING_KIT_MAPPING[kit]

    @abstractmethod
    def run(self, output_dir: str):
        pass


class GuppyBarcoder(Demultiplexer):

    SCORE = 90

    def __init__(self, fastq_dir: str, kit: str = "native96"):
        super().__init__(fastq_dir, kit)

    def run(
        self,
        output_dir: str,
        use_gpu: bool = True,
        recursive: bool = True,
        both_ends: bool = False,
        strict: bool = True,
        trim_barcodes: bool = False,
        dry_run: bool = False,
    ):

        _ = produce_dir(output_dir)

        # Construct command
        cmd = "guppy_barcoder"
        if use_gpu:
            cmd += " --device 'cuda:0'"
        cmd += f" --barcode_kits {self.barcode_kit}"
        cmd += f" --input_path {self.fastq_dir}"
        if recursive:
            cmd += " --recursive"
        if both_ends:
            cmd += " --require_barcodes_both_ends"
        if strict:
            cmd += f" --min_score_barcode_front {self.SCORE}"
            cmd += f" --min_score_barcode_rear {self.SCORE}"
        cmd += f" --save_path {output_dir}"
        if trim_barcodes:
            cmd += " --enable_trim_barcodes"
        cmd += " --compress_fastq"
        cmd += " --disable_pings"

        if dry_run:
            print(cmd)
            return

        # Run
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            if e.returncode == 127:
                print(
                    "Error `guppy_barcoder` was not found. Ensure it is installed and available in your $PATH."
                )
            else:
                print(
                    f"Running `guppy_barcoder` failed with exit code: {e.returncode}."
                )
