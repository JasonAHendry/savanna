import os
import subprocess
import warnings
from pathlib import Path
from savanna.util.dirs import produce_dir
from abc import ABC, abstractmethod


BARCODING_KITS = ["SQK-NBD114-96", "SQK-RBK114-96"]                        


class Demultiplexer(ABC):
    """
    Interface for demultiplexing tools

    """

    def __init__(self, fastq_dir: str, kit: str = BARCODING_KITS[0]):
        self.fastq_dir = fastq_dir
        self.barcode_kit = kit

    @abstractmethod
    def run(self, output_dir: str):
        pass

    @staticmethod
    def _run_cli_tool(cmd, tool: str):
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            if e.returncode == 127:
                print(
                    f"Error `{tool}` was not found."
                    + "Ensure it is installed and available in your $PATH."
                )
            else:
                print(f"Running `{tool}` failed with exit code: {e.returncode}.")



class GuppyBarcoder(Demultiplexer):
    """
    Wrapper for `guppy_barcoder`
    
    """

    # Identity score (as a percent) for barcode matches
    SCORE = 90

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
        """
        Run `guppy_barcoder`

        """

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
        self._run_cli_tool(cmd, tool="guppy_barcoder")


class DoradoDemux(Demultiplexer):
    """
    Wrapper for `dorado demux`
    
    """

    def _get_all_fastq(self):
        """
        Get all files ending with `.fastq` or `.fastq.gz` underneath
        the `.fastq` directory

        """
        fastqs = []
        for suffix in [".fastq", ".fastq.gz"]:
            fastqs.extend(Path(self.fastq_dir).rglob(f"**/*{suffix}"))
        return " ".join([str(f) for f in fastqs])                

    def run(
        self,
        output_dir: str,
        both_ends: bool = False,
        strict: bool = True,
        trim_barcodes: bool = False,
        dry_run: bool = False,
    ):
        """
        Run `guppy_barcoder`

        """

        _ = produce_dir(output_dir)
        reads = self._get_all_fastq()

        cmd = "dorado demux"
        cmd += f" --kit-name {self.barcode_kit}"
        if both_ends:
            cmd += f" --barcode-both-ends"
        if strict:
            warnings.warn("`--strict` is not yet implemented for `dorado demux`.")
        if not trim_barcodes:
            cmd += f" --no-trim"
        cmd += f" {output_dir}"
        cmd += f" {reads}"

        if dry_run:
            print(cmd)
            return
        
        # Run
        self._run_cli_tool(cmd, tool="dorado")



DEMUXER_COLLECTION = {
    "guppy": GuppyBarcoder,
    "dorado": DoradoDemux
}

