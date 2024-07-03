import subprocess
from abc import ABC, abstractmethod


class Basecaller(ABC):
    """
    Demultiplex POD5 data from a sequencing experiment;

    """

    def __init__(self, pod5_dir: str):
        self.pod5_dir = pod5_dir

    @abstractmethod
    def run(self, output_fastq: str):
        pass


class Dorado(Basecaller):
    def __init__(self, pod5_dir: str, model: str) -> None:
        super().__init__(pod5_dir)

        self.model = model
        self.model_path = self._get_model_path(model)

    @staticmethod
    def _get_model_path(model):
        """
        Temporary, will resolve a better approach based on how the
        models get downloaded

        """
        return f"dorado_models/dna_r10.4.1_e8.2_400bps_{model}@v4.2.0"

    def run(self, output_fastq: str, min_qscore: int = 10, dry_run: bool = True):
        """
        
        Want to implement something like this:
        dorado basecaller $BASECALL_MODEL $POD5_DIR \
        --device 'cuda:0' \
        --min-qscore 10 \
        --emit-fastq \
        --recursive > $FASTQ_PATH
        
        """

        # Construct command
        cmd = "dorado basecaller"
        cmd += f" {self.model_path} {self.pod5_dir}"
        cmd += " --device 'cuda:0'"
        cmd += f" --min-qscore {min_qscore}"
        cmd += " --emit-fastq"
        cmd += f" --recursive > {output_fastq}"

        if dry_run:
            print(cmd)
            return

        # Run
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            if e.returncode == 127:
                print(
                    "Error `dorado` was not found. Ensure it is installed and available in your $PATH."
                )
            else:
                print(f"Running `dorado` failed with exit code: {e.returncode}.")
