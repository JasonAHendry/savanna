import subprocess
from savanna.util.dirs import ROOT_DIR, produce_dir


class Dorado:
    """
    Basecalling using `dorado`

    Usage note:
        > "--no-trim  Skip trimming of barcodes, adapters, and primers. 
        If option is not chosen, trimming of all three is enabled."

    """

    MODEL_DIR = produce_dir(ROOT_DIR, "configs", "dorado")

    def __init__(self):
        pass

    def download_model(self, model: str):
        model_name = self._get_model_name(model)

        cmd = "dorado download"
        cmd += f" --model {model_name}"
        cmd += f" --directory {self.MODEL_DIR}"
        cmd += " --overwrite"
        self._run_cli_tool(cmd, "dorado")

    def _get_model_name(self, model: str):
        return f"dna_r10.4.1_e8.2_400bps_{model}@v4.2.0"

    def _get_model_path(self, model: str):
        return f"{self.MODEL_DIR}/{self._get_model_name(model)}"

    def run(
        self,
        model: str,
        pod5_dir: str,
        output_fastq: str,
        min_qscore: int = 9,
        dry_run: bool = False,
    ):

        model_path = self._get_model_path(model)

        cmd = "dorado basecaller"
        cmd += f" {model_path} {pod5_dir}"
        cmd += " --device 'cuda:0'"
        cmd += f" --min-qscore {min_qscore}"
        cmd += " --emit-fastq"
        cmd += " --no-trim"
        cmd += f" --recursive > {output_fastq}"

        if dry_run:
            print(cmd)
            return

        # Run
        self._run_cli_tool(cmd, "dorado")

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
