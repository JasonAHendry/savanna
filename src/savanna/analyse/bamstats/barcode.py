import json
from typing import List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from savanna.wrappers import samtools


class BamStats(BarcodeAnalysis):
    name = "bamstats"

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        reference: Reference = PlasmodiumFalciparum3D7(),
        make_plot: bool = True,
    ):
        self.reference = reference
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self):
        self.bam_path = (
            f"{self.barcode_dir}/bams/{self.barcode_name}.{self.reference.name}.bam"
        )
        return [self.bam_path]

    def _define_outputs(self) -> List[str]:
        self.output_json = f"{self.output_dir}/bamstats.{self.reference.name}.json"
        return [self.output_json]

    def _run(self):
        samtools.flagstats(input_bam=self.bam_path, output_json=self.output_json)
        dt = json.load(open(self.output_json, "r"))
        dt["barcode"] = self.barcode_name
        json.dump(dt, open(self.output_json, "w"))

    def _plot(self):
        pass
