from typing import List
from savanna.analyse._interfaces import BarcodeAnalysis
from savanna.util.dirs import ExperimentDirectories
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from savanna.wrappers import samtools


MIN_MAPQ = 40
EXCLUDE_FLAGS = "0x904"  # Only retain primary alignments


class FilterBAM(BarcodeAnalysis):
    name = "bams"  # we want to write to this folder

    def __init__(
        self,
        barcode_name: str,
        expt_dirs: ExperimentDirectories,
        reference: Reference = PlasmodiumFalciparum3D7(),
        make_plot: bool = True,
        min_mapq: int = MIN_MAPQ,
        exclude_flags: str = EXCLUDE_FLAGS
    ):
        self.reference = reference
        self.min_mapq = min_mapq
        self.exclude_flags = exclude_flags
        super().__init__(barcode_name, expt_dirs, make_plot)

    def _define_inputs(self) -> List[str]:
        self.input_bam = (
            f"{self.barcode_dir}/bams/{self.barcode_name}.{self.reference.name}.bam"
        )
        return [self.input_bam]

    def _define_outputs(self) -> List[str]:
        self.output_bam = self.input_bam.replace(".bam", ".filtered.bam")
        return [self.output_bam]

    def _run(self):
        """
        Notes on BAM Flags
        0x4     : segment unmapped
        0x100   : secondary alignment
        0x800   : chimeric / supplementary alignment
        """

        args = f" -q {self.min_mapq} -F {self.exclude_flags}"
        samtools.view(input_bam=self.input_bam, args=args, output_bam=self.output_bam)
        samtools.index(self.output_bam)

    def _plot(self):
        pass
