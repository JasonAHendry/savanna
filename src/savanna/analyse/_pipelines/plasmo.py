import logging

from savanna.download.references import PlasmodiumFalciparum3D7
from savanna.analyse.fastq.experiment import ExperimentFASTQCount
from savanna.analyse.call.experiment import ExperimentCallWithBcftools
from .interface import Pipeline


log = logging.getLogger()


class PlasmoPipeline(Pipeline):

    reference = PlasmodiumFalciparum3D7()

    def run(self):
        log.info(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()
        self._count_fastqs()
        self._map_to_reference(self.reference)
        self._calc_bamstat(self.reference)
        self._filter_bam(self.reference)
        self._calc_bedcov(self.reference)
        self._call_with_bcftools(self.reference)
