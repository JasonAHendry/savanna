import logging

from savanna.download.references import AnophelesGambiaePEST
from savanna.analyse.fastq.experiment import ExperimentFASTQCount
from savanna.analyse.call.experiment import ExperimentCallWithBcftools
from .interface import Pipeline


log = logging.getLogger()


class EntoPipeline(Pipeline):

    reference = AnophelesGambiaePEST()

    def run(self):
        log.info(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()
        self._count_fastqs()
        self._map_to_reference(self.reference)
        self._calc_bamstat(self.reference)
        self._calc_bedcov(self.reference)
        self._call_with_bcftools(self.reference)
