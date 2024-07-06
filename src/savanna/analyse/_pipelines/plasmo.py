import logging

from savanna.util.dirs import ROOT_DIR
from savanna.download.references import PlasmodiumFalciparum3D7
from savanna.analyse.bamfilt.experiment import ExperimentFilterBAM
from .interface import Pipeline


log = logging.getLogger()


class PlasmoPipeline(Pipeline):

    reference = PlasmodiumFalciparum3D7()
    parameter_path = f"{ROOT_DIR}/configs/parameters/plasmo-default.yml"

    def run(self):
        log.info(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()
        self._load_parameters(self.parameter_path)
        self._count_fastqs()
        self._map_to_reference(self.reference)
        
        log.info("Filtering BAM file")
        log.info(f"  Reference: {self.reference.name}")
        filter = ExperimentFilterBAM(
            self.expt_dirs, 
            self.metadata,
            self.reference,
            min_mapq=self.params["bam_filter"]["min_mapq"],
            exclude_flags=self.params["bam_filter"]["exclude_flags"],
            **self.kwargs
        )
        filter.run()
        log.info("Done.")

        self._calc_bedcov(self.reference)
        self._call_with_bcftools(self.reference)
