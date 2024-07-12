import logging

from savanna.util.dirs import ROOT_DIR
from savanna.download.references import PlasmodiumFalciparum3D7, HomoSapiens
from savanna.analyse.bamfilt.experiment import ExperimentFilterBAM
from savanna.analyse.map.experiment import ExperimentMapUnmappedToHSapiens
from savanna.analyse.call.callers import CALLER_COLLECTION
from savanna.analyse.call.experiment import ExperimentVariantCaller
from savanna.analyse.callfilt.experiment import ExperimentVcfFilter
from .interface import Pipeline


log = logging.getLogger()


class PlasmoPipeline(Pipeline):

    reference = PlasmodiumFalciparum3D7()
    parameter_path = f"{ROOT_DIR}/configs/parameters/default.yml"

    def run(self):
        log.info(f"Running {self.__class__.__name__}")

        # For convenience below
        core_args = {
            "expt_dirs": self.expt_dirs,
            "metadata": self.metadata,
            "reference": self.reference
        }

        self.reference.confirm_downloaded()
        self._load_parameters(self.parameter_path)
        self._count_fastqs()
        self._map_to_reference(self.reference)
        self._calc_bamstat(reference=self.reference)
        
        log.info("Filtering BAM file")
        log.info(f"  Reference: {self.reference.name}")
        filter = ExperimentFilterBAM(
            **core_args,
            min_mapq=self.params["bam_filter"]["min_mapq"],
            exclude_flags=self.params["bam_filter"]["exclude_flags"],
            **self.kwargs
        )
        filter.run()
        log.info("Done.")

        self._calc_bedcov(self.reference)

        # Iterate over variant calling methods
        log.info("Calling variants:")
        for caller_name, caller_params in self.params["call"].items():
            log.info(f"  Tool: {caller_name}")
            log.info(f"  Reference: {self.reference.name}")
            log.info(f"  Calling parameters: {caller_params}")
            Caller = CALLER_COLLECTION[caller_name]
            caller = Caller(
                fasta_path=self.reference.fasta_path,
                **caller_params
            )
            expt_caller = ExperimentVariantCaller(
                **core_args,
                caller=caller,
                **self.kwargs
            )
            expt_caller.run()

            log.info(f"Filtering variants: {self.params['call_filter']}")
            expt_filter = ExperimentVcfFilter(
                self.expt_dirs,
                self.metadata,
                self.regions,
                self.reference,
                caller_name=caller_name,
                **self.params['call_filter'],
                **self.kwargs
            )
            expt_filter.run()
            log.info("Done with tool.")
        log.info("Done with all variant calling.")
        
