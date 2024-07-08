import logging

from savanna.util.dirs import ROOT_DIR
from savanna.download.references import AnophelesGambiaePEST
from savanna.analyse.call.experiment import ExperimentVariantCaller
from .interface import Pipeline


log = logging.getLogger()


class EntoPipeline(Pipeline):

    reference = AnophelesGambiaePEST()
    parameter_path = f"{ROOT_DIR}/configs/parameters/default.yml"

    def run(self):
        log.info(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()
        self._load_parameters(self.parameter_path)
        self._count_fastqs()
        self._map_to_reference(self.reference)
        self._calc_bamstat(self.reference)
        self._filter_bam(self.reference)
        self._calc_bedcov(self.reference)
        
        log.info("Calling variants with bcftools")
        log.info(f"  Reference: {self.reference.name}")
        log.info(f"  Calling parameters: {self.params['call']['bcftools']}")
        log.info(f"  Filtering parameters: {self.params['call_filter']}")
        caller = ExperimentCallWithMethod(
            expt_dirs=self.expt_dirs, 
            metadata=self.metadata,
            reference=self.reference,
            regions=self.regions,
            caller_name="bcftools",
            min_depth=self.params["call_filter"]["min_depth"],
            min_qual=self.params["call_filter"]["min_qual"],
            to_biallelic=self.params["call_filter"]["biallelic_snps"],
            **self.kwargs,
            **self.params["call"]["bcftools"]
        )
        caller.run()
        log.info("Done.")

