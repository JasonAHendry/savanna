import logging

from savanna.util.dirs import ROOT_DIR
from savanna.download.references import PlasmodiumFalciparum3D7, HomoSapiens
from savanna.analyse.bamfilt.experiment import ExperimentFilterBAM
from savanna.analyse.map.experiment import ExperimentMapUnmappedToHSapiens
from savanna.analyse.call.experiment import ExperimentCallWithMethod
from .interface import Pipeline


log = logging.getLogger()


class PlasmoPipeline(Pipeline):

    reference = PlasmodiumFalciparum3D7()
    hs_reference = HomoSapiens()
    parameter_path = f"{ROOT_DIR}/configs/parameters/default.yml"

    def run(self):
        log.info(f"Running {self.__class__.__name__}")

        self.reference.confirm_downloaded()
        self._load_parameters(self.parameter_path)
        self._count_fastqs()
        self._map_to_reference(self.reference)

        log.info("Remapping unmapped reads to H.s. reference.")
        hs_mapper = ExperimentMapUnmappedToHSapiens(
            expt_dirs=self.expt_dirs,
            metadata=self.metadata,
            **self.kwargs
        )
        hs_mapper.run()
        log.info("Done.")

        self._calc_bamstat(reference=self.reference)
        self._calc_bamstat(reference=self.hs_reference)
        
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

