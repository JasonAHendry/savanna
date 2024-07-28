import logging

from savanna.util.dirs import ROOT_DIR
from savanna.download.references import Reference, PlasmodiumFalciparum3D7
from savanna.analyse.bamfilt.experiment import ExperimentFilterBAM
from savanna.analyse.call.callers import CALLER_COLLECTION
from savanna.analyse.call.experiment import ExperimentVariantCaller
from savanna.analyse.callfilt.experiment import ExperimentVcfFilter
from .interface import Pipeline


log = logging.getLogger()


class PlasmoPipeline(Pipeline):
    """
    I think what would be nicer would be to define
    methods for each analysis HERE;

    Steps that are generic and could be used across pipelines:
    * Put them in class Pipeline(ABC)

    Steps specific to this pipeline:
    * Put them in this pipeline...

    Key challenge is handling of parameters

    """

    reference = PlasmodiumFalciparum3D7()
    parameter_path = f"{ROOT_DIR}/configs/parameters/default.yml"

    def _filter_bam(self, reference: Reference):
        log.info("Filtering BAM file")
        log.info(f"  Reference: {self.reference.name}")
        filter = ExperimentFilterBAM(
            expt_dirs=self.expt_dirs,
            metadata=self.metadata,
            reference=reference,
            min_mapq=self.params["bam_filter"]["min_mapq"],
            exclude_flags=self.params["bam_filter"]["exclude_flags"],
            **self.kwargs
        )
        filter.run()
        log.info("Done.")

    def _call_variants(self):
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
                expt_dirs=self.expt_dirs,
                metadata=self.metadata,
                reference=self.reference,
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

    def run(self):
        """
        TODO: I believe this arrangment is a lot more readable,
        but need to look again at what I did in Nomadic

        """
        log.info("-"*80)
        log.info(f"Running {self.__class__.__name__}")
        log.info("-"*80)
        
        self.reference.confirm_downloaded()
        self._load_parameters(self.parameter_path)
        self._count_fastqs()
        self._map_to_reference(self.reference)
        self._calc_bamstat(self.reference)
        self._filter_bam(self.reference)
        self._calc_bedcov(self.reference)
        self._run_experiment_qc()
        self._call_variants()


