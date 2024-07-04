import logging

from savanna.util.logging_config import config_root_logger
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import REFERENCE_COLLECTION
from savanna.analyse.fastq.experiment import ExperimentFASTQCount
from savanna.analyse.map.experiment import (
    ExperimentMapToReference,
    ExperimentMapToPfalciparum,
    ExperimentMapUnmappedToHSapiens,
)
from savanna.analyse.bedcov.experiment import ExperimentBedCoverage
from savanna.analyse.pipelines import PIPELINE_COLLECTION

def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str, pipeline: str):
    """
    Entry point for core analysis pipeline

    """

    # SETUP LOGGING
    config_root_logger(f"results/{expt_name}/metadata/savanna.log", verbose=False)
    log = logging.getLogger("savanna")

    # PARSE INPUT
    log.info("Input parameters:")
    log.info(f"  Experiment Name: {expt_name}")
    log.info(f"  FASTQ (.fastq): {fastq_dir}")
    log.info(f"  Metadata (.csv): {metadata_csv}")
    log.info(f"  Regions (.bed): {region_bed}")
    log.info(f"  Pipeline: {pipeline}")
    log.info("Processing...")

    # PREPARE TO RUN
    metadata = MetadataTableParser(metadata_csv)
    regions = RegionBEDParser(region_bed)
    expt_dirs = ExperimentDirectories(expt_name, metadata, regions, fastq_dir)
    log.info(f"  Found {len(metadata.barcodes) - 1} barcodes to track.")
    log.info(f"  Found {regions.n_regions} regions of interest.")
    log.info(f"  Outputs will be written to: {expt_dirs.expt_dir}.")
    log.info("Done.\n")

    # SELECT PIPELINE AND RUN
    SelectedPipeline = PIPELINE_COLLECTION[pipeline]
    pipeline_runner = SelectedPipeline(
        expt_dirs,
        metadata,
        regions,
        only_barcode=None,
        only_summary=False,
        make_plot=False,
    )
    pipeline_runner.run()


    # reference_name = "Pf3D7"
    # REFERENCE_COLLECTION[reference_name].confirm_downloaded()

    # # TODO
    # # - Select and run the appropriate pipeline from CLI

    # log.info("Counting FASTQ files...")
    # analyse_fastq = ExperimentFASTQCount(expt_dirs, metadata)
    # analyse_fastq.run()

    # #log.info("Mapping to P.f...")
    # reference = REFERENCE_COLLECTION["AgPEST"]
    # mapper = ExperimentMapToReference(
    #     expt_dirs,
    #     metadata,
    #     reference
    # )
    # mapper.run()
    # map_pf = ExperimentMapToPfalciparum(expt_dirs, metadata)
    # map_pf.run()

    # # log.info("Mapping to H.s....")
    # # map_hs = ExperimentMapUnmappedToHSapiens(expt_dirs, metadata)
    # # map_hs.run()

    # log.info("Running BEDCOV analysis...")
    # # -> Now this also needs to take a Reference
    # bedcov = ExperimentBedCoverage(expt_dirs, metadata, regions, reference)
    # bedcov.run()

    # log.info("Done.")

