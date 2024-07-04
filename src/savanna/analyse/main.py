import logging

from savanna.util.logging_config import config_root_logger
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.download.references import REFERENCE_COLLECTION
from savanna.analyse.fastq.experiment import ExperimentFASTQCount
from savanna.analyse.map.experiment import (
    ExperimentMapToPfalciparum,
    ExperimentMapUnmappedToHSapiens,
)
from savanna.analyse.bedcov.experiment import ExperimentBedCoverage


def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str):
    """
    Entry point for core analysis pipeline

    """

    reference_name = "Pf3D7"

    # SETUP LOGGING
    config_root_logger(f"results/{expt_name}/metadata/savanna.log", verbose=False)
    log = logging.getLogger("savanna")

    # PARSE INPUT
    # log = LoggingFascade(logger_name="nomadic")  # , verbose=verbose)
    log.info("Input parameters:")
    log.info(f"  Experiment Name: {expt_name}")
    log.info(f"  FASTQ (.fastq): {fastq_dir}")
    log.info(f"  Metadata (.csv): {metadata_csv}")
    log.info(f"  Regions (.bed): {region_bed}")
    log.info(f"  Reference genome: {reference_name}")
    REFERENCE_COLLECTION[reference_name].confirm_downloaded()
    log.info("Processing...")

    # PREPARE TO RUN
    metadata = MetadataTableParser(metadata_csv)
    regions = RegionBEDParser(region_bed)
    expt_dirs = ExperimentDirectories(expt_name, metadata, regions, fastq_dir)
    log.info(f"  Found {len(metadata.barcodes) - 1} barcodes to track.")
    log.info(f"  Found {regions.n_regions} regions of interest.")
    log.info(f"  Outputs will be written to: {expt_dirs.expt_dir}.")
    log.info("Done.\n")

    # TODO
    # - Select and run the appropriate pipeline from CLI

    log.info("Counting FASTQ files...")
    analyse_fastq = ExperimentFASTQCount(expt_dirs, metadata)
    analyse_fastq.run()

    log.info("Mapping to P.f...")
    map_pf = ExperimentMapToPfalciparum(expt_dirs, metadata)
    map_pf.run()

    # log.info("Mapping to H.s....")
    # map_hs = ExperimentMapUnmappedToHSapiens(expt_dirs, metadata)
    # map_hs.run()

    log.info("Running BEDCOV analysis...")
    bedcov = ExperimentBedCoverage(expt_dirs, metadata, regions)
    bedcov.run()

    log.info("Done.")

