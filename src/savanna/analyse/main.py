import logging

from savanna.util.logging_config import config_root_logger
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser
from savanna.analyse.pipelines import PIPELINE_COLLECTION

def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str, pipeline: str, barcode: int, summary_only: bool):
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
    if barcode is not None:
        log.info(f"  Processing exclusively barcode: {barcode}")
    log.info(f"  Only producing summary: {summary_only}")
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
        barcode=barcode,
        summary_only=summary_only
    )
    pipeline_runner.run()

