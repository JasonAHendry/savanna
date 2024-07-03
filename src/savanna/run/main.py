import logging

from savanna.gather.references import REFERENCE_COLLECTION

from savanna.demux.demultiplexing import GuppyBarcoder
from savanna.run.experiment.fastq import ExperimentFASTQCount
from savanna.run.experiment.map import (
    ExperimentMapToPfalciparum,
    ExperimentMapUnmappedToHSapiens,
)
from savanna.run.experiment.bedcov import ExperimentBedCoverage

from savanna.util.logging_config import config_root_logger  # , LoggingFascade
from savanna.util.metadata import MetadataTableParser
from savanna.util.dirs import ExperimentDirectories
from savanna.util.regions import RegionBEDParser


def main(expt_name: str, fastq_dir: str, metadata_csv: str, region_bed: str):
    """
    Entry point for main pipeline

    """

    reference_name = "Pf3D7"

    # SETUP LOGGING
    # - Need to reorder main somewhat
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
    expt_dirs = ExperimentDirectories(expt_name, metadata, regions)
    log.info(f"  Found {len(metadata.barcodes) - 1} barcodes to track.")
    log.info(f"  Found {regions.n_regions} regions of interest.")
    log.info(f"  Outputs will be written to: {expt_dirs.expt_dir}.")
    log.info("Done.\n")

    # Basecall and demultiplex
    # --> Have moved into other subcommands
    # log.info("Demultiplexing...")
    # demuxer = GuppyBarcoder(fastq_dir, kit="rapid96")
    # demuxer.run(
    #     output_dir=expt_dirs.demux_dir,
    #     use_gpu=False,
    #     dry_run=False,
    #     strict=True,
    #     trim_barcodes=True
    # )
    # log.info("Done.\n")

    log.info("Counting FASTQ files...")
    analyse_fastq = ExperimentFASTQCount(expt_dirs, metadata)
    analyse_fastq.run()

    log.info("Mapping to P.f...")
    map_pf = ExperimentMapToPfalciparum(expt_dirs, metadata)
    map_pf.run()

    log.info("Mapping to H.s....")
    map_hs = ExperimentMapUnmappedToHSapiens(expt_dirs, metadata)
    map_hs.run()

    log.info("Running BEDCOV analysis...")
    bedcov = ExperimentBedCoverage(expt_dirs, metadata, regions)
    bedcov.run()

    log.info("Done.")

