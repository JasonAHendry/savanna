import logging
from datetime import datetime 
from savanna.util.dirs import produce_dir


def config_root_logger(log_dir: str, verbose=False):
    """
    Configure the root logger

    For this project, this will be the *only* logger that has handlers.
    All other loggers that are generated will propgate their logging
    records to the root logger for them to be handled.

    """

    # Formatting defaults
    DATE_FORMAT = "%Y-%m-%d %H:%M"
    STREAM_FORMAT = "%(message)s"
    FILE_FORMAT = "[%(asctime)s][%(levelname)s] %(message)s"

    # Set an overall level
    # logging.basicConfig(level=logging.INFO) # Incorrect, this will add another handler I don't want.
    logging.getLogger().setLevel(logging.DEBUG if verbose else logging.INFO)

    # Add console handler
    console_handler = logging.StreamHandler()
    console_formatter = logging.Formatter(STREAM_FORMAT)
    console_handler.setFormatter(console_formatter)
    logging.getLogger().addHandler(console_handler)  # adds to root

    # Add file handler
    log_dir = produce_dir(log_dir)
    log_path = f"{log_dir}/{datetime.today().strftime('%Y-%m-%d')}.log"
    file_handler = logging.FileHandler(log_path)
    file_formatter = logging.Formatter(FILE_FORMAT, DATE_FORMAT)
    file_handler.setFormatter(file_formatter)
    logging.getLogger().addHandler(file_handler)  # adds to root

