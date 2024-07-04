import os
import logging
from savanna.util.dirs import produce_dir


def config_root_logger(log_path: str, verbose=False):
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
    log_dir = produce_dir(os.path.dirname(log_path))
    file_handler = logging.FileHandler(log_path)
    file_formatter = logging.Formatter(FILE_FORMAT, DATE_FORMAT)
    file_handler.setFormatter(file_formatter)
    logging.getLogger().addHandler(file_handler)  # adds to root


# class LoggingFascade:
#     def config_root_logger():
#         pass

#     def get_logger():
#         pass


class LoggingFascade:
    """
    Interface with pythons logging module

    """

    # Formatting defaults
    DATE_FORMAT = "%Y-%m-%d %H:%M"
    STREAM_FORMAT = "%(message)s"
    FILE_FORMAT = "[%(asctime)s][%(levelname)s] %(message)s"

    def __init__(
        self, logger_name: str = "Default", verbose: bool = False, log_path: str = None
    ):
        """
        Instantiate the the default logger, create the console and optionally
        the file handler

        """

        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(logging.DEBUG if verbose else logging.INFO)

        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self._add_console_handler()

        self.log_path = log_path
        if self.log_path is not None:
            self._add_file_handler(log_path)

    def _add_console_handler(self):
        """
        Add a console handler

        """

        console_handler = logging.StreamHandler()
        console_formatter = logging.Formatter(self.STREAM_FORMAT)
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)

    def _add_file_handler(self, log_path: str):
        """
        Add a file handler

        """

        file_handler = logging.FileHandler(log_path)
        file_formatter = logging.Formatter(self.FILE_FORMAT, self.DATE_FORMAT)
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)

    def info(self, message: str):
        self.logger.info(message)

    def debug(self, message: str):
        self.logger.debug(message)

    def warning(self, message: str):
        self.logger.warning(message)

    def error(self, message: str):
        self.logger.error(message)

    def critical(self, message: str):
        self.logger.critical(message)
