# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import logging
import sys
from pathlib import Path
from typing import Optional
from colorama import init, Fore, Style

init(autoreset=True)

class ColoredLogger(logging.Logger):
    COLORS = {
        'DEBUG': Fore.BLUE,
        'INFO': Fore.GREEN,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.RED + Style.BRIGHT,
    }

    def __init__(self, name, level=logging.NOTSET):
        super().__init__(name, level)

    def _log(self, level, msg, args, exc_info=None, extra=None, stack_info=False):
        if level >= self.level:
            levelname = logging.getLevelName(level)
            if levelname in self.COLORS:
                msg = f"{self.COLORS[levelname]}{msg}{Style.RESET_ALL}"
        super()._log(level, msg, args, exc_info, extra, stack_info)

    def title(self, title):
        title_length = len(title)
        border = "*" * (title_length + 4)
        self.info(f"{Fore.CYAN}{Style.BRIGHT}{border}")
        self.info(f"{Fore.CYAN}{Style.BRIGHT}* {title} *")
        self.info(f"{Fore.CYAN}{Style.BRIGHT}{border}{Style.RESET_ALL}")

    def section(self, title):
        separator = "=" * len(title)
        self.info(f"{Fore.CYAN}{separator}")
        self.info(f"{Fore.CYAN}{title}")
        #self.info(f"{Fore.CYAN}{separator}{Style.RESET_ALL}")

    def subsection(self, title):
        separator = "-" * len(title)
        #self.info(f"{Fore.MAGENTA}{separator}")
        self.info(f"{Fore.MAGENTA}{title}")
        #self.info(f"{Fore.MAGENTA}{separator}{Style.RESET_ALL}")

def setup_logger(
    name: str,
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    format_string: Optional[str] = None,
) -> ColoredLogger:
    """
    Set up and configure a logger.

    Parameters
    ----------
    name : str
        Name of the logger.
    level : int, optional
        Logging level. Default is logging.INFO.
    log_file : Path, optional
        Path to the log file. If not provided, logs will only be printed to console.
    format_string : str, optional
        Custom format string for log messages. If not provided, a default format will be used.

    Returns
    -------
    logging.Logger
        Configured logger instance.

    Examples
    --------
    >>> logger = setup_logger("my_module")
    >>> logger.info("This is an info message")
    >>> logger.error("This is an error message")

    >>> debug_logger = setup_logger("debug_module", level=logging.DEBUG, 
    ...                             log_file=Path("debug.log"))
    >>> debug_logger.debug("This is a debug message")
    """
    logging.setLoggerClass(ColoredLogger)
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not format_string:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    console_format = "%(message)s"
    console_formatter = logging.Formatter(console_format)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    if log_file:
        if not format_string:
            format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        file_formatter = logging.Formatter(format_string)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger

if __name__ == "__main__":
    basic_logger = setup_logger("basic_logger")
    basic_logger.debug("This is a basic debug log")
    basic_logger.info("This is a basic info log")
    basic_logger.warning("This is a basic warning log")
    basic_logger.error("This is a basic error log")
    basic_logger.critical("This is a basic critical log")

    advanced_logger = setup_logger(
        "advanced_logger",
        level=logging.DEBUG,
        log_file=Path("advanced.log"),
        format_string="%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s"
    )
    advanced_logger.debug("This is an advanced debug log")
    advanced_logger.info("This is an advanced info log")
    advanced_logger.warning("This is an advanced warning log")
    advanced_logger.error("This is an advanced error log")
    advanced_logger.critical("This is an advanced critical log")