"""
This module provides utilities for setting up and managing loggers with enhanced functionality, 
including colored console output and support for logging to files. It is designed to simplify 
logging configuration and improve readability of log messages in the ColBuilder pipeline.
"""

# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import logging
import sys
import os
import re
from pathlib import Path
from datetime import datetime
import threading
from typing import Optional, Dict, Any, Union
from colorama import init, Fore, Style

init(autoreset=True)

# Define custom log levels for section headers
SECTION = 25
SUBSECTION = 24
TITLE = 26

# Configure the root logger
logging.addLevelName(SECTION, "SECTION")
logging.addLevelName(SUBSECTION, "SUBSECTION")
logging.addLevelName(TITLE, "TITLE")

# Global variables for single log file
_log_file_path = None
_log_file_lock = threading.Lock()
_log_file_announced = False
_console_handler = None
_file_handler = None

# ANSI color code regex pattern for stripping colors
ANSI_ESCAPE = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

class ColoredLogger(logging.Logger):
    """Enhanced logger with colored output and special formatters."""
    
    COLORS = {
        'DEBUG': Fore.CYAN,
        'INFO': Fore.GREEN,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.MAGENTA + Style.BRIGHT,
        'SECTION': Fore.BLUE + Style.BRIGHT,
        'SUBSECTION': Fore.BLUE,
        'TITLE': Fore.MAGENTA + Style.BRIGHT
    }

    def __init__(self, name, level=logging.NOTSET):
        super().__init__(name, level)

    def title(self, title):
        """Display a prominent title in logs."""
        title_length = 80
        border = "*" * (title_length + 4)
        self.log(TITLE, f"{border}")
        self.log(TITLE, f"* {title} *")
        self.log(TITLE, f"{border}")

    def section(self, title):
        """Display a section heading in logs."""
        separator = "=" * 80
        self.log(SECTION, f"{separator}")
        self.log(SECTION, f"{title}")
        self.log(SECTION, f"{separator}")

    def subsection(self, title):
        """Display a subsection heading in logs."""
        separator = "-" * 80
        self.log(SUBSECTION, f"{title}")
        self.log(SUBSECTION, f"{separator}")

class ConsoleFormatter(logging.Formatter):
    """
    Formatter for console output with colors
    """
    COLORS = {
        'DEBUG': Fore.CYAN,
        'INFO': Fore.GREEN,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.MAGENTA + Style.BRIGHT,
        'SECTION': Fore.BLUE + Style.BRIGHT,
        'SUBSECTION': Fore.MAGENTA,
        'TITLE': Fore.MAGENTA + Style.BRIGHT
    }

    def format(self, record):
        levelname = record.levelname
        msg = record.msg
        
        if sys.stdout.isatty() and levelname in self.COLORS:
            if levelname in ('SECTION', 'TITLE', 'SUBSECTION'):
                record.msg = f"{self.COLORS[levelname]}{msg}{Style.RESET_ALL}"
            else:
                if '%(levelname)s' in self._fmt:
                    record.levelname = f"{self.COLORS[levelname]}{levelname}{Style.RESET_ALL}"
        
        result = super().format(record)
        
        record.levelname = levelname
        record.msg = msg
        
        return result

class FileFormatter(logging.Formatter):
    """
    Formatter for file output without colors
    """
    def format(self, record):
        if isinstance(record.msg, str):
            record.msg = ANSI_ESCAPE.sub('', record.msg)
        
        return super().format(record)

# Initialize global handlers
def initialize_handlers(log_dir: Optional[Path] = None, level: int = logging.INFO, debug: bool = False):
    """
    Initialize global console and file handlers.
    
    Parameters:
        log_dir: Directory for log files
        level: Logging level for console
        debug: Whether to enable debug mode
    """
    global _console_handler, _file_handler, _log_file_path, _log_file_announced
    
    if _console_handler is not None and _file_handler is not None:
        return
    
    _console_handler = logging.StreamHandler(sys.stdout)
    _console_handler.setFormatter(ConsoleFormatter('%(message)s'))
    _console_handler.setLevel(level)
    
    with _log_file_lock:
        if log_dir is not None and _log_file_path is None:
            log_dir.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            _log_file_path = log_dir / f'colbuilder_{timestamp}.log'
            
            try:
                _file_handler = logging.FileHandler(_log_file_path)
                _file_handler.setFormatter(FileFormatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
                _file_handler.setLevel(logging.DEBUG)  
                
                if not _log_file_announced:
                    print(f"Logging to file: {_log_file_path}")
                    _log_file_announced = True
            except (IOError, PermissionError) as e:
                print(f"Could not create log file at {_log_file_path}: {e}")
                print("Continuing with console logging only")
                _file_handler = None

# Global registry of loggers to avoid creating duplicates
_loggers = {}

def setup_logger(
    name: str,
    level: int = logging.INFO,
    log_dir: Optional[Path] = None,
    debug: bool = False,
) -> ColoredLogger:
    """
    Set up and configure a logger with console and file output.
    Uses a single log file for all loggers in the application.

    Parameters:
        name: Name of the logger
        level: Logging level (default: INFO)
        log_dir: Directory for log files
        debug: Enable debug mode

    Returns:
        Configured logger instance
    """
    global _console_handler, _file_handler
    
    if name in _loggers:
        return _loggers[name]
    
    initialize_handlers(log_dir, level, debug)
    
    logging.setLoggerClass(ColoredLogger)
    logger = logging.getLogger(name)
    
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    if debug or os.environ.get('COLBUILDER_DEBUG', '0') == '1':
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(level)
    
    logger.propagate = False
    
    if _console_handler is not None:
        logger.addHandler(_console_handler)
    
    if _file_handler is not None:
        logger.addHandler(_file_handler)
    
    _loggers[name] = logger
    return logger

def initialize_root_logger(debug: bool = False, log_dir: Optional[Path] = None) -> logging.Logger:
    """
    Initialize the root logger for the application.
    Should be called once at the start of the application.
    
    Parameters:
        debug: Enable debug mode
        log_dir: Directory for log files
        
    Returns:
        Root logger instance
    """
    if debug:
        os.environ['COLBUILDER_DEBUG'] = '1'
    
    return setup_logger("colbuilder", debug=debug, log_dir=log_dir)

def get_system_info() -> Dict[str, Any]:
    """
    Collect system information for debugging purposes.
    
    Returns:
        Dictionary with system information
    """
    import platform
    import sys
    
    return {
        "platform": platform.platform(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "cwd": os.getcwd(),
        "environment_variables": {
            k: v for k, v in os.environ.items() 
            if k.startswith(('COLBUILDER_', 'PYTHONPATH', 'PATH'))
        }
    }

def log_system_info(logger: logging.Logger) -> None:
    """
    Log system information for debugging.
    
    Parameters:
        logger: Logger instance to use
    """
    info = get_system_info()
    logger.debug("System Information:")
    for key, value in info.items():
        if isinstance(value, dict):
            logger.debug(f"{key}:")
            for k, v in value.items():
                logger.debug(f"  {k}: {v}")
        else:
            logger.debug(f"{key}: {value}")

def log_exception(logger: logging.Logger, exception: Exception) -> None:
    """
    Log an exception with traceback.
    
    Parameters:
        logger: Logger instance
        exception: Exception to log
    """
    import traceback
    logger.error(f"Exception: {type(exception).__name__}: {str(exception)}")
    tb_lines = traceback.format_exception(
        type(exception), 
        exception, 
        exception.__traceback__
    )
    for line in tb_lines:
        logger.debug(line.rstrip())

def reset_logging():
    """
    Reset the logging system. Useful for testing.
    """
    global _log_file_path, _console_handler, _file_handler, _loggers, _log_file_announced
    
    if _console_handler:
        _console_handler.close()
    if _file_handler:
        _file_handler.close()
    
    _log_file_path = None
    _console_handler = None
    _file_handler = None
    _loggers = {}
    _log_file_announced = False
    
    root = logging.getLogger()
    while root.handlers:
        root.removeHandler(root.handlers[0])