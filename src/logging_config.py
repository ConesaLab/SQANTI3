import logging
import logging.config
import json
import os
from src.config import __version__

class InfoFilter:
    """Allow only LogRecords whose severity levels are below ERROR."""

    def __call__(self, log):
        if log.levelno != logging.INFO:
            return 1
        else:
            return 0

class OnlyInfo:
    """Only allow info level logs to be printed"""
    
    def __call__(self, log):
        if log.levelno == logging.INFO:
            return 1
        else:
            return 0
        
MAIN_LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'default_formatter': {
            'format': '[%(levelname)s:%(asctime)s] %(message)s'
        },
        'simple_formatter': {
            'format': '%(message)s'
        },
        'error_formatter': {
            'format': '%(levelname)s: %(message)s'
        },
        'process_formatter': {
            'format': '%(levelname)s:\n%(message)s'
        }
    },
    'filters': {
        'info_filter': {
            '()': InfoFilter
        },
        'only_info': {
            '()': OnlyInfo
        }
    },
    'handlers': {
        'main_script_handler': {
            'class': 'logging.StreamHandler',
            'filters': ['only_info'],
            'formatter': 'simple_formatter',
        },
        'error_handler': {
            'class': 'logging.StreamHandler',
            'level': 'DEBUG',
            'filters': ['info_filter'],
            'formatter': 'error_formatter',
        },
        "process_handler": {
            "class": "logging.FileHandler",
            "filename": "hey",
            "delay": True,
            "mode": "w",
            "formatter": "process_formatter",
            "encoding": "utf8"
        }

    },
    'loggers': {
        'main_logger': {
            'handlers': ['main_script_handler', 'error_handler'],
            'level': 'INFO',
            'propagate': True
        },
        'art_logger': {
            'handlers': ['main_script_handler'],
            'level': 'INFO',
            'propagate': False
        },
        'process_logger': {
            'handlers': ['process_handler'],
            'level': 'INFO',
            'propagate': False
        }
    }
}

logging.config.dictConfig(MAIN_LOGGING_CONFIG)
main_logger = logging.getLogger('main_logger')
art_logger = logging.getLogger('art_logger')

def save_module_logger_info(logpath, module, level, json_file):
    """
    Saves the updated logger configuration to a JSON file.
    Args:
        logpath (str): Path to the directory for log files.
        module (str): Name of the module.
        level (str): Logging level (e.g., 'DEBUG', 'INFO').
        json_file (str): Path to the JSON configuration file.
    """
    with open(json_file, 'r') as f:
        config = json.load(f)
    try:
        os.remove(f"{logpath}/{module}_module.log")
    except FileNotFoundError:
        os.makedirs(logpath, exist_ok=True)
        pass
    config['handlers']['module_file_handler']['filename'] = f"{logpath}/{module}_module.log"
    config['loggers']['module_logger']['level'] = level

    with open(json_file, 'w') as f:
        json.dump(config, f, indent=4)
     

def sqanti_art():
    message= f"""
    ░██████╗░██████╗░░█████╗░███╗░░██╗████████╗██╗██████╗░
    ██╔════╝██╔═══██╗██╔══██╗████╗░██║╚══██╔══╝██║╚════██╗
    ╚█████╗░██║██╗██║███████║██╔██╗██║░░░██║░░░██║░█████╔╝
    ░╚═══██╗╚██████╔╝██╔══██║██║╚████║░░░██║░░░██║░╚═══██╗
    ██████╔╝░╚═██╔═╝░██║░░██║██║░╚███║░░░██║░░░██║██████╔╝
    ╚═════╝░░░░╚═╝░░░╚═╝░░╚═╝╚═╝░░╚══╝░░░╚═╝░░░╚═╝╚═════╝░
    Version {__version__}
    """
    return message

def qc_art():
    message= f"""
    =====================
      ░██████╗░░█████╗░
      ██╔═══██╗██╔══██╗
      ██║██╗██║██║░░╚═╝
      ╚██████╔╝██║░░██╗
      ░╚═██╔═╝░╚█████╔╝
      ░░░╚═╝░░░░╚════╝░
    =====================
    """
    return message

def filter_art():
    message= f"""
    ================================================     
      ███████╗██╗██╗░░░░░████████╗███████╗██████╗░
      ██╔════╝██║██║░░░░░╚══██╔══╝██╔════╝██╔══██╗
      █████╗░░██║██║░░░░░░░░██║░░░█████╗░░██████╔╝
      ██╔══╝░░██║██║░░░░░░░░██║░░░██╔══╝░░██╔══██╗
      ██║░░░░░██║███████╗░░░██║░░░███████╗██║░░██║
      ╚═╝░░░░░╚═╝╚══════╝░░░╚═╝░░░╚══════╝╚═╝░░╚═╝
    ================================================
    """
    return message

def rescue_art():
    message= f"""
    =====================================================
      ██████╗░███████╗░██████╗░█████╗░██╗░░░██╗███████╗
      ██╔══██╗██╔════╝██╔════╝██╔══██╗██║░░░██║██╔════╝
      ██████╔╝█████╗░░╚█████╗░██║░░╚═╝██║░░░██║█████╗░░
      ██╔══██╗██╔══╝░░░╚═══██╗██║░░██╗██║░░░██║██╔══╝░░
      ██║░░██║███████╗██████╔╝╚█████╔╝╚██████╔╝███████╗
      ╚═╝░░╚═╝╚══════╝╚═════╝░░╚════╝░░╚═════╝░╚══════╝
    =====================================================
    """
    return message

def get_logger_info(logger):
  print(f"Logger Name: {logger.name}")
  print(f"Level: {logging.getLevelName(logger.level)}")
  print("Handlers:")
  for handler in logger.handlers:
      print(f"  - {handler}")
      print(f"    Formatter: {handler.formatter}")
      print(f"    Level: {logging.getLevelName(handler.level)}")
  if logger.propagate:
      print("Propagates to parent logger")
  else:
      print("Does not propagate to parent logger")  
    
