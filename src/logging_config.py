import logging
import logging.config
import json
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
            "filename": "none",
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
            'propagate': True
        },
        'process_logger': {
            'handlers': ['process_handler'],
            'level': 'INFO',
            'propagate': True
        }
    }
}

logging.config.dictConfig(MAIN_LOGGING_CONFIG)
main_logger = logging.getLogger('main_logger')
art_logger = logging.getLogger('art_logger')

def save_module_logger_info(logpath, module, level, json_file):
    with open(json_file, 'r') as f:
        config = json.load(f)

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