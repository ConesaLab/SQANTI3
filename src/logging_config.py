import logging
import logging.config
from .config import __version__

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
        
# TODO: Adapt the config so the log is written to the terminal directly
MY_LOGGING_CONFIG = {
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
        'stream_handler': {
            'class': 'logging.StreamHandler',
            'formatter': 'default_formatter',
        },
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
        }
    },
    'loggers': {
        'main_logger': {
            'handlers': ['main_script_handler', 'error_handler'],
            'level': 'INFO',
            'propagate': True
        },
        'QC_logger': {
            'handlers': ['stream_handler'],
            'level': 'INFO',
            'propagate': True
        },
    }
}

logging.config.dictConfig(MY_LOGGING_CONFIG)
main_logger = logging.getLogger('main_logger')
logger = logging.getLogger('QC_logger')

def sqanti_art(main_logger):
    message= f"""
    ░██████╗░██████╗░░█████╗░███╗░░██╗████████╗██╗██████╗░
    ██╔════╝██╔═══██╗██╔══██╗████╗░██║╚══██╔══╝██║╚════██╗
    ╚█████╗░██║██╗██║███████║██╔██╗██║░░░██║░░░██║░█████╔╝
    ░╚═══██╗╚██████╔╝██╔══██║██║╚████║░░░██║░░░██║░╚═══██╗
    ██████╔╝░╚═██╔═╝░██║░░██║██║░╚███║░░░██║░░░██║██████╔╝
    ╚═════╝░░░░╚═╝░░░╚═╝░░╚═╝╚═╝░░╚══╝░░░╚═╝░░░╚═╝╚═════╝░
    Version {__version__}
    """
    print(message)
