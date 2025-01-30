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
        },
        'error_qc_handler': {
            'class': 'logging.FileHandler',
            'level': 'WARNING',
            'filename': '%(output_dir)s',
            'formatter': 'error_formatter'
        }
    },
    'loggers': {
        'main_logger': {
            'handlers': ['main_script_handler', 'error_handler'],
            'level': 'INFO',
            'propagate': True
        },
        'qc_logger': {
            'handlers': ['stream_handler', 'error_qc_handler'],
            'level': 'INFO',
            'propagate': True
        },
    }
}

logging.config.dictConfig(MY_LOGGING_CONFIG)
main_logger = logging.getLogger('main_logger')
qc_logger = logging.getLogger('qc_logger')
def setup_logger(output_dir,module,config=MY_LOGGING_CONFIG):
    log_file= f'{output_dir}/log/{module}_module.err'
    config['handlers']['error_qc_handler']['filename']= log_file
    logging.config.dictConfig(config)
    return logging.getLogger('module_logger')

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
    print(message)

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
    print(message)

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
    print(message)

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
    print(message)