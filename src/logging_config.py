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
            'filename': 'qc_module.err',
            'formatter': 'error_formatter'
        },
        'process_handler': {
            'class': 'logging.FileHandler',
            'filename': 'none',
            'mode': 'w',
            'formatter': 'process_formatter',
            'encoding': 'utf8'
        }
    },
    'loggers': {
        'main_logger': {
            'handlers': ['main_script_handler', 'error_handler'],
            'level': 'INFO',
            'propagate': True
        },
        'qc_logger': {
            'handlers': ['stream_handler'], # TODO: Find a way to dinamically assing the file handler
            'level': 'DEBUG',
            #TODO: Set log level info as default but changeable from the wrapper
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

logging.config.dictConfig(MY_LOGGING_CONFIG)
main_logger = logging.getLogger('main_logger')
qc_logger = logging.getLogger('qc_logger')
art_logger = logging.getLogger('art_logger')

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
    art_logger(message)

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