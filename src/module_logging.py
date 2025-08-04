import logging
import logging.config
import json
import os

data_path = os.path.join(os.path.dirname(__file__), 'data','module_logger_config.json')

with open(data_path,'r') as f:
    MODULE_LOGGING_CONFIG = json.load(f)

# Create the loggers
logging.config.dictConfig(MODULE_LOGGING_CONFIG)
qc_logger = logging.getLogger('module_logger')
filter_logger = logging.getLogger('module_logger')
rescue_logger = logging.getLogger('module_logger')
rescue_logger.propagate = False

def message(text,logger):
    size = max(len(text)+4, 50)
    line = "-" * size
    centered_text = f"{text:^{size}}"
    logger.info(line)
    logger.info(centered_text)
    logger.info(line)

def update_logger(logger,dir,module,log_lvl):
    """
    Update the logger configuration to use a new log file.
    """
    logPath= os.path.join(dir, 'logs', f"sqanti3_{module}.log")
    # Get the current log file handler
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            # Update the log file path
            handler.baseFilename = os.path.join(logPath)

            logger.setLevel(log_lvl)
            break
    if not os.path.exists(logPath):
        os.makedirs(f"{dir}/logs", exist_ok=True)
    else:
        os.remove(logPath)