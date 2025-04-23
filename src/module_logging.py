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