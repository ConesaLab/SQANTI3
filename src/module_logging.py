import logging
import logging.config
import json
import os

with open('src/data/module_logger_config.json','r') as f:
    MODULE_LOGGING_CONFIG = json.load(f)

os.makedirs(os.path.dirname(MODULE_LOGGING_CONFIG['handlers']['module_file_handler']['filename']), exist_ok=True)
logging.config.dictConfig(MODULE_LOGGING_CONFIG)
qc_logger = logging.getLogger('module_logger')
filter_logger = logging.getLogger('module_logger')
rescue_logger = logging.getLogger('module_logger')
