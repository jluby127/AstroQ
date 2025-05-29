import os
import logging
from kpfcc import driver

__version__='2.1.0'
_ROOT = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.WARNING,
)
formatter = logging.Formatter('%(levelname)s - %(name)s - %(message)s')
file_handler = logging.FileHandler('my_log_file.log')
file_handler.setFormatter(formatter)

logger.handlers.clear()
logger.addHandler(file_handler)

# eval(f'logging{loglevel}')