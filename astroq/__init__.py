"""
AstroQ: Optimized observation scheduling for astronomical observations.
"""

# Standard library imports
import logging
import os

# Local imports
from astroq import driver

__version__ = '2.1.0'
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')

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