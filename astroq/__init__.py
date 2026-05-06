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
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()],
)

logger.setLevel(logging.INFO)