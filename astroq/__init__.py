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

# Configure logging to capture everything
logging.basicConfig(
    level=logging.INFO,  # Lower level to capture more messages
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('astroq.log'),  # Better filename
        logging.StreamHandler()  # Also show on console
    ]
)

# Set specific logger levels
logger.setLevel(logging.INFO)

# Redirect stdout and stderr to also go to the log file
import sys
from io import StringIO

class TeeLogger:
    def __init__(self, original_stream, log_file):
        self.original_stream = original_stream
        self.log_file = log_file
        self.buffer = StringIO()
    
    def write(self, text):
        self.original_stream.write(text)
        self.log_file.write(text)
        self.log_file.flush()
    
    def flush(self):
        self.original_stream.flush()
        self.log_file.flush()

# Open log file for stdout/stderr redirection
log_file = open('astroq.log', 'a', encoding='utf-8')

# Redirect stdout and stderr to both console and log file
sys.stdout = TeeLogger(sys.stdout, log_file)
sys.stderr = TeeLogger(sys.stderr, log_file)

# eval(f'logging{loglevel}')