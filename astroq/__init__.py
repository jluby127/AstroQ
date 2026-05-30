"""
AstroQ: Optimized observation scheduling for astronomical observations.
"""

# Standard library imports
import logging
import os

# Local imports
from astroq import driver

__version__ = "2.1.0"
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

logger = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler()],
)

logger.setLevel(logging.INFO)

# Gurobi prints solver progress to stdout from its C library directly. The
# ``gurobipy`` Python logger re-emits the same content via Python ``logging``,
# which our root handler then formats with a timestamp -- producing duplicate
# lines. Silence the Python-side copy and keep only the raw Gurobi output.
logging.getLogger("gurobipy").setLevel(logging.WARNING)
logging.getLogger("gurobipy").propagate = False
