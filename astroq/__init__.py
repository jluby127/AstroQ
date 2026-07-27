"""
AstroQ: Optimized observation scheduling for astronomical observations.
"""

# Standard library imports
import logging
import warnings

# Third-party imports
from astropy.utils.exceptions import AstropyWarning
from erfa import ErfaWarning
from tables.exceptions import DataTypeWarning

# Booleans persisted via h5py (e.g. show_gurobi_output in splan.py)
# are stored as H5T_ENUM, which PyTables does not
# recognize. PyTables scans root attributes on every pd.read_hdf() and emits a
# DataTypeWarning for each unrecognized attribute. The data still round-trips
# correctly via h5py; silence the cosmetic warning here. Installed before any
# astroq submodule import so the filter is in place when h5 files are first read.
warnings.filterwarnings("ignore", category=DataTypeWarning)

# Nightly script generation (hirescps.starlist.format_hires_row,
# kpfcc.starlist.pm_correcter) propagates catalog RA/Dec to the observation date via
# SkyCoord.apply_space_motion without parallax. ERFA's pmsafe then warns
# "distance overridden" once per target; coordinates are still correct.
warnings.filterwarnings("ignore", category=ErfaWarning)

# Alt/az transforms use IERS Earth-orientation tables. When a night falls outside
# the tabulated (or predictive) range, Astropy falls back to the 50-yr mean polar
# motion (~arcsec). That is fine for slot scheduling; silence the cosmetic warning.
warnings.filterwarnings(
    "ignore",
    message=r".*polar motions for times.*IERS data is valid.*",
    category=AstropyWarning,
)

# Prefer a fresh IERS table when online so near-future nights stay in-range.
try:
    from astropy.utils.iers import IERS_Auto

    IERS_Auto.open()
except Exception:
    pass

# Local imports
from astroq import driver  # noqa: E402

__version__ = "2.1.0"

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
