"""AstroQ plotting package: semester figures, night figures, and HTML tables."""

from astroq.plot._common import (
    gray,
    clear,
    labelsize,
    hours_per_night,
    _football_cache_dir,
    programs_ledger_for_plot,
)
from astroq.plot.context import (
    PlotData,
    PlotSelection,
    RequestView,
    build_plot_data,
)
from astroq.plot.semester import (
    get_birdseye,
    get_cof,
    get_completion_histogram_by_weight,
    get_completion_vs_target_name,
    get_football,
    get_rawobs,
    get_request_frame,
    get_tau_inter_line,
    get_timebar,
    get_timebar_by_program,
)
from astroq.ttp.plot import get_ladder
from astroq.plot.night import get_script_plan
from astroq.plot.tables import (
    dataframe_to_html,
    nightplan_table_to_html,
    request_frame_to_html,
)

__all__ = [
    "PlotData",
    "PlotSelection",
    "RequestView",
    "build_plot_data",
    "get_birdseye",
    "get_cof",
    "get_completion_histogram_by_weight",
    "get_completion_vs_target_name",
    "get_football",
    "get_ladder",
    "get_rawobs",
    "get_request_frame",
    "get_script_plan",
    "get_tau_inter_line",
    "get_timebar",
    "get_timebar_by_program",
    "dataframe_to_html",
    "nightplan_table_to_html",
    "request_frame_to_html",
    "_football_cache_dir",
    "programs_ledger_for_plot",
]
