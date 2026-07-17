"""Night-plan Plotly figures (ladder, script plan)."""

from astroq.plot._common import *  # noqa: F403

def _min_to_utc_hhmm(night_start, minutes):
    if night_start is None or pd.isna(minutes):
        return ""
    return (night_start + TimeDelta(float(minutes) * 60, format="sec")).isot[11:16]


def _floor_utc_hour(dt):
    return dt.replace(minute=0, second=0, microsecond=0)


def _ladder_utc_ticks(night_start_time, x_min, x_max):
    """Whole-hour UTC tick positions within ``[x_min, x_max]`` (minutes from night start)."""
    start_dt = night_start_time.to_datetime()
    utc_tickvals, utc_ticktext = [], []
    t = _floor_utc_hour(start_dt)
    while True:
        offset_min = (t - start_dt).total_seconds() / 60.0
        if offset_min > x_max + 1e-9:
            break
        if offset_min >= x_min - 1e-9:
            utc_tickvals.append(offset_min)
            utc_ticktext.append(t.strftime("%H:%M"))
        t += timedelta(hours=1)
    return utc_tickvals, utc_ticktext


def _ladder_minute_axis_ticks(x_min, x_max, interval=60):
    """Tick positions for minutes-since-start (0, 60, 120, …), not UTC-aligned."""
    first = 0 if x_min <= 0 else int(np.ceil(x_min / interval)) * interval
    tickvals = []
    v = float(first)
    while v <= x_max + 1e-9:
        if v >= x_min - 1e-9:
            tickvals.append(v)
        v += interval
    return tickvals, [str(int(round(v))) for v in tickvals]


def _add_ladder_night_boundary(fig, x, utc_hhmm, *, side):
    """Vertical marker at night start or end with UTC time label."""
    fig.add_shape(
        type="line",
        x0=x,
        x1=x,
        y0=0,
        y1=1,
        xref="x",
        yref="paper",
        line=dict(color="black", width=1.5),
        layer="above",
    )
    if side == "start":
        label = f"Night start {utc_hhmm} (UT)"
        xanchor = "left"
    else:
        label = f"Night end {utc_hhmm} (UT)"
        xanchor = "right"
    fig.add_annotation(
        x=x,
        y=0.99,
        xref="x",
        yref="paper",
        text=f"<b>{label}</b>",
        showarrow=False,
        yanchor="top",
        xanchor=xanchor,
        font=dict(size=11, color="black"),
    )


def _night_duration_min(model):
    """Night length in minutes (works on live and HDF5-reloaded models)."""
    dur = getattr(model, "dur_min", None)
    if dur is not None:
        return float(dur)
    stats = getattr(model, "stats", None) or {}
    if "dur_min" in stats:
        return float(stats["dur_min"])
    night_start = getattr(model, "night_start", None)
    night_end = getattr(model, "night_end", None)
    if night_start is not None and night_end is not None:
        return (night_end.jd - night_start.jd) * 24 * 60
    return 600.0


def _night_timeline_segments(model):
    """Chronological (idle, visit, slew) segments from 0 .. dur_min."""
    scheduled = model.schedule[~model.schedule["is_anchor"]]
    scheduled = scheduled[scheduled["scheduled"]].sort_values("order")
    dur = _night_duration_min(model)
    cursor = 0.0
    segments = []
    for _, row in scheduled.iterrows():
        t0, t1 = float(row["t_start"]), float(row["t_end"])
        slew = float(row["t_slew"]) if pd.notna(row["t_slew"]) else 0.0
        if t0 > cursor + 1e-9:
            segments.append(("idle", cursor, t0))
        segments.append(("visit", t0, t1))
        if slew > 0:
            segments.append(("slew", t1, t1 + slew))
        cursor = t1 + slew
    if cursor < dur - 1e-9:
        segments.append(("idle", cursor, dur))
    return segments


def _night_aggregate_segments(model):
    """Contiguous visit, slew, idle totals (left-aligned) for the aggregate summary row."""
    totals = {"visit": 0.0, "slew": 0.0, "idle": 0.0}
    for kind, x0, x1 in _night_timeline_segments(model):
        totals[kind] += x1 - x0
    cursor = 0.0
    segments = []
    for kind in ("visit", "slew", "idle"):
        width = totals[kind]
        if width <= 1e-9:
            continue
        segments.append((kind, cursor, cursor + width))
        cursor += width
    return segments


_VISIT_COLOR = "lightgreen"
_SLEW_COLOR = "#FFE4B5"
_IDLE_COLOR = "#FFB3B3"
_ACCESS_LINE_COLOR = "lightgreen"
_ACCESS_LINE_WIDTH = 2
_ACCESS_FILL_COLOR = "white"
_LADDER_BG = "#f4f4f4"
_SEGMENT_COLORS = {"visit": _VISIT_COLOR, "slew": _SLEW_COLOR, "idle": _IDLE_COLOR}


def _add_ladder_required_visit_bar(fig, *, y, x0, visit_len, showlegend=False):
    """Cross-hatched bar for an unscheduled request's required visit duration."""
    x1 = x0 + visit_len
    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y, y, y, y, y],
            mode="lines",
            fill="toself",
            fillcolor="white",
            fillpattern=dict(
                shape="/",
                bgcolor="white",
                fgcolor="#666666",
                fgopacity=0.75,
                size=8,
                solidity=0.4,
            ),
            line=dict(color="#888888", width=1),
            hoverinfo="skip",
            showlegend=showlegend,
            name="Required visit",
            legendgroup="required_visit",
        )
    )


def _synthetic_ladder_row(columns):
    """Blank ladder row with NaN numerics and empty UTC strings."""
    row = {}
    for col in columns:
        if col == "_row_kind":
            continue
        if col.endswith("(UTC)") or col == "Scheduled (UTC)":
            row[col] = ""
        elif col == "is_scheduled":
            row[col] = False
        elif col == "Target":
            row[col] = " "
        else:
            row[col] = np.nan
    return row


def get_ladder(data, tonight_start_time):
    """Produce a plotly figure which illustrates the night plan solution.

    Args:
        data (obj): a TTP data object containing the schedule information

    Returns:
        fig (plotly figure): a plotly figure illustrating the night plan solution.
    """

    from astroq.ttp.plot import schedule_to_ladder_frame, _as_model

    model = _as_model(data)
    orderData = schedule_to_ladder_frame(model)
    if orderData.empty:
        orderData = pd.DataFrame(
            columns=[
                "unique_id",
                "Target",
                "Earliest Start",
                "Latest Finish",
                "Start Exposure",
                "Stop Exposure",
                "Visit Length (min)",
                "Slew to Next (min)",
                "Scheduled (min. from start)",
            ]
        )
    if "Slew to Next (min)" not in orderData.columns:
        orderData["Slew to Next (min)"] = 0.0

    night_start = getattr(model, "night_start", None)
    if night_start is not None and len(orderData):
        orderData["Scheduled (UTC)"] = [
            _min_to_utc_hhmm(night_start, se) if se > 0 else ""
            for se in orderData["Start Exposure"]
        ]
        orderData["Earliest start (UTC)"] = [
            _min_to_utc_hhmm(night_start, t) for t in orderData["Earliest Start"]
        ]
        orderData["Latest finish (UTC)"] = [
            _min_to_utc_hhmm(night_start, t) for t in orderData["Latest Finish"]
        ]
    elif len(orderData):
        orderData["Scheduled (UTC)"] = ""
        orderData["Earliest start (UTC)"] = ""
        orderData["Latest finish (UTC)"] = ""

    on_sky = model.schedule[~model.schedule["is_anchor"]]
    n_unscheduled = int((~on_sky["scheduled"]).sum())

    # reverse so the plot flows top -> bottom with time; after reversal,
    # the lowest indices (bottom of plot) hold the unscheduled block.
    orderData = orderData.iloc[::-1].reset_index(drop=True)
    orderData["_row_kind"] = "target"

    n_before_insert = len(orderData)
    summary_y = None
    aggregate_y = None
    if n_unscheduled > 0 and n_unscheduled < n_before_insert:
        unsched_header = _synthetic_ladder_row(orderData.columns)
        unsched_header["unique_id"] = "__unsched_header__"
        unsched_header["Target"] = "Unscheduled targets"
        unsched_header["_row_kind"] = "section_header"

        aggregate = _synthetic_ladder_row(orderData.columns)
        aggregate["unique_id"] = "__aggregate__"
        aggregate["Target"] = " "
        aggregate["_row_kind"] = "aggregate"

        summary = _synthetic_ladder_row(orderData.columns)
        summary["unique_id"] = "__summary__"
        summary["Target"] = "All scheduled targets"
        summary["_row_kind"] = "summary"

        sched_header = _synthetic_ladder_row(orderData.columns)
        sched_header["unique_id"] = "__sched_header__"
        sched_header["Target"] = "Scheduled targets"
        sched_header["_row_kind"] = "section_header"

        orderData = pd.concat(
            [
                orderData.iloc[:n_unscheduled],
                pd.DataFrame([unsched_header]),
                pd.DataFrame([aggregate]),
                pd.DataFrame([summary]),
                orderData.iloc[n_unscheduled:],
                pd.DataFrame([sched_header]),
            ],
            ignore_index=True,
        )
        aggregate_y = n_unscheduled + 1
        summary_y = n_unscheduled + 2

    # Hide scatter markers on synthetic rows
    mask = orderData["_row_kind"] != "target"
    orderData.loc[mask, "Scheduled (min. from start)"] = np.nan

    # One categorical slot per dataframe row (Target names are not unique).
    orderData["_ladder_y"] = orderData.index.astype(str)

    plot_height = max(400, 40 * len(orderData) + 200)
    categories = orderData["_ladder_y"].tolist()
    y_ticktext = [
        "" if kind in ("section_header", "summary", "aggregate") else target
        for target, kind in zip(orderData["Target"], orderData["_row_kind"])
    ]
    fig = px.scatter(
        orderData,
        x="Scheduled (min. from start)",
        y="_ladder_y",
        title="Night Plan",
        width=800,
        height=plot_height,
    )
    fig.update_traces(
        customdata=np.column_stack(
            [
                orderData["Target"],
                orderData["Scheduled (UTC)"].fillna(""),
                orderData["Earliest Start"],
                orderData["Earliest start (UTC)"].fillna(""),
                orderData["Latest Finish"],
                orderData["Latest finish (UTC)"].fillna(""),
                orderData["Visit Length (min)"],
                orderData["Slew to Next (min)"],
            ]
        ),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "Scheduled (min. from start): %{x:.1f}<br>"
            "Scheduled (UTC): %{customdata[1]}<br>"
            "Earliest start (min. from start): %{customdata[2]:.1f}<br>"
            "Earliest start (UTC): %{customdata[3]}<br>"
            "Latest finish (min. from start): %{customdata[4]:.1f}<br>"
            "Latest finish (UTC): %{customdata[5]}<br>"
            "Visit Length (min): %{customdata[6]:.1f}<br>"
            "Slew to Next (min): %{customdata[7]:.1f}"
            "<extra></extra>"
        ),
        marker=dict(size=0, opacity=0),
    )
    fig.update_layout(
        margin=dict(l=160),
        plot_bgcolor=_LADDER_BG,
        paper_bgcolor="white",
        yaxis_title="",
        yaxis=dict(
            categoryorder="array",
            categoryarray=categories,
            tickmode="array",
            tickvals=categories,
            ticktext=y_ticktext,
        ),
    )
    # x-axis titles/ticks applied after x_max is known (minutes bottom, UTC top)
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor=_VISIT_COLOR,
        line=dict(width=0),
        showlegend=True,
        name="Visit",
    )
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor=_SLEW_COLOR,
        line=dict(width=0),
        showlegend=True,
        name="Slew",
    )
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor=_IDLE_COLOR,
        line=dict(width=0),
        showlegend=True,
        name="Idle",
    )
    fig.add_shape(
        type="rect",
        x0=-100,
        x1=-80,
        y0=-0.5,
        y1=0.5,
        fillcolor=_ACCESS_FILL_COLOR,
        line=dict(color=_ACCESS_LINE_COLOR, width=_ACCESS_LINE_WIDTH),
        showlegend=True,
        name="Accessible",
    )

    new_already_processed = []
    ifixer = 0
    required_visit_bars = []
    for i in range(len(orderData)):
        if orderData["_row_kind"].iloc[i] != "target":
            continue
        if orderData["unique_id"][i] not in new_already_processed:
            indices = [
                k
                for k in range(len(orderData))
                if orderData["unique_id"][k] == orderData["unique_id"][i]
            ]
            for j in range(len(indices)):
                if j == 0:
                    fig.add_shape(
                        type="rect",
                        x0=orderData["Earliest Start"][indices[j]],
                        x1=orderData["Latest Finish"][indices[j]],
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor=_ACCESS_FILL_COLOR,
                        line=dict(color=_ACCESS_LINE_COLOR, width=_ACCESS_LINE_WIDTH),
                        showlegend=False,
                    )
                start_exp = float(orderData["Start Exposure"][indices[j]])
                visit_len = float(orderData["Visit Length (min)"][indices[j]])
                is_scheduled = bool(orderData["is_scheduled"].iloc[indices[j]])
                if is_scheduled and start_exp > 0:
                    fig.add_shape(
                        type="rect",
                        x0=start_exp,
                        x1=start_exp + visit_len,
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor=_VISIT_COLOR,
                        line=dict(width=0),
                    )
                elif not is_scheduled and visit_len > 0:
                    earliest = float(orderData["Earliest Start"][indices[j]])
                    if not np.isnan(earliest):
                        required_visit_bars.append(
                            (
                                orderData["_ladder_y"].iloc[indices[j]],
                                earliest,
                                visit_len,
                            )
                        )
                slew = float(orderData["Slew to Next (min)"][indices[j]])
                if is_scheduled and slew > 0:
                    fig.add_shape(
                        type="rect",
                        x0=orderData["Stop Exposure"][indices[j]],
                        x1=orderData["Stop Exposure"][indices[j]] + slew,
                        y0=i + ifixer - 0.5,
                        y1=i + ifixer + 0.5,
                        fillcolor=_SLEW_COLOR,
                        line=dict(width=0),
                    )
            new_already_processed.append(orderData["unique_id"][i])
        else:
            ifixer -= 1

    for bar_idx, (target, x0, visit_len) in enumerate(required_visit_bars):
        _add_ladder_required_visit_bar(
            fig,
            y=target,
            x0=x0,
            visit_len=visit_len,
            showlegend=(bar_idx == 0),
        )

    if aggregate_y is not None:
        for kind, x0, x1 in _night_aggregate_segments(model):
            fig.add_shape(
                type="rect",
                x0=x0,
                x1=x1,
                y0=aggregate_y - 0.5,
                y1=aggregate_y + 0.5,
                fillcolor=_SEGMENT_COLORS[kind],
                line=dict(width=0),
                showlegend=False,
            )

    if summary_y is not None:
        for kind, x0, x1 in _night_timeline_segments(model):
            fig.add_shape(
                type="rect",
                x0=x0,
                x1=x1,
                y0=summary_y - 0.5,
                y1=summary_y + 0.5,
                fillcolor=_SEGMENT_COLORS[kind],
                line=dict(width=0),
                showlegend=False,
            )

    for idx in range(len(orderData)):
        row = orderData.iloc[idx]
        if row["_row_kind"] not in ("section_header", "summary"):
            continue
        label = row["Target"]
        fig.add_annotation(
            y=idx,
            xref="paper",
            x=0,
            xanchor="right",
            text=f"<b>{label}</b>",
            showarrow=False,
            font=dict(size=13, color="black"),
        )

    night_end = getattr(model, "night_end", None)
    fallback_end = None
    if night_start is not None and night_end is not None:
        fallback_end = (night_end.jd - night_start.jd) * 24 * 60
    elif hasattr(model, "dur_min") and model.dur_min is not None:
        fallback_end = float(model.dur_min)
    elif getattr(model, "stats", None) and "dur_min" in model.stats:
        fallback_end = float(model.stats["dur_min"])
    elif len(orderData) > 0:
        target_rows = orderData[orderData["_row_kind"] == "target"]
        end_times = (
            target_rows["Start Exposure"]
            + target_rows["Visit Length (min)"]
            + target_rows["Slew to Next (min)"]
        )
        fallback_end = float(end_times.max())
    else:
        fallback_end = 600.0

    x_min = 0.0
    x_max = fallback_end
    if tonight_start_time is not None:
        utc_tickvals, utc_ticktext = _ladder_utc_ticks(
            tonight_start_time, x_min, x_max
        )
    else:
        utc_tickvals, utc_ticktext = [], []

    min_tickvals, min_ticktext = _ladder_minute_axis_ticks(x_min, x_max)

    for x_line in utc_tickvals:
        fig.add_shape(
            type="line",
            x0=x_line,
            x1=x_line,
            y0=0,
            y1=1,
            xref="x",
            yref="paper",
            line=dict(color="white", width=1),
            layer="below",
        )

    if tonight_start_time is not None:
        night_end_min = fallback_end
        start_utc = tonight_start_time.isot[11:16]
        if night_end is not None:
            end_utc = night_end.isot[11:16]
        else:
            end_utc = _min_to_utc_hhmm(tonight_start_time, night_end_min)
        _add_ladder_night_boundary(fig, 0.0, start_utc, side="start")
        _add_ladder_night_boundary(fig, night_end_min, end_utc, side="end")

    y_ref = orderData["_ladder_y"].iloc[-1] if len(orderData) else ""
    fig.add_trace(
        go.Scatter(
            x=[x_min, x_max],
            y=[y_ref, y_ref],
            mode="markers",
            marker=dict(size=0.001, opacity=0),
            showlegend=False,
            hoverinfo="skip",
            xaxis="x2",
        )
    )
    fig.update_layout(
        xaxis=dict(
            title="time since start [min]",
            range=[x_min, x_max],
            tickmode="array",
            tickvals=min_tickvals,
            ticktext=min_ticktext,
            showgrid=False,
        ),
        xaxis2=dict(
            title=dict(text="time [UTC]", standoff=0),
            overlaying="x",
            side="top",
            range=[x_min, x_max],
            tickmode="array",
            tickvals=utc_tickvals,
            ticktext=utc_ticktext,
            showgrid=False,
            showline=True,
            mirror=True,
        ),
    )

    return fig


def get_script_plan(night_planner):
    """Generate script plan DataFrame from semester planner and night planner objects.

    This function reads the request_selected.csv file from the semester planner's output directory,
    merges it with the night planner's solution data, and returns a properly formatted DataFrame
    with the same column structure as the original get_script_plan function.

    Args:
        night_planner: NightPlanner object containing solution attribute

    Returns:
        final_df (pd.DataFrame): a formatted observing plan DataFrame
    """

    # Read the request_selected.csv file from the semester planner's output directory
    request_selected_path = os.path.join(
        night_planner.output_directory, "request_selected.csv"
    )

    if not os.path.exists(request_selected_path):
        raise FileNotFoundError(
            f"request_selected.csv not found at {request_selected_path}"
        )

    # Read the request_selected.csv file
    request_selected_df = pd.read_csv(request_selected_path)
    solution = night_planner.solution
    on_sky = solution.schedule[~solution.schedule["is_anchor"]]
    scheduled = on_sky[on_sky["scheduled"]].sort_values("order")

    merged_df = request_selected_df.merge(
        scheduled[["unique_id", "t_start", "t_earliest_start", "t_latest_finish"]],
        on="unique_id",
        how="inner",
    )
    merged_df = merged_df.rename(
        columns={
            "t_start": "Start Exposure",
            "t_earliest_start": "Earliest Start",
            "t_latest_finish": "Latest Finish",
        }
    )

    # Select and reorder only the specific columns requested
    # desired_columns = [
    #     'Start Exposure', 'unique_id', 'target', 'program_code', 'ra', 'dec',
    #     'exptime', 'n_exp', 'n_intra_max', 'tau_intra', 'weather_band_1', 'weather_band_2', 'weather_band_3', 'teff',
    #     'jmag', 'Vmag', 'epoch', 'gaia_id', 'First Available', 'Last Available'
    # ]
    desired_columns = [
        "Earliest Start",
        "Start Exposure",
        "Latest Finish",
        "unique_id",
        "target",
        "program_code",
        "ra",
        "dec",
        "exptime",
        "n_exp",
        "n_intra_max",
        "tau_intra",
        "jmag",
        "Vmag",
    ]

    # Keep only the columns that exist in the merged dataframe
    available_columns = [col for col in desired_columns if col in merged_df.columns]

    # Reorder columns to match the desired structure
    final_df = merged_df[available_columns].copy()

    # Round numeric fields to appropriate decimal places
    if "ra" in final_df.columns:
        # Ensure ra is numeric before rounding, handle 'None' strings
        final_df["ra"] = final_df["ra"].replace("None", pd.NA)
        final_df["ra"] = pd.to_numeric(final_df["ra"], errors="coerce").round(1)

    if "dec" in final_df.columns:
        # Ensure dec is numeric before rounding, handle 'None' strings
        final_df["dec"] = final_df["dec"].replace("None", pd.NA)
        final_df["dec"] = pd.to_numeric(final_df["dec"], errors="coerce").round(1)

    if "jmag" in final_df.columns:
        # Ensure jmag is numeric before rounding, handle 'None' strings
        final_df["jmag"] = final_df["jmag"].replace("None", pd.NA)
        final_df["jmag"] = pd.to_numeric(final_df["jmag"], errors="coerce").round(1)

    if "Vmag" in final_df.columns:
        final_df["Vmag"] = final_df["Vmag"].replace("None", pd.NA)
        final_df["Vmag"] = pd.to_numeric(final_df["Vmag"], errors="coerce").round(1)

    # if 'teff' in final_df.columns:
    #     # Ensure teff is numeric before rounding, handle 'None' strings
    #     final_df['teff'] = final_df['teff'].replace('None', pd.NA)
    #     final_df['teff'] = pd.to_numeric(final_df['teff'], errors='coerce').round(0)

    # Convert time fields from "minutes from start of night" to HST timestamps
    try:
        # Get the night start time from the night planner
        from astroq.nplan import get_nightly_times_from_allocation
        from astropy.time import TimeDelta

        night_start_time, _ = get_nightly_times_from_allocation(
            night_planner.allocation_file, night_planner.current_day
        )

        # Convert the time columns to HST timestamps
        if "Start Exposure" in final_df.columns:
            final_df["Start Exposure"] = final_df["Start Exposure"].apply(
                lambda x: (
                    str(TimeDelta(x * 60, format="sec") + night_start_time)[11:16]
                    if pd.notna(x)
                    else ""
                )
            )

        if "Earliest Start" in final_df.columns:
            final_df["Earliest Start"] = final_df["Earliest Start"].apply(
                lambda x: (
                    str(TimeDelta(x * 60, format="sec") + night_start_time)[11:16]
                    if pd.notna(x)
                    else ""
                )
            )

        if "Latest Finish" in final_df.columns:
            final_df["Latest Finish"] = final_df["Latest Finish"].apply(
                lambda x: (
                    str(TimeDelta(x * 60, format="sec") + night_start_time)[11:16]
                    if pd.notna(x)
                    else ""
                )
            )

    except Exception as e:
        print(f"Warning: Could not convert time fields to HST timestamps: {e}")
        print("Time fields will remain as minutes from start of night")

    # Handle missing values and 'None' strings
    final_df = final_df.replace(["", "NoGaiaName", "None"], pd.NA)

    # Ensure DataFrame is clean and properly structured for DataTables
    final_df = final_df.reset_index(drop=True)
    # Remove duplicate column names if any exist
    final_df = final_df.loc[:, ~final_df.columns.duplicated(keep="first")]
    # Fill NaN values with empty strings to ensure consistent structure
    final_df = final_df.fillna("")
    # Ensure all columns have consistent data types (convert objects to strings)
    for col in final_df.columns:
        if final_df[col].dtype == "object":
            final_df[col] = (
                final_df[col]
                .astype(str)
                .replace("nan", "")
                .replace("None", "")
                .replace("", "")
            )

    return final_df


REQUEST_FRAME_COLUMNS = [
    "target",
    "unique_id",
    "program_code",
    "ra",
    "dec",
    "exptime",
    "n_exp",
    "n_inter_max",
    "tau_inter",
    "n_intra_max",
    "n_intra_min",
    "tau_intra",
    "weather_band_1",
    "weather_band_2",
    "weather_band_3",
    "inactive",
    "comments",
]
BOOLEAN_COLUMNS = {
    "weather_band_1": "Band1",
    "weather_band_2": "Band2",
    "weather_band_3": "Band3",
    "inactive": "Inactive",
}
REQUEST_FRAME_DISPLAY_NAMES = {
    "target": "Target",
    "unique_id": "ID",
    "program_code": "Program",
    "ra": "RA",
    "dec": "Dec",
    "exptime": "ExpTime",
    "comments": "Comments",
}
# Tooltips shown when hovering over column headers.
REQUEST_FRAME_COLUMN_TOOLTIPS = {
    "Star": "Name of the star",
    "ID": "Keck OB database unique ID",
    "Program": "Program Code",
    "RA": "RA in decimal degrees",
    "Dec": "Declination in decimal degrees",
    "ExpTime": "Exposure time in seconds",
    "n_exp": "Number of Exposures per Visit",
    "n_inter_max": "Maximum number of unique nights to observe the star",
    "tau_inter": "The minimum inter-night cadence between unique night observations",
    "n_intra_max": "The desired number of visits to the star in each night it is observed",
    "n_intra_min": "The accepted minimum number of visits to the star in each night it is observed",
    "tau_intra": "The minimum intra-night cadence between visits within a night in hours",
    "Band1": "Allowed to observe in Band1?",
    "Band2": "Allowed to observe in Band2?",
    "Band3": "Allowed to observe in Band3?",
    "Inactive": "Is the star set to inactive?",
    "Comments": "Observer notes (e.g. from Keck star list)",
}


_REQUEST_BAND_COLS = ("Band1", "Band2", "Band3", "Inactive")
_REQUEST_NO_PAD_COLS = (
    "n_inter_max",
    "tau_inter",
    "n_intra_max",
    "n_intra_min",
    "tau_intra",
)
# Numeric columns (post-rename): RA(3), Dec(4), ExpTime(5), n_exp(6),
# n_inter_max(7), tau_inter(8), n_intra_max(9), n_intra_min(10), tau_intra(11).
_REQUEST_NUMERIC_COLS = [3, 4, 5, 6, 7, 8, 9, 10, 11]

