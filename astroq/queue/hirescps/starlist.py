"""HIRES-CPS nightly starlist writer.

``write_starlist`` is invoked via :meth:`HIRESCPS.write_starlist` from
:func:`astroq.driver.plan_night` to emit the per-night MAGIQ-format script
(``script_<date>_<version>.txt``) consumed by the Keck observer.
"""

import os

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u


def write_starlist(
    frame,
    schedule,
    night_start_time,
    filler_stars,
    current_day,
    outputdir,
    version="nominal",
    all_active_requests=None,
):
    """
    Generate the nightly script in the correct format.

    Args:
        frame (dataframe): the request.csv in dataframe format for just the targets that were selected to be observed tonight
        schedule (dataframe): TTP ``schedule`` DataFrame (parallel to ``nodes``)
        night_start_time (astropy time object): Beginning of observing interval
        filler_stars (array): star names of the stars added in the bonus round
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file
        version (str): a tag for thescript (e.g. nominal, slowdown, backups, etc)
        all_active_requests (pd.DataFrame | None): full active request frame for
            the semester. When provided, a ``BACKUPS`` section is appended
            listing every active request along with an ``obs=N/M`` token
            (past nights observed / requested) read off the
            ``past_nights_observed`` column.

    Returns:
        list[str]: the lines written to the script file.
    """
    frame["target"] = frame["target"].astype(str)

    on_sky = schedule[~schedule["is_anchor"]]
    scheduled_df = on_sky[on_sky["scheduled"]].sort_values("order")
    extras_df = on_sky[~on_sky["scheduled"]]

    total_exptime = 0
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(
        outputdir, "script_{}_{}.txt".format(current_day, version)
    )

    lines = []
    for _, srow in scheduled_df.iterrows():
        uid = str(srow["unique_id"])
        filler_flag = uid in filler_stars
        row = frame.loc[frame["unique_id"] == uid]
        row.reset_index(inplace=True)
        total_exptime += float(row["exptime"].iloc[0])

        start_exposure_hst = str(
            TimeDelta(srow["t_start"] * 60, format="sec") + night_start_time
        )[11:16]
        first_available_hst = str(
            TimeDelta(srow["t_early"] * 60, format="sec") + night_start_time
        )[11:16]
        last_available_hst = str(
            TimeDelta(srow["t_late"] * 60, format="sec") + night_start_time
        )[11:16]
        lines.append(
            format_hires_row(
                row,
                start_exposure_hst,
                first_available_hst,
                last_available_hst,
                current_day,
                filler_flag=filler_flag,
            )
        )

    lines.append("")
    lines.append("X" * 45 + "EXTRAS" + "X" * 45)
    lines.append("")

    for _, erow in extras_df.iterrows():
        uid = str(erow["unique_id"])
        filler_flag = uid in filler_stars
        row = frame.loc[frame["unique_id"] == uid]
        row.reset_index(inplace=True)
        first_available_hst = str(
            TimeDelta(erow["t_early"] * 60, format="sec") + night_start_time
        )[11:16]
        last_available_hst = str(
            TimeDelta(erow["t_late"] * 60, format="sec") + night_start_time
        )[11:16]
        lines.append(
            format_hires_row(
                row,
                "24:00",
                first_available_hst,
                last_available_hst,
                current_day,
                filler_flag,
                True,
            )
        )

    if all_active_requests is not None:
        backup_df = all_active_requests.copy()
        backup_df["_ra_float"] = pd.to_numeric(backup_df["ra"], errors="coerce")
        backup_df["_vmag_float"] = pd.to_numeric(
            backup_df.get("Vmag", pd.Series(dtype=float)), errors="coerce"
        )
        backup_df = backup_df.sort_values("_ra_float", kind="mergesort")

        def emit_block(header, sub_df):
            lines.append("")
            lines.append(header)
            lines.append("")
            for _, req_row in sub_df.iterrows():
                uid = req_row["unique_id"]
                n_done_raw = req_row.get("past_nights_observed", 0)
                n_done = int(n_done_raw) if pd.notna(n_done_raw) else 0
                n_req_raw = req_row.get("n_inter_max", 0)
                n_req = int(n_req_raw) if pd.notna(n_req_raw) else 0
                obs_token = f"obs={n_done}/{n_req}"

                row = sub_df.loc[sub_df["unique_id"] == uid].head(1)
                row = row.reset_index()
                lines.append(
                    format_hires_row(
                        row,
                        None,
                        None,
                        None,
                        current_day,
                        filler_flag=False,
                        extra=False,
                        omit_timing=True,
                        obs_token=obs_token,
                    )
                )

        emit_block(
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 2026A - Requests - All XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
            backup_df,
        )
        emit_block(
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 2026A - Requests - V < 8 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
            backup_df[backup_df["_vmag_float"] < 8],
        )

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, "w") as f:
        f.write("\n".join(lines))
    return lines


def format_hires_row(
    row,
    obs_time,
    first_available,
    last_available,
    current_day,
    filler_flag=False,
    extra=False,
    omit_timing=False,
    obs_token=None,
):
    """
    Format request data in the specific way needed for the script (relates to the Keck "Magiq"
    software's data ingestion requirements).

    Args:
        row (dataframe): a single row from the requests sheet dataframe
        obs_time (str): the timestamp of the night to begin the exposure according to the TTP.
                        In format HH:MM in HST timezone
        first_available (str): the timestamp of the night where the star is first accessible.
                                In format HH:MM in HST timezone.
        last_available (str): the timestamp of the night where the star is last accessible.
                                In format HH:MM in HST timezone.
        filler_flag (boolean): True of the target was added in the bonus round
        extra (boolean): is this an "extra" target
        omit_timing (boolean): if True, drop the trailing
            ``<obs_time> <first_available> <last_available>`` triplet (used for
            the BACKUPS section, where targets are not tied to a specific time).
            ``obs_time`` / ``first_available`` / ``last_available`` may be passed
            as ``None`` in this mode.
        obs_token (str | None): optional trailing token (e.g. ``"obs=3/10"``)
            appended at the very end of the line. Used by the BACKUPS section to
            show past-vs-requested observation counts.

    Returns:
        line (str): the properly formatted string to be included in the script file
    """

    equinox = "2000"
    # Treat missing/NaN proper motions as zero so apply_space_motion stays finite.
    pmra_raw = row.get("pmra", pd.Series([0.0])).iloc[0] if "pmra" in row else 0.0
    pmdec_raw = row.get("pmdec", pd.Series([0.0])).iloc[0] if "pmdec" in row else 0.0
    pmra = 0.0 if pd.isna(pmra_raw) else float(pmra_raw)
    pmdec = 0.0 if pd.isna(pmdec_raw) else float(pmdec_raw)

    current_time = Time(current_day)
    coord = SkyCoord(
        ra=row["ra"].iloc[0] * u.deg,
        dec=row["dec"].iloc[0] * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=Time(f"J{equinox}"),
    )
    new_coord = coord.apply_space_motion(new_obstime=current_time)
    updated_ra = new_coord.ra.to_string(
        unit=u.hourangle, sep=" ", pad=True, precision=1
    )
    updated_dec = new_coord.dec.to_string(
        unit=u.deg, sep=" ", pad=True, precision=0
    )
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    # Annotate the script line with the observation epoch when proper motion was actually
    # applied, so the observer can tell the coords are propagated (frame/equinox is still J2000).
    epoch_token = (
        f"epoch={current_time.jyear:.1f}"
        if (pmra != 0.0 or pmdec != 0.0)
        else None
    )

    target_str = str(row["target"].iloc[0])
    namestring = " " * (16 - len(target_str[:16])) + target_str[:16]

    # Handle missing columns with default values
    vmag_val = row.get("Vmag", [15.0])[0] if "Vmag" in row else 15.0

    try:
        vmag_val = float(vmag_val) if vmag_val is not None else 15.0
    except (ValueError, TypeError):
        vmag_val = 25.0

    exposurestring = (
        " " * (4 - len(str(int(row["exptime"].iloc[0]))))
        + str(int(row["exptime"].iloc[0]))
        + "/"
        + str(int(row["maxtime"].iloc[0]))
        + " " * (4 - len(str(int(row["maxtime"].iloc[0]))))
    )

    ofstring = "1of" + str(int(row["n_intra_max"].iloc[0]))

    numstring = str(int(row["n_exp"].iloc[0])) + "x"
    vmagstring = (
        "vmag="
        + str(np.round(float(vmag_val), 1))
        + " " * (4 - len(str(np.round(float(vmag_val), 1))))
    )

    programstring = row["program_code"].iloc[0]
    priostring = row["priority"].iloc[0]
    deckerstring = row["decker"].iloc[0]
    cellstring = row["cell in/out?"].iloc[0]
    exp_meter_thresholdstring = row["exp_meter_threshold"].iloc[0]

    if extra == False:
        timestring2 = str(obs_time)
    else:
        # designate a nonsense time
        timestring2 = "24:00"

    line = (
        namestring
        + " "
        + updated_ra
        + " "
        + updated_dec
        + " "
        + str(equinox)
        + " "
        + vmagstring
        + " "
        + exposurestring
        + " "
        + exp_meter_thresholdstring
        + " "
        + deckerstring
        + " "
        + numstring
        + " "
        + cellstring
        + " "
        + priostring
        + " CC "
        + programstring
    )

    if not omit_timing:
        line += " " + timestring2 + " " + first_available + " " + last_available

    if epoch_token is not None:
        line += " " + epoch_token

    # Handle missing Observing Notes column
    observing_notes = (
        row.get("Observing Notes", [""])[0] if "Observing Notes" in row else ""
    )
    if observing_notes and not pd.isnull(observing_notes):
        line += " " + str(observing_notes)

    if obs_token is not None:
        line += " " + str(obs_token)

    return line
