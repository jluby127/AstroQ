"""KPF-CC nightly starlist writer.

``write_starlist`` is invoked via :meth:`KPFCC.write_starlist` from
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
    Generate the nightly script in the format required by the Keck "Magiq" software.
    Backwards compatable to pre-KPF-CC observing.

    Args:
        frame (dataframe): the request_frame of just the targets that were selected to be observed tonight
        schedule (dataframe): TTP ``schedule`` DataFrame (parallel to ``nodes``)
        night_start_time (astropy time object): Beginning of observing interval
        filler_stars (array): star names of the stars added in the bonus round
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file
        version (str): a tag for thescript (e.g. nominal, slowdown, backups, etc)

    Returns:
        lines (str): the script file as a string
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
            format_kpf_row(
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
            format_kpf_row(
                row,
                "24:00",
                first_available_hst,
                last_available_hst,
                current_day,
                filler_flag,
                True,
            )
        )

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, "w") as f:
        f.write("\n".join(lines))
    print(
        "Total Open Shutter Time Scheduled: "
        + str(np.round((total_exptime / 3600), 2))
        + " hours"
    )
    return lines


def format_kpf_row(
    row,
    obs_time,
    first_available,
    last_available,
    current_day,
    filler_flag=False,
    extra=False,
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

    Returns:
        line (str): the properly formatted string to be included in the script file
    """

    equinox = "2000"
    pmra = row.get("pmra", pd.Series([0.0])).iloc[0] if "pmra" in row else 0.0
    pmdec = row.get("pmdec", pd.Series([0.0])).iloc[0] if "pmdec" in row else 0.0
    updated_ra, updated_dec = pm_correcter(
        row["ra"].iloc[0], row["dec"].iloc[0], pmra, pmdec, current_day, equinox=equinox
    )
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    target_str = str(row["target"].iloc[0])
    namestring = " " * (16 - len(target_str[:16])) + target_str[:16]

    jmag_val = row.get("jmag", [15.0])[0] if "jmag" in row else 15.0
    gmag_val = row.get("gmag", [15.0])[0] if "gmag" in row else 15.0
    teff_val = row.get("teff", [5000])[0] if "teff" in row else 5000
    gaia_id_val = row.get("gaia_id", ["UNKNOWN"])[0] if "gaia_id" in row else "UNKNOWN"

    # Convert to float safely, with fallback to defaults
    try:
        jmag_val = float(jmag_val) if jmag_val is not None else 15.0
    except (ValueError, TypeError):
        jmag_val = 25.0

    try:
        gmag_val = float(gmag_val) if gmag_val is not None else 15.0
    except (ValueError, TypeError):
        gmag_val = 25.0

    try:
        teff_val = float(teff_val) if teff_val is not None else 5000
    except (ValueError, TypeError):
        teff_val = 0.0

    jmagstring = (
        "jmag="
        + str(np.round(float(jmag_val), 1))
        + " " * (4 - len(str(np.round(float(jmag_val), 1))))
    )
    exposurestring = (
        " " * (4 - len(str(int(row["exptime"].iloc[0]))))
        + str(int(row["exptime"].iloc[0]))
        + "/"
        + str(int(row["exptime"].iloc[0]))
        + " " * (4 - len(str(int(row["exptime"].iloc[0]))))
    )

    ofstring = "1of" + str(int(row["n_intra_max"].iloc[0]))
    scstring = "sc=" + "T"

    numstring = str(int(row["n_exp"].iloc[0])) + "x"
    gmagstring = (
        "gmag="
        + str(np.round(float(gmag_val), 1))
        + " " * (4 - len(str(np.round(float(gmag_val), 1))))
    )
    teffstr = "Teff=" + str(int(teff_val)) + " " * (4 - len(str(int(teff_val))))

    gaiastring = str(gaia_id_val) + " " * (25 - len(str(gaia_id_val)))
    programstring = row["program_code"].iloc[0]

    if filler_flag:
        # All targets added in round 2 bonus round are lower priority
        priostring = "p3"
    else:
        priostring = "p1"

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
        + jmagstring
        + " "
        + exposurestring
        + " "
        + ofstring
        + " "
        + scstring
        + " "
        + numstring
        + " "
        + gmagstring
        + " "
        + teffstr
        + " "
        + gaiastring
        + " CC "
        + priostring
        + " "
        + programstring
        + " "
        + timestring2
        + " "
        + first_available
        + " "
        + last_available
    )

    # Handle missing Observing Notes column
    observing_notes = (
        row.get("Observing Notes", [""])[0] if "Observing Notes" in row else ""
    )
    if observing_notes and not pd.isnull(observing_notes):
        line += " " + str(observing_notes)

    return line


def pm_correcter(ra, dec, pmra, pmdec, current_day, equinox="2000"):
    """
    Update a star's coordinates due to proper motion.

    Args:
        ra (float): RA in degrees
        dec (float): Dec in degrees
        pmra (float): proper motion in RA (mas/yr), including cos(Dec)
        pmdec (float): proper motion in Dec (mas/yr)
        equinox (str): original epoch (e.g. '2000.0')
        current_day (str): date to which to propagate (e.g. '2025-04-30')

    Returns:
        formatted_ra (str), formatted_dec (str): updated coordinates as strings
    """
    start_time = Time(f"J{equinox}")
    current_time = Time(current_day)
    coord = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        pm_ra_cosdec=pmra * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=start_time,
    )
    new_coord = coord.apply_space_motion(new_obstime=current_time)
    formatted_ra = new_coord.ra.to_string(
        unit=u.hourangle, sep=" ", pad=True, precision=1
    )
    formatted_dec = new_coord.dec.to_string(unit=u.deg, sep=" ", pad=True, precision=0)

    return formatted_ra, formatted_dec
