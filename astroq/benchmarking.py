"""
Module for preparing the benchmark test as used in Lubin et al. 2025.

"""

import numpy as np
import pandas as pd

# Programs from Lubin et al. 2025; 1-indexed to match paper notation.
# tau_inter: min nights between visits (days)
# exptime: seconds per visit (program 6 overwritten with ns * 300 at runtime)
# n_inter_max: max number of nights observed
# n_intra_max: max visits per night
PROGRAMS = pd.DataFrame(
    [
        (1, 300, 40, 1),    # 1: APF-50
        (3, 600, 20, 1),    # 2: TKS
        (10, 1200, 10, 1),  # 3: bi-weekly
        (1, 300, 8, 5),     # 4: multi-visit
        (1, 300, 40, 1),    # 5: RA-constrained (Kepler field)
        (0, 0, 1, 1),       # 6: single shot
    ],
    columns=["tau_inter", "exptime", "n_inter_max", "n_intra_max"],
    index=pd.RangeIndex(1, 7, name="program_num"),
)


def getDec(maxDec=90, minDec=-20):
    """
    Randomly draw a declination from a cosine distribution between two values.
    Default min/max are chosen for favorable viewing from Hawaii.
    """
    mincosdec = np.cos((90 + maxDec) * (np.pi / 180.0))
    maxcosdec = np.cos((90 + minDec) * (np.pi / 180.0))
    cosdec = np.random.uniform(mincosdec, maxcosdec)
    return (np.arccos(cosdec) * (180.0 / np.pi)) - 90


def build_toy_model_from_paper(ns, hours_per_program=80, seed=24):
    """
    Generate the nominal request set used for performance testing in Lubin et al. 2025.

    Args:
        ns (int): slots needed to complete each single-shot observation in program 6.
        hours_per_program (float): hours allocated to each program. Paper value: 80.
        seed (int): RNG seed. Paper value: 24.

    Returns:
        pandas.DataFrame: the request frame, equivalent to request.csv.
    """
    np.random.seed(seed)

    programs = PROGRAMS.copy()
    programs.loc[6, "exptime"] = ns * 300
    programs["n_stars"] = (
        (hours_per_program * 3600.0)
        / (programs.exptime * programs.n_inter_max * programs.n_intra_max)
    ).round().astype(int)

    # Cross-join programs to per-star rows; preserves program order (1..6).
    df = (
        programs.loc[programs.index.repeat(programs.n_stars)]
        .reset_index()
        .drop(columns="n_stars")
    )

    # Per-star RA/Dec. Non-Kepler stars draw RA from [0, 360) deg with the
    # 12h-18h band excluded (favorable for a B semester).
    exclude_start, exclude_end = 12 * 15, 18 * 15
    p_below_exclusion = exclude_start / (360 - (exclude_end - exclude_start))

    ra = np.empty(len(df))
    dec = np.empty(len(df))
    for i, prog in enumerate(df["program_num"].to_numpy()):
        if prog == 5:
            ra[i] = np.random.uniform(19 * 15, 19.66 * 15)
            dec[i] = np.random.uniform(40, 50)
        else:
            ra[i] = (
                np.random.uniform(0, exclude_start)
                if np.random.random() < p_below_exclusion
                else np.random.uniform(exclude_end, 360)
            )
            dec[i] = getDec()

    df["ra"] = ra
    df["dec"] = dec

    df.insert(0, "starname", [f"Star{i:04d}" for i in range(len(df))])
    df.insert(1, "program_code", "Program" + df["program_num"].astype(str))
    df = df.drop(columns="program_num")

    # Final schema. n_intra_min == n_intra_max except program 4 (5 -> 3).
    df["n_exp"] = 1
    df["n_intra_min"] = np.where(df["n_intra_max"] == 5, df["n_intra_max"] - 2, df["n_intra_max"])
    df["tau_intra"] = np.where(df["n_intra_max"] == 1, 0, 1)
    df["unique_id"] = df["starname"]
    df["minimum_elevation"] = 33.0
    df["minimum_moon_separation"] = 33.0
    df["inactive"] = False

    return df[
        [
            "starname", "program_code", "unique_id",
            "ra", "dec",
            "exptime", "n_exp",
            "n_intra_max", "n_intra_min", "tau_intra",
            "n_inter_max", "tau_inter",
            "minimum_elevation", "minimum_moon_separation",
            "inactive",
        ]
    ]
