"""
Module for preparing all data specific to from Keck Observatory's custom made Observing Block (OB) database.
This is specific to the KPF-CC program and the observatory's infrastructure as way to power the prep kpfcc command.
New observatories should write their own module to connect to a new "prep <your observatory>" command.
"""

# Standard library imports
import io
import json
import logging
import os
# Third-party imports
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import astroplan as apl
import astropy.coordinates as apy
import requests
import re
from bs4 import BeautifulSoup

# Local imports
from astroq.access import Access

logs = logging.getLogger(__name__)

# specify the google sheet URLs for each 2026B program here. 
PROGRAM_URLS_2026B = ['https://docs.google.com/spreadsheets/d/1cjmWsht6d_Q2OrM5mhDz3mHwcPjQhYzm346Oyn3muM8/edit?usp=sharing']
# specify the allocation and programs manual file for the 2026B program here. 
ALLOCATION_MANUAL_2026B = 'allocation_hires_cps_2026B.csv'
PROGRAMS_MANUAL_2026B = 'programs_hires_cps_2026B.csv'

# Column definitions for requests tab (output: no start/stop)
REQUEST_COLS = [
    'program_code', 'starname', 'unique_id', 'ra', 'dec', 'exptime', 'maxtime', 'n_exp',
    'n_inter_max', 'tau_inter', 'n_intra_max', 'n_intra_min', 'tau_intra',
    'minimum_elevation', 'minimum_moon_separation', 'weather_band_1',
    'weather_band_2', 'weather_band_3', 'gaia_id', 'teff', 'jmag', 'Vmag',
    'pmra', 'pmdec', 'epoch', 'exp_meter_threshold', 'inactive', 'decker', 'cell in/out?', 'priority', 'scriptline'
]

# Columns to read from requests tab (includes start/stop for building customs)
REQUEST_COLS_READ = REQUEST_COLS + ['start', 'stop']

# Column definitions for custom dataframe (built from start/stop on requests)
CUSTOM_COLS = ['unique_id', 'starname', 'start', 'stop']

# Backward compatibility alias
cols = REQUEST_COLS


def _parse_bracket_array(s):
    """
    Parse a string like '[2026-02-01 12:00, 2026-03-01 12:00]' into a list of strings.
    Handles '[]', single element '[2026-02-01 12:00]', and multiple. Returns [] if invalid.
    """
    if s is None or not isinstance(s, str):
        return []
    s = s.strip()
    if not s or len(s) < 2 or s[0] != '[' or s[-1] != ']':
        return []
    inner = s[1:-1].strip()
    if not inner:
        return []
    return [part.strip() for part in inner.split(',') if part.strip()]


def _customs_from_requests_df(req_df):
    """
    Build customs DataFrame from requests DataFrame that has start/stop columns.
    start/stop are strings like '[2026-02-01 12:00, 2026-03-01 12:00]'. Pairs by index.
    """
    if req_df is None or req_df.empty:
        return pd.DataFrame(columns=CUSTOM_COLS)
    for col in ("start", "stop"):
        if col not in req_df.columns:
            raise ValueError(f"requests DataFrame missing required column: {col}")
    rows = []
    for _, r in req_df.iterrows():
        uid = r.get('unique_id', '')
        star = r.get('starname', '')
        starts = _parse_bracket_array(r.get('start'))
        stops = _parse_bracket_array(r.get('stop'))
        n = min(len(starts), len(stops))
        for i in range(n):
            rows.append({'unique_id': uid, 'starname': star, 'start': starts[i], 'stop': stops[i]})
    return pd.DataFrame(rows, columns=CUSTOM_COLS) if rows else pd.DataFrame(columns=CUSTOM_COLS)


def _extract_sheet_id(url):
    """Extract Google Sheet ID from a share/edit or publish URL. Raises ValueError if not found."""
    url = (url or "").strip()
    prefix = "/spreadsheets/d/"
    i = url.find(prefix)
    if i == -1:
        raise ValueError(f"No Google Sheet ID in URL: {url[:80]}...")
    start = i + len(prefix)
    end = len(url)
    for j in range(start, len(url)):
        if url[j] in "/?":
            end = j
            break
    sheet_id = url[start:end].strip()
    if not sheet_id:
        raise ValueError(f"Empty sheet ID in URL: {url[:80]}...")
    return sheet_id


def _parse_export_csv(text, required_cols, skip_rows):
    """
    Parse CSV from Google's export (same format as File > Download > CSV).
    Header row is at skip_rows (0-based). Raises on parse error or missing columns.
    """
    df = pd.read_csv(io.StringIO(text), skiprows=skip_rows)
    df.columns = [str(c).strip() for c in df.columns]
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {sorted(missing)}")
    return df[required_cols].copy()


def _pull_sheet_via_public_csv(sheet_id, skip_rows=0):
    """
    Fetch requests tab via public CSV export (no credentials). Customs are built
    from the start/stop columns on the requests tab. Raises if no valid tab found.
    """
    last_missing_msg = None
    for gid in range(10):
        url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid}"
        try:
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if e.response.status_code in (400, 404):
                continue
            raise
        text = resp.text.strip()
        if not text or text.lstrip().startswith("<!") or "<html" in text[:200].lower():
            continue
        try:
            df = _parse_export_csv(text, REQUEST_COLS_READ, skip_rows)
        except ValueError as e:
            last_missing_msg = str(e)
            continue
        requests_df = df[REQUEST_COLS].copy()
        if 'Vmag' in requests_df.columns:
            requests_df = requests_df.rename(columns={'Vmag': 'gmag'})
        custom_df = _customs_from_requests_df(df)
        return requests_df, custom_df
    msg = (
        f"Sheet {sheet_id}: no tab had required 'requests' columns. "
        f"Tried export gid 0..{gid}."
    )
    if last_missing_msg:
        msg += f" Last tab checked: {last_missing_msg}"
    raise ValueError(msg)


def pull_requests(sheet_urls, credentials_path=None, skip_rows=0):
    """
    Pull request and custom data from a list of Google Sheet URLs.

    Reads only the "requests" tab. That tab must have REQUEST_COLS plus "start" and
    "stop" columns. start/stop hold bracket arrays like "[2026-02-01 12:00, 2026-03-01 12:00]";
    customs are built from those (one custom row per start/stop pair). The returned
    requests DataFrame does not include start/stop.

    Args:
        sheet_urls (list of str): List of Google Sheets URLs (Share link).
        credentials_path (str, optional): Path to service account JSON. If None,
            uses GOOGLE_APPLICATION_CREDENTIALS env var.
        skip_rows (int, optional): Rows to skip before the header row. Use 2 if
            your column names are on row 3 (first two rows are comments). Default 0.

    Returns:
        tuple: (requests_df, custom_df) - requests use REQUEST_COLS; custom_df has
        unique_id, starname, start, stop (one row per start/stop pair from requests).
    """
    path = credentials_path or os.environ.get("GOOGLE_APPLICATION_CREDENTIALS")
    use_public = not path or not os.path.isfile(path)

    if use_public:
        request_dfs = []
        custom_dfs = []
        for url in sheet_urls:
            url = (url or "").strip()
            if not url:
                continue
            sheet_id = _extract_sheet_id(url)
            req_df, custom_df = _pull_sheet_via_public_csv(sheet_id, skip_rows=skip_rows)
            request_dfs.append(req_df)
            custom_dfs.append(custom_df)
        requests_df = pd.concat(request_dfs, ignore_index=True) if request_dfs else pd.DataFrame(columns=REQUEST_COLS)
        custom_df = pd.concat(custom_dfs, ignore_index=True) if custom_dfs else pd.DataFrame(columns=CUSTOM_COLS)
        # Convert ra (HH:MM:SS.ss) and dec (+/-DD:MM:SS.s) from sexagesimal to decimal degrees
        if not requests_df.empty and "ra" in requests_df.columns and "dec" in requests_df.columns:
            c = SkyCoord(ra=requests_df["ra"].astype(str), dec=requests_df["dec"].astype(str), unit=(u.hourangle, u.deg))
            requests_df = requests_df.copy()
            requests_df["ra"] = c.ra.deg
            requests_df["dec"] = c.dec.deg
        return requests_df, custom_df

    try:
        import gspread
        from google.oauth2.service_account import Credentials
    except ImportError as e:
        raise ImportError(
            "pull_requests with credentials requires gspread and google-auth. "
            "Install with: pip install gspread google-auth"
        ) from e

    scopes = [
        "https://www.googleapis.com/auth/spreadsheets.readonly",
        "https://www.googleapis.com/auth/drive.readonly",
    ]
    creds = Credentials.from_service_account_file(path, scopes=scopes)
    client = gspread.authorize(creds)

    request_dfs = []
    custom_dfs = []

    for url in sheet_urls:
        url = (url or "").strip()
        if not url:
            continue
        wb = client.open_by_url(url)
        ws_req = wb.worksheet("requests")
        header_row = skip_rows + 1  # 1-based; skip_rows=2 -> header on row 3
        req_records = ws_req.get_all_records(head=header_row)
        if not req_records:
            raise ValueError(f"Sheet {url[:60]}... 'requests' tab is empty.")
        req_df = pd.DataFrame(req_records)
        req_df.columns = [c.strip() for c in req_df.columns]
        missing_req = set(REQUEST_COLS_READ) - set(req_df.columns)
        if missing_req:
            raise ValueError(
                f"Sheet 'requests' tab missing required columns (need start/stop): {sorted(missing_req)}"
            )
        request_dfs.append(req_df[REQUEST_COLS])
        custom_dfs.append(_customs_from_requests_df(req_df))

    requests_df = pd.concat(request_dfs, ignore_index=True) if request_dfs else pd.DataFrame(columns=REQUEST_COLS)
    custom_df = pd.concat(custom_dfs, ignore_index=True) if custom_dfs else pd.DataFrame(columns=CUSTOM_COLS)
    # Convert ra (HH:MM:SS.ss) and dec (+/-DD:MM:SS.s) from sexagesimal to decimal degrees
    if not requests_df.empty and "ra" in requests_df.columns and "dec" in requests_df.columns:
        c = SkyCoord(ra=requests_df["ra"].astype(str), dec=requests_df["dec"].astype(str), unit=(u.hourangle, u.deg))
        requests_df = requests_df.copy()
        requests_df["ra"] = c.ra.deg
        requests_df["dec"] = c.dec.deg
    if 'Vmag' in requests_df.columns:
        requests_df = requests_df.rename(columns={'Vmag': 'gmag'})
    return requests_df, custom_df


def login_JUMP():
    login_url = 'https://jump.caltech.edu/user/login/'
    s = requests.session()
    s.get(login_url)
    csrftoken = s.cookies['csrftoken']
    # you'll need to add your credentials for username and password
    payload = {'action':'login', 'username':os.environ['KPFCC_JUMP_USERNAME'], 'password':os.environ['KPFCC_JUMP_PASSWORD'],
               'csrfmiddlewaretoken': csrftoken}
    new_login = s.post(login_url, data = payload, headers = dict(Referer = login_url))
    return s

def get_database_explorer(name, path_for_csv, url='https://jump.caltech.edu/explorer/', links=[]):
    # log into JUMP and go to DataBase page
    session = login_JUMP()
    response = session.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    # find the table of queries
    table = soup.find("tbody", attrs={"class":"list"})
    # find the correct row by the query name
    for row in table.find_all("tr"):
        tab = row.find("td", attrs={"class":"name"})
        try:
            x = tab.a.string
        except:
            pass
        else:
            # once it finds the match, save the download link
            if x == name:
                for l in row.find_all("a", href = re.compile('download')):
                    links.append('/'.join(url.split('/')[:3])+l.get('href'))
    # make sure it finds it before saving
    if links != []:
        response = session.get(links[0])
        # pretty sure you can just read this directly into pandas if that's what you want
        open(path_for_csv, 'wb').write(response.content)
    else:
        print(" It isn't finding that query name, try again. (needs to be exact)")
    session.keep_alive = False
    return

def get_hires_past_history(path_to_csv):
    name = 'HIRES_Queue_Past_History_for_CC'
    # comment this line out when playing with synthetic schedules
    get_database_explorer(name, path_to_csv)
    print("All KPF observations pulled from Jump. Saved to csv: " + path_to_csv)
    
    data = pd.read_csv(path_to_csv)
    
    data = data.rename(columns={
        "starname": "target",
        "program_name": "semid"
    })

    # Add dummy columns
    data["exposure_start_time"] = data["timestamp"]
    data["observer"] = ""

    # Create unique_id column
    data["id"] = data["target"]
    
    # make ints
    data["exposure_time"] = data["exposure_time"].astype(int)

    data.to_csv(path_to_csv, index=False)
    print("Data cleaned. Done.")
    

def write_starlist(frame, solution_frame, night_start_time, extras, filler_stars, current_day,
                    outputdir, version='nominal'):
    """
    Generate the nightly script in the correct format.

    Args:
        frame (dataframe): the request.csv in dataframe format for just the targets that were selected to be observed tonight
        solution_frame (dataframe): the solution attribute from the TTP model.plotly object
        night_start_time (astropy time object): Beginning of observing interval
        extras (array): starnames of "extra" stars (those not fit into the script)
        filler_stars (array): star names of the stars added in the bonus round
        current_day (str): today's date in format YYYY-MM-DD
        outputdir (str): the directory to save the script file
        version (str): a tag for thescript (e.g. nominal, slowdown, backups, etc)

    Returns:
        None
    """
    # Cast starname column to strings to ensure proper matching
    frame['starname'] = frame['starname'].astype(str)
    
    # Cast extras star names to strings to ensure proper matching
    extras['Starname'] = extras['Starname'].astype(str) if isinstance(extras, pd.DataFrame) else [str(star) for star in extras['Starname']]
    
    total_exptime = 0
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    script_file = os.path.join(outputdir,'script_{}_{}.txt'.format(current_day, version))

    lines = []
    for i, item in enumerate(solution_frame['Starname']):
        filler_flag = solution_frame['Starname'][i] in filler_stars
        row = frame.loc[frame['unique_id'] == solution_frame['Starname'][i]]
        row.reset_index(inplace=True)
        total_exptime += float(row['exptime'].iloc[0])

        start_exposure_hst = str(TimeDelta(solution_frame['Start Exposure'][i]*60,format='sec') + \
                                                night_start_time)[11:16]
        first_available_hst = str(TimeDelta(solution_frame['First Available'][i]*60,format='sec')+ \
                                                night_start_time)[11:16]
        last_available_hst = str(TimeDelta(solution_frame['Last Available'][i]*60,format='sec') + \
                                                night_start_time)[11:16]
        lines.append(format_hires_row(row, start_exposure_hst, first_available_hst,last_available_hst,
                                        current_day, filler_flag = filler_flag))

    lines.append('')
    lines.append('X' * 45 + 'EXTRAS' + 'X' * 45)
    lines.append('')

    for j in range(len(extras['Starname'])):
        if extras['Starname'][j] in filler_stars:
            filler_flag = True
        else:
            filler_flag = False
        row = frame.loc[frame['unique_id'] == extras['Starname'][j]]
        row.reset_index(inplace=True)

        lines.append(format_hires_row(row, '24:00', extras['First Available'][j],
                    extras['Last Available'][j], current_day, filler_flag, True))

    # add buffer lines to end of file
    lines.append("")
    lines.append("")

    with open(script_file, 'w') as f:
        f.write('\n'.join(lines))
    # logs.info("Total Open Shutter Time Scheduled: " + str(np.round((total_exptime/3600),2)) + " hours")
    return lines

def format_hires_row(row, obs_time, first_available, last_available, current_day,
                    filler_flag = False, extra=False):
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

    equinox = '2000'
    # Handle missing pmra/pmdec columns with default values
    pmra = row.get('pmra', pd.Series([0.0])).iloc[0] if 'pmra' in row else 0.0
    pmdec = row.get('pmdec', pd.Series([0.0])).iloc[0] if 'pmdec' in row else 0.0
    updated_ra, updated_dec = pm_correcter(row['ra'].iloc[0], row['dec'].iloc[0],
                                pmra, pmdec, current_day, equinox=equinox)
    if updated_dec[0] != "-":
        updated_dec = "+" + updated_dec

    starname_str = str(row['starname'].iloc[0])
    namestring = ' '*(16-len(starname_str[:16])) + starname_str[:16]

    # Handle missing columns with default values
    gmag_val = row.get('gmag', [15.0])[0] if 'gmag' in row else 15.0
    
    try:
        gmag_val = float(gmag_val) if gmag_val is not None else 15.0
    except (ValueError, TypeError):
        gmag_val = 25.0
    
    exposurestring = (' '*(4-len(str(int(row['exptime'].iloc[0])))) + \
        str(int(row['exptime'].iloc[0])) + '/' + \
        str(int(row['maxtime'].iloc[0])) + ' '* \
        (4-len(str(int(row['maxtime'].iloc[0])))))

    ofstring = ('1of' + str(int(row['n_intra_max'].iloc[0])))

    numstring = str(int(row['n_exp'].iloc[0])) + "x"
    gmagstring = 'vmag=' + str(np.round(float(gmag_val),1)) + \
                                                ' '*(4-len(str(np.round(float(gmag_val),1))))

    programstring = row['program_code'].iloc[0]
    priostring = row['priority'].iloc[0]
    deckerstring = row['decker'].iloc[0]
    cellstring = row['cell in/out?'].iloc[0]
    exp_meter_thresholdstring = row['exp_meter_threshold'].iloc[0]

    if extra == False:
        timestring2 = str(obs_time)
    else:
        # designate a nonsense time
        timestring2 = "24:00"

    line = (namestring + ' ' + updated_ra + ' ' + updated_dec + ' ' + str(equinox) + ' '
                + gmagstring + ' ' + exposurestring + ' ' + exp_meter_thresholdstring + ' ' + deckerstring +  ' '
                + numstring + ' ' + cellstring + ' '+ priostring + ' CC '+ programstring + ' ' + timestring2 +
                         ' ' + first_available  + ' ' + last_available )

    # Handle missing Observing Notes column
    observing_notes = row.get('Observing Notes', [''])[0] if 'Observing Notes' in row else ''
    if observing_notes and not pd.isnull(observing_notes):
        line += (' ' + str(observing_notes))

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
    start_time = Time(f'J{equinox}')
    current_time = Time(current_day)
    coord = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        pm_ra_cosdec=pmra * u.mas/u.yr,
        pm_dec=pmdec * u.mas/u.yr,
        obstime=start_time
    )
    new_coord = coord.apply_space_motion(new_obstime=current_time)
    formatted_ra = new_coord.ra.to_string(unit=u.hourangle, sep=' ', pad=True, precision=1)
    formatted_dec = new_coord.dec.to_string(unit=u.deg, sep=' ', pad=True, precision=0)

    return formatted_ra, formatted_dec