# Import the relevant packages
from socket import gethostname
from paramiko import SSHClient
from scp import SCPClient
import psycopg2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u

from slack import RTMClient, WebClient
import numpy as np
import logging
import matplotlib.pyplot as plt
import io
import os
import sys
from astropy.time import Time
import traceback
import time

from penquins import Kowalski

from functions_misc import make_triplet, plot_triplet, get_cutouts
from functions_misc import get_dust_info, plot_lc
from functions_misc import get_xmatch_clu_glade, get_bgal_ebv

def touch(fname):
    try:
        os.utime(fname, None)
    except OSError:
        open(fname, 'a').close()

def run_on_event(channel_id, program_ids=[1,2], bypass=False,
                 only_caltech=False,
                 no_plots=False,
                 start_day=None,
                 end_day=None,
                 outdir=None):

    #response = web_client.chat_postMessage(
    #    channel=channel_id,
    #    text='testing')

    message = []
    message.append("Hi {0}! You are interested in ztfrest scanning, right? Let me get right on that for you.".format(os.environ["USER"]))

    print("\n".join(message))

    # Start with full list, hard rejects will be removed HERE
    scoring_df = pd.read_sql_query("SELECT name FROM candidate where hard_reject is NULL", con)

    # Define the thresholds
    thresh = {'rise': {'g': -1.0, 'r': -1., 'i': -0.5, 'sign_select': '<', 'sign_reject': '>'},
              'fade': {'g': 0.30, 'r': 0.30, 'i': 0.30, 'sign_select': '>', 'sign_reject': '<'}
              }

    # Define the filters (list)
    list_filters = ['g', 'r', 'i']

    # Rise, fade, or both? (list, e.g. ['rise', 'fade'])
    list_rise_fade = ['fade']

    scores = {'rise_select': 5,
              'rise_pen': 0,
              'fade_select': 10,
              'fade_pen': 0,
             }

    message = []
    #message.append(f"Working with {len(scoring_df)} candidates")
    message.append(f"Considering the following filter(s): {list_filters}")
    message.append("---")
    message.append("Selected thresholds:")
    for k1 in thresh.keys():
        for k2 in thresh[k1].keys():
            message.append(f"{k1} {k2}: {thresh[k1][k2]}")
    message.append("---")
    message.append("Scoring points:")
    for k in scores.keys():
        message.append(f"{k}: {scores[k]}")

    # Filter for rise and fade for ALERTS
    for rf in list_rise_fade:
        for f in list_filters:
            rf_filt = pd.read_sql_query(f"SELECT name FROM candidate WHERE index_{rf}_{f} {thresh[rf]['sign_select']} {thresh[rf][f]}",con).values
            #message.append(f"{rf}_{f}_filt_alerts: {len(rf_filt)}" )
    
            # Assign points if condition is met, otherwise 0
            scoring_df[f'{rf}_{f}_filt_alerts'] = [scores[f"{rf}_select"] if name in rf_filt else 0 for name in scoring_df['name']]
    
    # Filter for rise and fade for FORCED PHOTOMETRY
    for rf in list_rise_fade:
        for f in list_filters:
            rf_filt = pd.read_sql_query(f"SELECT name FROM candidate WHERE index_{rf}_forced_{f} {thresh[rf]['sign_select']} {thresh[rf][f]}",con).values
            #message.append(f"{rf}_{f}_filt_forced: {len(rf_filt)}" )
    
            # Assign points if condition is met, otherwise 0
            scoring_df[f'{rf}_{f}_filt_forced'] = [scores[f"{rf}_select"] if name in rf_filt else 0 for name in scoring_df['name']]
    
    # Filter for rise and fade for STACKED FORCED PHOTOMETRY
    for rf in list_rise_fade:
        for f in list_filters:
            rf_filt = pd.read_sql_query(f"SELECT name FROM candidate WHERE index_{rf}_stack_{f} {thresh[rf]['sign_select']} {thresh[rf][f]}",con).values
            #message.append(f"{rf}_{f}_filt_stack: {len(rf_filt)}" )
    
            # Assign points if condition is met, otherwise 0
            scoring_df[f'{rf}_{f}_filt_stack'] = [scores[f"{rf}_select"] if name in rf_filt else 0 for name in scoring_df['name']]

    # ### Penalize slow rise or fade (alerts, forced phot, stacked forced phot)
    # Penalize slow rise and fade for ALERTS
    for rf in list_rise_fade:
        for f in list_filters:
            rf_filt = pd.read_sql_query(f"SELECT name FROM candidate WHERE index_{rf}_{f} {thresh[rf]['sign_reject']} {thresh[rf][f]}",con).values
            #message.append(f"{rf}_{f}_pen_alerts: {len(rf_filt)}" )
    
            # Assign points if condition is met, otherwise 0
            scoring_df[f'{rf}_{f}_pen_alerts'] = [scores[f"{rf}_pen"] if name in rf_filt else 0 for name in scoring_df['name']]
    
    # Penalize slow rise and fade for FORCED PHOTOMETRY
    for rf in list_rise_fade:
        for f in list_filters:
            rf_filt = pd.read_sql_query(f"SELECT name FROM candidate WHERE index_{rf}_forced_{f} {thresh[rf]['sign_reject']} {thresh[rf][f]}",con).values
            #message.append(f"{rf}_{f}_pen_forced: {len(rf_filt)}" )
    
            # Assign points if condition is met, otherwise 0
            scoring_df[f'{rf}_{f}_pen_forced'] = [scores[f"{rf}_pen"] if name in rf_filt else 0 for name in scoring_df['name']]

    # Penalize slow rise and fade for STACKED FORCED PHOTOMETRY
    for rf in list_rise_fade:
        for f in list_filters:
            rf_filt = pd.read_sql_query(f"SELECT name FROM candidate WHERE index_{rf}_stack_{f} {thresh[rf]['sign_reject']} {thresh[rf][f]}",con).values
            #message.append(f"{rf}_{f}_pen_stack: {len(rf_filt)}" )

            # Assign points if condition is met, otherwise 0
            scoring_df[f'{rf}_{f}_pen_stack'] = [scores[f"{rf}_pen"] if name in rf_filt else 0 for name in scoring_df['name']]

    # Penalize long duration transients - TOTAL
    duration_pen = pd.read_sql_query("SELECT name FROM candidate WHERE duration_tot > 14", con).drop_duplicates('name').values
    #message.append('duration_pen: ' + str(len(duration_pen)))
    
    scoring_df['duration_pen'] = [-100 if name in duration_pen else 0 for name in scoring_df['name']]
    
    # Penalize long duration transients - INDIVIDUAL FILTERS
    duration_dict = {"g": 10, "r": 12, "i": 14}
    for f in ["g", "r", "i"]:
        duration_pen = pd.read_sql_query(f"SELECT name FROM candidate     WHERE duration_{f} > {duration_dict[f]}", con).drop_duplicates('name').values
        #message.append(f'duration_pen_{f}: ' + str(len(duration_pen)))
        scoring_df[f'duration_pen_{f}'] = [-100 if name in duration_pen else 0 for name in scoring_df['name']]
    
    # ##  Crossmatch scoring
    # Have ranges in assigning/penalizing points based on rising rate -> fast +1, slow-1, between 0. Ranges defined by RCF SN results, and/or 2018cow
    # 
    # Points for candidates with more than 6 detections
    # 
    # Point there is a non-PSF LS source with phot_z less than x
    # 
    # Penalize if PS object within 1.5 arcsec
    # 
    # Point is PS source with sgscore less than x and mag brighter than y
    # 
    # ### Present in either CLU or GLADE 
    # NOTE: GLADE crossmatch is not yet implemented, WIP.
    
    # Filter for match in either CLU or GLADE 
    galaxy_match_filt = pd.read_sql_query("SELECT name FROM crossmatch WHERE (clu_dist_kpc < 100) and (clu_distmpc > 10)", con).drop_duplicates('name').values
    #message.append('galaxy_match_filt: ' + str(len(galaxy_match_filt)))
    
    # Now the CLU crossmatching is giving zero points, change as desired. 
    scoring_df['galaxy_match_filt'] = [1 if name in galaxy_match_filt else 0 for name in scoring_df['name']]

    # Create dataframe with the results
    result_df = pd.DataFrame([])
    result_df['name'] = scoring_df['name']
    result_df['sum'] = scoring_df.sum(axis=1)
    #message.append(str(result_df.sort_values(by='sum', ascending=False)[:12]))
 
    # Define a scoring threshold
    score_thresh = 1
    #message.append(f"There are {len(result_df[result_df['sum'] > score_thresh])} candidates above the scoring threshold of {score_thresh} in the database")
    
    print("\n".join(message))

    # Plotting cell…
    list_names = result_df.sort_values(by='sum', ascending=False)[result_df['sum'] > score_thresh]['name']
 
    # Select only recent stuff
    list_recent = pd.read_sql_query(f"SELECT name FROM lightcurve WHERE (jd > {start_day.jd}) AND (jd < {end_day.jd}) ", con).drop_duplicates('name')
    list_names = list(n for n in list_names if n in list(list_recent['name']))

    # What if the list is empty
    if len(list_names) == 0:
        message = []
        message.append(f"No candidates with Time.now().jd - recency_thresh and score > {score_thresh}")
        print("\n".join(message))
        return
    # Reverse-sort in alphabetical order
    list_names.sort(reverse=True)

    # Get CLU and GLADE crossmatch information
    clu, glade = get_xmatch_clu_glade(list_names, con, cur)
    
    # Get coords, Galactic latitude and E(B-V)
    bgal_ebv = get_bgal_ebv(list_names, con, cur)

    # Get the linear decay indexes
    names_str = "'" + "','".join(list(list_names)) + "'"
    indexes = pd.read_sql_query(f"SELECT name, \
index_fade_g, index_fade_r, index_fade_i, \
index_fade_forced_g, index_fade_forced_r, index_fade_forced_i, \
index_fade_stack_g, index_fade_stack_r, index_fade_stack_i \
FROM candidate WHERE name IN ({names_str})", con)
    
    list_out = []
    gal_out = []

    # Get all the program ids
    names_str = "'" + "','".join(list(list_names)) + "'"
    pid_alerts = pd.read_sql_query(f"SELECT name, programid FROM lightcurve \
WHERE name in ({names_str}) and magpsf < 50",con)
    pid_alerts.rename(columns={"magpsf": "mag"})
    pid_forced = pd.read_sql_query(f"SELECT name, programid FROM lightcurve_forced \
WHERE name in ({names_str}) and mag < 50",con)
    pid_all = pd.concat([pid_alerts, pid_forced])

    if len(list_names) == 0:
        message.append("No candidates are present!")
        print("\n".join(message))
        return

    for name in list_names:
        message = []
        message.append(name)

        pid_cand = pid_all[pid_all['name'] == name]
        # Do not show if there are no programid=3 detections
        if only_caltech is True and not pid_cand.empty:
            if 3 in set(pid_cand['programid'].values):
                pass
            else:
                message.append("No Caltech detections, skipping")
                message.append("------")
                message.append(" ")
                print("\n".join(message))
                continue

        clu_crossmatch = clu[clu['name'] == name]
        if clu_crossmatch.empty:
            message.append("No CLU crossmatch")
        else:
            message.append("CLU crossmatch:")
            message.append(str(clu[clu['name'] == name]))
        message.append(f"Coordinates: RA, Dec = {'{:.6f}'.format(float(bgal_ebv[bgal_ebv['name'] == name]['ra']))}, {'{:.5f}'.format(float(bgal_ebv[bgal_ebv['name'] == name]['dec']))}")
        # Fade rates
        ti = indexes[indexes['name'] == name]
        # Replace None with NaN
        ti = ti.fillna(value=np.nan)
        message.append(f"Fade rate alerts: \
g: {'{:.2f}'.format(ti['index_fade_g'].values[0])}, \
r: {'{:.2f}'.format(ti['index_fade_r'].values[0])}, \
i: {'{:.2f}'.format(ti['index_fade_i'].values[0])} mag/d")
        message.append(f"Fade rate forced phot: \
g: {'{:.2f}'.format(ti['index_fade_forced_g'].values[0])}, \
r: {'{:.2f}'.format(ti['index_fade_forced_r'].values[0])}, \
i: {'{:.2f}'.format(ti['index_fade_forced_i'].values[0])} mag/d")
        message.append(f"Fade rate stacked forced phot: \
g: {'{:.2f}'.format(ti['index_fade_stack_g'].values[0])}, \
r: {'{:.2f}'.format(ti['index_fade_stack_r'].values[0])}, \
i: {'{:.2f}'.format(ti['index_fade_stack_i'].values[0])} mag/d")
        message.append(f"Score: {result_df[result_df['name'] == name]['sum'].values[0]}")
        message.append(f"Extinction: E(B-V) = {'{:.2f}'.format(float(bgal_ebv[bgal_ebv['name'] == name]['ebv']))}")
        message.append(f"Galactic latitude: b_Gal = {'{:.2f}'.format(float(bgal_ebv[bgal_ebv['name'] == name]['b_gal']))}")

        list_out.append(name)
        gal_out.append(float(bgal_ebv[bgal_ebv['name'] == name]['b_gal']))
        print("\n".join(message))

        forcedbool = 0
        stackbool = 0
        forced = True
        stack = False
        if forced is True:
            forcedbool = 1
            table = 'lightcurve_forced'
        elif stack is True:
            stackbool = 1
            table = 'lightcurve_stacked'
        table = 'lightcurve_forced'
 
        lc = pd.read_sql_query(f"SELECT jd, mag, mag_unc, filter, limmag, programid FROM {table} WHERE name = '{name}'", con)
        lc["forced"] = np.ones(len(lc))

        alerts = pd.read_sql_query(f"SELECT jd, magpsf, sigmapsf, filter, programid FROM lightcurve WHERE name = '{name}'", con)

        d_alerts = {'jd': alerts['jd'],
                    'mag': alerts['magpsf'],
                    'mag_unc': alerts['sigmapsf'],
                    'filter': alerts['filter'],
                    'limmag': np.ones(len(alerts))*99.0,
                    'forced': np.zeros(len(alerts)),
                    'programid': alerts['programid']}
        lc_alerts = pd.DataFrame(data=d_alerts)

        jd_forced = list(lc["jd"])
        for j in set(alerts['jd'].values):
            # If the time difference between the alert and any 
            # forced phot is <5s, remove the alert
            if not (lc_alerts.empty) and not (lc.empty):
                if np.min(np.abs(lc['jd'].values - j)) < 5./60/60/24.:
                    lc_alerts.drop(lc_alerts[lc_alerts['jd'] == j].index,
                                   inplace=True)
            lc = lc.append([lc_alerts], ignore_index=True)
        idx = np.where(np.in1d(lc['programid'],program_ids))[0]
        lc = lc.iloc[idx]
        lc.to_csv(f"{outdir}/lc_{name}_forced{forcedbool}_stacked{stackbool}.csv")

    message = []
    message.append(f"Found {len(list_out)} candidates")
    message.append(str(list_out))

    print("\n".join(message))

def doscp(path_local, path_remote, secrets):
    """Copy the output to schoty"""
    ssh = SSHClient()
    ssh.load_system_host_keys()

    # Fetch connection info from secrets
    hostname = secrets['schoty_host'][0]
    username = secrets['schoty_user'][0]
    password = secrets['schoty_pwd'][0]

    ssh.connect(hostname=hostname,
    username=username,
    password=password)

    # SCPCLient takes a paramiko transport as its only argument
    scp = SCPClient(ssh.get_transport())

    scp.put(path_local, path_remote, recursive=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default="candidates")
    parser.add_argument("-c", "--channel", type=str, default="partnership")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-oc", "--onlycaltech", action="store_true", default=False)
    parser.add_argument("-np", "--noplots", action="store_true", default=False)
    parser.add_argument("-s", "--start_day", type=str, default=(Time.now()-7*u.day).isot.split("T")[0])
    parser.add_argument("-e", "--end_day", type=str, default=Time.now().isot.split("T")[0])
    parser.add_argument("-scp", "--scp", action="store_true", default=False, help="scp to schoty")

    cfg = parser.parse_args()

    baseoutdir = os.path.join(cfg.outdir, cfg.channel)
    if not os.path.isdir(baseoutdir):
        os.makedirs(baseoutdir)

    if cfg.channel == 'partnership':
        channel = 'C01B3ME3GEQ'
        program_ids = [1,2]
    elif cfg.channel == 'caltech':
        channel = 'C01AUCVUKTP'
        program_ids = [1,2,3]
    elif cfg.channel == 'test':
        channel = 'G01A2AUV8Q2'
        program_ids = [1,2,3]
    else:
        print('Sorry, I do not know that channel...')
        exit(0)

    # Read the database secrets
    info = ascii.read('./db_access.csv', format='csv')
    info_db = info[info['db'] == 'db_kn_rt_user']
    if gethostname() == 'usnik':
        host = 'localhost'
    else:
        host = info_db['host'][0]
    db_kn = f"host={host} dbname={info_db['dbname'][0]} port={info_db['port'][0]} user={info_db['user'][0]} password={info_db['password'][0]}"
    
    # Connect to the database
    con = psycopg2.connect(db_kn)
    cur = con.cursor()
    print(f"Connected to the '{info_db['dbname'][0]}' database")
    
    # Read the secrets for kowalski access
    secrets = ascii.read('secrets.csv', format='csv')
    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]

    start_day = Time(cfg.start_day + ' 00:00:00.000',format='iso').jd
    end_day = Time(cfg.end_day + ' 00:00:00.000',format='iso').jd

    outdir = os.path.join(baseoutdir, '%d-%d' % (start_day, end_day)) 
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    run_on_event(channel, program_ids=program_ids, bypass=True,
                 only_caltech=cfg.onlycaltech,
                 no_plots=cfg.noplots,
                 start_day=Time(start_day, format='jd'),
                 end_day=Time(end_day,format='jd'),
                 outdir=outdir)
    if cfg.scp is True:
        path_schoty = "/scr2/ztfrest/ZTF/ztfrest/"
        path_remote = os.path.join(path_schoty, baseoutdir) 
        doscp(outdir, path_remote, secrets)
