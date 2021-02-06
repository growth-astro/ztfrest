
# Import the relevant packages
import psycopg2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
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

slack_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt")
with open(slack_token, "r") as f:
    access_token = f.read()
web_client = WebClient(token=access_token)

def convert_to_hex(val, delimiter=':', force_sign=False):
    """
    Converts a numerical value into a hexidecimal string

    Parameters:
    ===========
    - val:           float
                     The decimal number to convert to hex.

    - delimiter:     string
                     The delimiter between hours, minutes, and seconds
                     in the output hex string.

    - force_sign:    boolean
                     Include the sign of the string on the output,
                     even if positive? Usually, you will set this to
                     False for RA values and True for DEC

    Returns:
    ========
    A hexadecimal representation of the input value.
    """
    s = np.sign(val)
    s_factor = 1 if s > 0 else -1
    val = np.abs(val)
    degree = int(val)
    minute = int((val  - degree)*60)
    second = (val - degree - minute/60.0)*3600.
    if degree == 0 and s_factor < 0:
        return '-00{2:s}{0:02d}{2:s}{1:.2f}'.format(minute, second, delimiter)
    elif force_sign or s_factor < 0:
        deg_str = '{:+03d}'.format(degree * s_factor)
    else:
        deg_str = '{:02d}'.format(degree * s_factor)
    return '{0:s}{3:s}{1:02d}{3:s}{2:.2f}'.format(deg_str, minute, second, delimiter)

def run_on_event(channel_id, bypass=False,
                 lco_programs=None):

    thread_ts = time.time()

    #response = web_client.chat_postMessage(
    #    channel=channel_id,
    #    text='testing')

    if not bypass:
        try:
            converations = web_client.conversations_list(
                channel=channel_id
            )
            channel_slack_id = channel_id
            for converation in converations:
                for chan in converation.data["channels"]:
                    if chan["name"] == channel_id.replace("#",""):
                        channel_slack_id = chan["id"]
        
            delay_thresh = 15.0
        
            payload = web_client.conversations_history(
                channel=channel_slack_id,
                oldest=thread_ts-delay_thresh
            )
        except:
            return   

        if len(payload["messages"]) == 0:
            return
   
        doGCN, gcn_action = False, 'write'
        for mess in payload["messages"]:
            print(mess)
            message_ts = float(mess["ts"])
            if np.abs(message_ts - thread_ts) > delay_thresh:
                continue
            txt = mess['text']
            txtsplit = list(filter(None,txt.split(" ")))
            if len(txtsplit) == 0: continue
            if txtsplit[0] == "write":
                doTrigger = True
                if len(txtsplit) == 2:
                    name = txtsplit[1]
                elif len(txtsplit) == 3:
                    name = txtsplit[1]
                    trigger_action = txtsplit[2]
                elif len(txtsplit) == 4:
                    name = txtsplit[1]
                    trigger_action = txtsplit[2]
                    filts = txtsplit[3]
            user = mess['user']
        if not doTrigger:
            return
    else:
        user, message_ts = 'test', thread_ts
        # options: gp,rp,ip,zs,Y
        name, gcn_action = 'ZTF21aagwbjr', 'write'

    message = []
    message.append("Hi <@{0}>! You are interested in ztfrest GCNs, right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking into %s for you" % name)

    #web_client.chat_postMessage(
    #    channel=channel_id,
    #    text="\n".join(message)
    #)

    df = pd.read_sql_query(f"SELECT * from candidate WHERE name = '{name}'", con)
    ra, dec = df["ra"].values[0], df["dec"].values[0]
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg) 

    ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

    messages = []
    messages.append(f"We report the discovery of the fast optical transient {name}/????? with the Zwicky Transient Facility (ZTF, Bellm et al. 2019, Graham et al. 2019) at coordinates (J2000, <0.5”):")
    messages.append("\n")
    messages.append(f"RA = %s (%.5fd)" % (ra_hex, ra))

    if dec > 0:
        messages.append(f"Dec = +%s (+%.5fd)" % (dec_hex, dec))
    else:
        messages.append(f"Dec = %s (%.5fd)" % (dec_hex, dec))
  
    messages.append("\n")
    messages.append(f"{name} was found by the new “ZTF Realtime Search and Triggering” (ZTF-ReST) project, which aims at near real-time identification of compelling kilonova candidates in ZTF data using the methods described in Andreoni et al. (2020d), independently of gravitational-wave or gamma-ray triggers.")

    lc = pd.read_sql_query(f"SELECT jd, mag, mag_unc, filter, limmag, programid FROM lightcurve_forced WHERE name = '{name}'", con)
    lc["forced"] = np.ones(len(lc))

    mag = np.array(lc["mag"])
    detections = np.where(mag != 99)[0]
    idlast = detections[-1]
    idfirst = detections[0]
    nondetections = np.where(mag == 99)[0]
    idnon = nondetections[-1]
    for index, row in lc.iterrows():
        tt = Time(row["jd"], format='jd')
        ttstring = ":".join(tt.iso.split(":")[:-1])

        if index == idfirst:
            first_detection_jd = row["jd"]
            first_detection_tt = ttstring
            first_detection_mag = row["mag"]
            first_detection_magerr = row["mag_unc"]
            first_detection_filt = row["filter"]
        elif index == idlast:
            final_detection_jd = row["jd"]
            final_detection_tt = ttstring
            final_detection_mag = row["mag"]
            final_detection_magerr = row["mag_unc"]
            final_detection_filt = row["filter"]
        elif index == idnon:
            final_nondetection_jd = row["jd"]
            final_nondetection_tt = ttstring

    dtmax, dymax, filtmax = np.inf, 0, 'n'
    filtval = {}
    for filt in ["g", "r", "i"]:
        idx = np.where((lc["filter"] == filt) & (mag != 99))[0]

        if len(idx) < 1: continue
        filtval[filt] = lc["mag"][idx[0]]

        if len(idx) < 2: continue
        dy = np.abs(lc["mag"][idx[0]] - lc["mag"][idx[-1]])
        dt = np.abs(lc["jd"][idx[0]] - lc["jd"][idx[-1]])
        
        if dt < dtmax:
            filtmax = filt
            dtmax = dt
            dymax = dy

    colordiff = filtval["g"] - filtval["r"]
    if colordiff > 0:
        color = 'red'
    else:
        color = 'blue'

    messages.append("\n")
    messages.append(f"{name} was first detected on {first_detection_tt} UT, hereafter labelled T_det. It faded by ~%.1f mag in {filtmax}-band in the %.1f days since T_det. The transient was last detected on {final_detection_tt} UT at {final_detection_filt} = %.2f ± %.2f mag. Stringent upper limits constrain the transient onset time to be within ~%.1f day from T_det. The color of the transient appears to be {color}, with g-r~%.1f at T_det. The Galactic extinction on the line of sight is low, with E(B-V)=%.2f mag (Planck Collaboration et al., 2015)." % (dymax, dtmax, final_detection_mag, final_detection_magerr,first_detection_jd-final_nondetection_jd, colordiff,get_dust_info(coord)))
    messages.append("\n")
    messages.append("In the table below, we report photometry obtained on images processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC (Masci et al. 2019).")
    messages.append("\n")
    messages.append("-------------------------------------")
    messages.append("Date (UT) | mag | emag | band")
    messages.append("-------------------------------------")

    for index, row in lc.iterrows():
        tt = Time(row["jd"], format='jd')
        ttstring = ":".join(tt.iso.split(":")[:-1])

        if index == idnon:
            lineform = "%s | > %.2f | - | %s" % (ttstring, row["limmag"], row["filter"])
        elif (not row["mag"] == 99):
            lineform = "%s | %.2f | %.2f | %s" % (ttstring, row["mag"], row["mag_unc"], row["filter"])
        else:
            continue

        messages.append(lineform)

    messages.append("-------------------------------------")
    messages.append("\n")
    messages.append(f"{name} is located off the Galactic plane, with Galactic latitude b_Gal = %.1f deg." % (coord.galactic.b.deg))

    messages.append("\n")
    messages.append("Based on observations obtained with the Samuel Oschin Telescope 48-inch and the 60-inch Telescope at the Palomar Observatory as part of the Zwicky Transient Facility project. ZTF is supported by the National Science Foundation under Grant No. AST-2034437 and a collaboration including Caltech, IPAC, the Weizmann Institute for Science, the Oskar Klein Center at Stockholm University, the University of Maryland, Deutsches Elektronen-Synchrotron and Humboldt University, the TANGO Consortium of Taiwan, the University of Wisconsin at Milwaukee, Trinity College Dublin, Lawrence Livermore National Laboratories, and IN2P3, France. Operations are conducted by COO, IPAC, and UW.\n")

    print("\n".join(messages))

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(messages)
    )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", type=str, default="partnership")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument('--lco-programs', dest='lco_programs',
                        type=str, required=False,
                        default='NOAO2020B-005,TOM2020A-008')

    cfg = parser.parse_args()

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
    db_kn = f"host={info_db['host'][0]} dbname={info_db['dbname'][0]} port={info_db['port'][0]} user={info_db['user'][0]} password={info_db['password'][0]}"
    
    # Connect to the database
    con = psycopg2.connect(db_kn)
    cur = con.cursor()
    print(f"Connected to the '{info_db['dbname'][0]}' database")
    
    # Read the secrets for kowalski access
    secrets = ascii.read('secrets.csv', format='csv')
    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]

    if cfg.debug:
        run_on_event(channel, bypass=True, lco_programs=cfg.lco_programs)
        exit(0)

    while True:
        #try:
        print('Looking for some triggering to do!')
        run_on_event(channel, lco_programs=cfg.lco_programs)
        #except:
        #    pass
        time.sleep(5)
