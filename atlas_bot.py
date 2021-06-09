
# Import the relevant packages
import psycopg2
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
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
import io
import os
import sys
from astropy.time import Time
import traceback
import time

from atlas import submit_forced_photometry 

slack_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt")
with open(slack_token, "r") as f:
    access_token = f.read()
web_client = WebClient(token=access_token)

atlas_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".atlasforced_secret_key.txt")
with open(atlas_token, "r") as f:
    atlas_token = f.read().rstrip()

def upload_fig(fig, user, filename, channel_id):
    imgdata = io.BytesIO()
    fig.savefig(imgdata, format='png', dpi=600, transparent=True)
    imgdata.seek(0)
    #wc = WebClient(token=bot_access_token)
    web_client.files_upload(
        file=imgdata.getvalue(),
        filename=filename,
        channels=channel_id,
        text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, filename)
    )
    #fig.close()

def run_on_event(channel_id, bypass=False):

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
   
        doTrigger, trigger_action = False, 'forced'
        for mess in payload["messages"]:
            print(mess)
            message_ts = float(mess["ts"])
            if np.abs(message_ts - thread_ts) > delay_thresh:
                continue
            txt = mess['text']
            txtsplit = list(filter(None,txt.split(" ")))
            if len(txtsplit) == 0: continue
            if txtsplit[0] == "atlas":
                doTrigger = True
                if len(txtsplit) == 2:
                    name = txtsplit[1]
                elif len(txtsplit) == 3:
                    name = txtsplit[1]
                    trigger_action = txtsplit[2]
            user = mess['user']
        if not doTrigger:
            return
    else:
        user, message_ts = 'test', thread_ts
        # options: gp,rp,ip,zs,Y
        name, trigger_action = 'ZTF21abdhqle', 'trigger'

    message = []
    message.append("Hi <@{0}>! You are interested in ATLAS forced photometry, right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking into %s for you" % name)

    print(message)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    df = pd.read_sql_query(f"SELECT * from candidate WHERE name = '{name}'", con)
    ra, dec = df["ra"].values[0], df["dec"].values[0]

    outfile = submit_forced_photometry(name, ra, dec, atlas_token)
    data_out = pd.read_csv(outfile)

    mjd, mag, magerr, uJy,duJy, maglim = data_out['MJD'], data_out['m'], data_out['dm'], data_out['uJy'], data_out['duJy'], data_out['mag5sig']
    filters = data_out['F']
    idx1 = np.where(filters == 'o')[0]
    idx2 = np.where(filters == 'c')[0]
    snr = uJy/duJy
    m_3sig = -2.5*np.log10(3 * duJy) + 23.9
    idx_det = np.where(snr >= 3)[0]
    idx_nondet = np.where(snr < 3)[0]

    timetmp = mjd-mjd[0]
    fig, ax1 = plt.subplots(1, 1)
    idx = np.intersect1d(idx1, idx_det)
    ax1.errorbar(timetmp[idx],mag[idx],magerr[idx],fmt='o',color='orange')
    idx = np.intersect1d(idx2, idx_det)
    ax1.errorbar(timetmp[idx],mag[idx],magerr[idx],fmt='o',color='cyan') 
    idx = np.intersect1d(idx1, idx_nondet)
    ax1.errorbar(timetmp[idx],m_3sig[idx],fmt='v',color='orange')
    idx = np.intersect1d(idx2, idx_nondet)
    ax1.errorbar(timetmp[idx],m_3sig[idx],fmt='v',color='cyan')
    ax1.set_xlabel('Time [days] mjd=%.5f' % mjd[0])
    ax1.set_ylabel('Magnitude [ab]')
    ax1.set_xlim([min(timetmp), max(timetmp)])
    ax1.invert_yaxis()
    upload_fig(fig, user, "atlas.png", channel_id)
    plt.close(fig)

    message = []
    message.append("Complete!")
    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", type=str, default="partnership")
    parser.add_argument("-d", "--debug", action="store_true", default=False)

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
        run_on_event(channel, bypass=True)
        exit(0)

    while True:
        #try:
        print('Looking for some triggering to do!')
        run_on_event(channel)
        #except:
        #    pass
        time.sleep(5)
