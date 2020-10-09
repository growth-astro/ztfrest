
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
   
        doTrigger, trigger_action = False, 'trigger'
        for mess in payload["messages"]:
            print(mess)
            message_ts = float(mess["ts"])
            if np.abs(message_ts - thread_ts) > delay_thresh:
                continue
            txt = mess['text']
            txtsplit = list(filter(None,txt.split(" ")))
            if len(txtsplit) == 0: continue
            if txtsplit[0] == "trigger":
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
        name, trigger_action = 'ZTF20acgnelh', 'trigger'

    message = []
    message.append("Hi <@{0}>! You are interested in ztfrest triggering, right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking into %s for you" % name)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    df = pd.read_sql_query(f"SELECT * from candidate WHERE name = '{name}'", con)
    ra, dec = df["ra"].values[0], df["dec"].values[0]

    alerts = pd.read_sql_query(f"SELECT jd, magpsf, sigmapsf, filter, programid FROM lightcurve WHERE name = '{name}'", con)
    d_alerts = {'jd': alerts['jd'],
                'mag': alerts['magpsf'],
                'mag_unc': alerts['sigmapsf'],
                'filter': alerts['filter'],
                'limmag': np.ones(len(alerts))*99.0,
                'forced': np.zeros(len(alerts)),
                'programid': alerts['programid']}
    lc_alerts = pd.DataFrame(data=d_alerts)
    medmags = {}

    message = []
    for filt in ["g", "r", "i"]:
        idx = np.where(lc_alerts["filter"] == filt)[0]
        medmags[filt] = np.median(lc_alerts["mag"][idx])
        message.append("Median magnitude for %s-band: %.3f" % (filt, medmags[filt]))

    #web_client.chat_postMessage(
    #    channel=channel_id,
    #    text="\n".join(message)
    #)

    message = []
    if trigger_action in ["status", "delete"]:
        print('Checking LCO for existing observations...')

        # LCO sometime over next 2 weeks
        tstart = Time.now()
        tend = Time.now() + TimeDelta(14*u.day)
        tstart = str(tstart.isot).replace("T"," ")
        tend = str(tend.isot).replace("T"," ")

        #Read the secrets
        lco_secrets = ascii.read('../lco/secrets.csv', format = 'csv')
        PROPOSAL_ID = lco_secrets['PROPOSAL_ID'][0]
        API_TOKEN = lco_secrets['API_TOKEN'][0]

        lco_programs = lco_programs.split(",")

        status = {-1: "CANCELED",
                  0: "WINDOW_EXPIRED",
                  1: "PENDING",
                  2: "COMPLETED"}

        from lco import check_observations, delete_observation
        obs = check_observations(API_TOKEN, lco_programs=lco_programs)
        for key in obs.keys():
            if key == name:
                if trigger_action == "status":
                    message.append('Observations of %s...' % name)
                    message.append('Observation ID: %d...' % obs[key]["obsid"])
                    message.append('Completion: %s...' % status[obs[key]["completed"]])
                    message.append('Filters: %s' % (','.join(obs[key]["observations"])))
                    message.append('Program: %s' % obs[key]["program"])
                elif trigger_action == "delete":
                    delete_observation(obs[key]["obsid"], API_TOKEN) 
                    message.append('Deleted Observation ID: %d...' % obs[key]["obsid"]) 

    elif trigger_action == "trigger":

        print('Triggering LCO...')

        # LCO sometime over next 2 weeks
        tstart = Time.now()
        tend = Time.now() + TimeDelta(2*u.day)
        tstart = str(tstart.isot).replace("T"," ")
        tend = str(tend.isot).replace("T"," ")

        #Read the secrets
        lco_secrets = ascii.read('../lco/secrets.csv', format = 'csv')
        PROPOSAL_ID = lco_secrets['PROPOSAL_ID'][0]
        API_TOKEN = lco_secrets['API_TOKEN'][0]

        from lco import submit_photometric_observation
        from lco import submit_spectroscopic_observation

        obsid = submit_photometric_observation(name, ra, dec,
                                              PROPOSAL_ID, API_TOKEN,
                                              tstart=tstart, tend=tend,
                                              exposure_time = 300,
                                              doSubmission=True)
        message.append('View the observing request: https://observe.lco.global/requestgroups/{}/'.format(obsid))

        #submit_spectroscopic_observation(name, ra, dec,
        #                                 PROPOSAL_ID, API_TOKEN,
        #                                 tstart=tstart, tend=tend,
        #                                 exposure_time = 300,
        #                                 doSubmission=True)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
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
        #time.sleep(15)
