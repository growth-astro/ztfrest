
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

from forced import do_forcephot

slack_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt")
with open(slack_token, "r") as f:
    access_token = f.read()
web_client = WebClient(token=access_token)

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

def run_on_event(channel_id, bypass=False,
                 program_ids=[1,2], outdir="forcedphot"):

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
   
        doTrigger, trigger_action, name, ra, dec = False, 'forced', 'ZTF21abinaiu', 310.480103, 11.70041
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
                if len(txtsplit) == 4:
                    name = txtsplit[1]
                    ra = float(txtsplit[2])
                    dec = float(txtsplit[3])
            user = mess['user']
        if not doTrigger:
            return
    else:
        user, message_ts = 'test', thread_ts
        # options: gp,rp,ip,zs,Y
        name, ra, dec = 'ZTF21abinaiu', 310.480103, 11.70041

    message = []
    message.append("Hi <@{0}>! You are interested in forced photometry, right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking into %s for you" % name)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    tstart = Time.now() - TimeDelta(100*u.day)
    tend = Time.now()    

    fig = do_forcephot(outdir, name, ra, dec,
                       tstart.jd, tend.jd,
                       ncpu=1, limit=100, do_mcmc=True,
                       program_ids=program_ids, verbose=False)
    upload_fig(fig, user, "forced.png", channel_id)

    message = []
    message.append('Forced photometry complete')

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )    

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", type=str, default="partnership")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-o", "--outdir", type=str, default="forcedphot")

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

    if cfg.debug:
        run_on_event(channel, bypass=True, program_ids=program_ids, outdir=cfg.outdir)
        exit(0)

    while True:
        #try:
        print('Looking for forced photometry to do!')
        run_on_event(channel, program_ids=program_ids)
        #except:
        #    pass
        time.sleep(5)
