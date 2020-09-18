
# Import the relevant packages
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

slack_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt")
with open(slack_token, "r") as f:
    access_token = f.read()
web_client = WebClient(token=access_token)

def run_on_event(channel_id):

    #response = web_client.chat_postMessage(
    #    channel=channel_id,
    #    text='testing')

    converations = web_client.conversations_list(
        channel=channel_id
    )
    channel_slack_id = channel_id
    for converation in converations:
        for chan in converation.data["channels"]:
            if chan["name"] == channel_id.replace("#",""):
                channel_slack_id = chan["id"]

    message = []
    message.append("Hi everyone! Processing for the day is done. Now to do some scanning, right?")

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--timestamp", type=str, default=str(Time.now()))
    parser.add_argument("-c", "--channel", type=str, default="partnership")
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

    run_on_event(channel)
