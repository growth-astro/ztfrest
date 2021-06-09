#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 15:05:27 2017

@author: kburdge
"""

import os, sys
import time
import optparse
import pandas as pd
import numpy as np
import glob

from io import StringIO

import requests

def submit_forced_photometry(name, ra, declination,
                             API_TOKEN, outputDir = './atlas'):

    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)

    outfile = os.path.join(outputDir,'%s.csv' % name)
    if os.path.isfile(outfile):
        return outfile

    BASEURL = "https://fallingstar-data.com/forcedphot"

    headers = {'Authorization': f'Token {API_TOKEN}', 'Accept': 'application/json'}
    task_url = None
    while not task_url:
        with requests.Session() as s:
            # alternative to token auth
            # s.auth = ('USERNAME', 'PASSWORD')
            resp = s.post(f"{BASEURL}/queue/", headers=headers, data={
                'ra': ra, 'dec': declination, 'mjd_min': 59200., 'send_email': False})
    
            if resp.status_code == 201:  # successfully queued
                task_url = resp.json()['url']
                print(f'The task URL is {task_url}')
            elif resp.status_code == 429:  # throttled
                message = resp.json()["detail"]
                print(f'{resp.status_code} {message}')
                t_sec = re.findall(r'available in (\d+) seconds', message)
                t_min = re.findall(r'available in (\d+) minutes', message)
                if t_sec:
                    waittime = int(t_sec[0])
                elif t_min:
                    waittime = int(t_min[0]) * 60
                else:
                    waittime = 10
                print(f'Waiting {waittime} seconds')
                time.sleep(waittime)
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.text)
                sys.exit()
    
    
    result_url = None
    taskstarted_printed = False
    while not result_url:
        with requests.Session() as s:
            resp = s.get(task_url, headers=headers)
    
            if resp.status_code == 200:  # HTTP OK
                if resp.json()['finishtimestamp']:
                    result_url = resp.json()['result_url']
                    print(f"Task is complete with results available at {result_url}")
                elif resp.json()['starttimestamp']:
                    if not taskstarted_printed:
                        print(f"Task is running (started at {resp.json()['starttimestamp']})")
                        taskstarted_printed = True
                    time.sleep(2)
                else:
                    print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
                    time.sleep(4)
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.text)
                sys.exit()
    
    with requests.Session() as s:
        textdata = s.get(result_url, headers=headers).text
    
        # if we'll be making a lot of requests, keep the web queue from being
        # cluttered (and reduce server storage usage) by sending a delete operation
        # s.delete(task_url, headers=headers).json()
    
    dfresult = pd.read_csv(StringIO(textdata.replace("###", "")), delim_whitespace=True)
    dfresult.to_csv(outfile)

    return outfile


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Submit LCO.')
    parser.add_argument('--name', dest='name', type=str,
                        required=True, help='Object Name')
    parser.add_argument('--ra', dest='ra', type=float,
                        required=True, help='Right Ascension')
    parser.add_argument('--declination', dest='declination', type=float,
                        required=True, help='Declination')
    parser.add_argument('--API_TOKEN', dest='API_TOKEN', type=str,
                        required=True,
                        help="ATLAS API Token")

    args = parser.parse_args()

    outfile = submit_forced_photometry(args.name, args.ra, args.declination,
                                       args.API_TOKEN)

