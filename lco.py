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

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm

from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta

import requests

def check_candidate_status(requestgroup):
 
    program = requestgroup['proposal']
    obsid = requestgroup["id"]
    target = requestgroup['requests'][0]["configurations"][0]["target"]["name"]
    instconfigs = requestgroup['requests'][0]["configurations"][0]["instrument_configs"]
    obs = []
    for instconfig in instconfigs:
        if 'filter' in instconfig['optical_elements']:
            obs.append(instconfig['optical_elements']['filter'])
        elif 'slit' in instconfig['optical_elements']:
            obs.append(instconfig['optical_elements']['slit'])
    completed = -1
    if requestgroup['state'] == "COMPLETED":
        completed = 2
    elif requestgroup['state'] == "PENDING":
        completed = 1
    elif requestgroup['state'] == "WINDOW_EXPIRED":
        completed = 0

    return target, completed, obs, program, obsid

def check_observations(API_TOKEN, lco_programs=None):

    targets = {}
    response_link = 'https://observe.lco.global/api/requestgroups/?'
    while response_link is not None:
  
        response = requests.get(response_link,
            headers={'Authorization': 'Token {}'.format(API_TOKEN)})
 
        # Make sure the API call was successful
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            print('Request failed: {}'.format(response.content))
            raise exc
    
        # Loop over the returned RequestGroups and print some information about them
        for requestgroup in response.json()['results']:
            target, completed, obs, program, obsid = check_candidate_status(requestgroup)
            if lco_programs is not None:
                if program not in lco_programs: continue

            if not target in targets:
                targets[target] = {}

            targets[target][obsid] = {}
            targets[target][obsid]["completed"] = completed
            targets[target][obsid]["observations"] = obs
            targets[target][obsid]["program"] = program
            targets[target][obsid]["obsid"] = obsid

        response_link = response.json()["next"]

    return targets

def submit_spectroscopic_observation(objname, ra, declination,
                                     PROPOSAL_ID, API_TOKEN,
                                     tstart, tend,
                                     exposure_time=300,
                                     doSubmission=False):

    # Constraints used for scheduling this observation
    constraints = {
        'max_airmass': 2,
        'min_lunar_distance': 30
    }

    # The target of the observation
    target = {
        'name': objname,
        'type': 'ICRS',
        'ra': ra,
        'dec': declination,
        'epoch': 2000
    }

    # The telescope class that should be used for this observation
    location = {
        'telescope_class': '2m0'
    }

    configurations = [
    {
        'type': 'LAMP_FLAT',
        'instrument_type': '2M0-FLOYDS-SCICAM',
        'constraints': constraints,
        'target': target,
        'acquisition_config': {},
        'guiding_config': {
            'mode': 'OFF',
            'optional': False},
        'instrument_configs': [
            {
                'exposure_time': 50,
                'exposure_count': 1,
                'rotator_mode': 'VFLOAT',
                'optical_elements': {
                    'slit': 'slit_1.6as'
                }
            }
        ]
    },
    {   
        'type': 'ARC',
        'instrument_type': '2M0-FLOYDS-SCICAM',
        'constraints': constraints,
        'target': target,
        'acquisition_config': {},
        'guiding_config': {
            'mode': 'OFF',
            'optional': False},
        'instrument_configs': [
            {
                'exposure_time': 60,
                'exposure_count': 1,
                'rotator_mode': 'VFLOAT',
                'optical_elements': {
                    'slit': 'slit_1.6as'
                }
            }
        ]
    },
    {
        'type': 'SPECTRUM',
        'instrument_type': '2M0-FLOYDS-SCICAM',
        'constraints': constraints,
        'target': target,
        'acquisition_config': {
            'mode': 'WCS'
        },
        'guiding_config': {
            'mode': 'ON',
            'optional': False
        },
        'instrument_configs': [
            {
                'exposure_time': exposure_time,
                'exposure_count': 1,
                'rotator_mode': 'VFLOAT',
                'optical_elements': {
                    'slit': 'slit_1.6as'
                }
            }
        ]
    },
    {
        'type': 'ARC',
        'instrument_type': '2M0-FLOYDS-SCICAM',
        'constraints': constraints,
        'target': target,
        'acquisition_config': {},
        'guiding_config': {
            'mode': 'OFF',
            'optional': False},
        'instrument_configs': [
            {
                'exposure_time': 60,
                'exposure_count': 1,
                'rotator_mode': 'VFLOAT',
                'optical_elements': {
                    'slit': 'slit_1.6as'
                }
            }
        ]
    },
    {
        'type': 'LAMP_FLAT',
        'instrument_type': '2M0-FLOYDS-SCICAM',
        'constraints': constraints,
        'target': target,
        'acquisition_config': {},
        'guiding_config': {
            'mode': 'OFF',
            'optional': False},
        'instrument_configs': [
            {
                'exposure_time': 50,
                'exposure_count': 1,
                'rotator_mode': 'VFLOAT',
                'optical_elements': {
                    'slit': 'slit_1.6as'
                }
            }
        ]
    }]

    windows = [{
        'start': tstart,
        'end': tend
    }]

    # The full RequestGroup, with additional meta-data
    requestgroup = {
        'name': objname,
        'proposal': PROPOSAL_ID,
        'ipp_value': 0.95,
        'operator': 'SINGLE',
        'observation_type': 'NORMAL',
        'requests': [{
            'configurations': configurations,
            'windows': windows,
            'location': location,
        }]
    }

    if doSubmission:
        # Submit the fully formed RequestGroup
        response = requests.post(
            'https://observe.lco.global/api/requestgroups/',
            headers={'Authorization': 'Token {}'.format(API_TOKEN)},
            json=requestgroup  # Make sure you use json!
        )

        # Make sure this api call was successful
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            print('Request failed: {}'.format(response.content))
            raise exc

        requestgroup_dict = response.json()  # The API will return the newly submitted RequestGroup as json

        # Print out the url on the portal where we can view the submitted request
        print('View the observing request: https://observe.lco.global/requestgroups/{}/'.format(requestgroup_dict['id']))

        print('Sleeping for 5 seconds')
        time.sleep(5)

def submit_photometric_observation(objname, ra, declination,
                                   PROPOSAL_ID, API_TOKEN,
                                   tstart, tend,
                                   filters=['gp','rp'],
                                   exposure_time=300,
                                   doSubmission=False):

    # Constraints used for scheduling this observation
    constraints = {
        'max_airmass': 2,
        'min_lunar_distance': 10
    }

    # The target of the observation
    target = {
        'name': objname,
        'type': 'ICRS',
        'ra': ra,
        'dec': declination,
        'epoch': 2000
    }
   
    # The configurations for this request. In this example we are taking 2 exposures with different filters.
    configurations = []
    for filt in filters:
        configurations.append({'type': 'EXPOSE',
                              'instrument_type': '1M0-SCICAM-SINISTRO',
                              'constraints': constraints,
                              'target': target,
                              'acquisition_config': {},
                              'guiding_config': {},
                              'instrument_configs': [
                                  { 
                                      'exposure_time': exposure_time,
                                      'exposure_count': 1,
                                      'optical_elements':
                                          {
                                              'filter': filt
                                          }
                                   }
                              ]
                          })
 
    windows = [{
        'start': tstart,
        'end': tend
    }]

    # The telescope class that should be used for this observation
    location = {
        'telescope_class': '1m0'
    }
    
    # The full RequestGroup, with additional meta-data
    requestgroup = {
            'name': '%s' % (objname),  # The title
            'proposal': PROPOSAL_ID,
            'ipp_value': 0.95,
            'operator': 'SINGLE',
            'observation_type': 'NORMAL',
            'requests': [{
                'configurations': configurations,
                'windows': windows,
                'location': location,
            }]
        }

    if doSubmission:
        # Submit the fully formed RequestGroup
        response = requests.post(
            'https://observe.lco.global/api/requestgroups/',
            headers={'Authorization': 'Token {}'.format(API_TOKEN)},
            json=requestgroup  # Make sure you use json!
        )
            
        # Make sure this api call was successful
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            print('Request failed: {}'.format(response.content))
            raise exc
        
        requestgroup_dict = response.json()  # The API will return the newly submitted RequestGroup as json
        
        # Print out the url on the portal where we can view the submitted request
        print('View the observing request: https://observe.lco.global/requestgroups/{}/'.format(requestgroup_dict['id']))

        print('Sleeping for 5 seconds')
        time.sleep(5)

    return requestgroup_dict['id']

def delete_observation(uid, API_TOKEN):
    
    response = requests.post(
        'https://observe.lco.global/api/requestgroups/{}/cancel/'.format(uid),
        headers={'Authorization': 'Token {}'.format(API_TOKEN)}
    )

    # Make sure this api call was successful
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print('Request failed: {}'.format(response.content))
        raise exc

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Submit LCO.')
    parser.add_argument('--name', dest='name', type=str,
                        required=True, help='Object Name')
    parser.add_argument('--ra', dest='ra', type=float,
                        required=True, help='Right Ascension')
    parser.add_argument('--declination', dest='declination', type=float,
                        required=True, help='Declination')
    parser.add_argument('--filt', dest='filt', type=str, required=True,
                        help='Filter set', default = 'gp,rp')
    parser.add_argument('--exposure_time', dest='exposure_time', type=float,
                        required=False, help='Exposure Time')
    parser.add_argument('--date-start', dest='date_start', type=str,
                        required=False,
                        help="Start date of the query, in ISO format. \
                        Example: '2017-08-17 12:41:04.4'", default=None)
    parser.add_argument('--date-end', dest='date_end', type=str,
                        required=False,
                        help="End date of the query, in ISO format. \
                        Example: '2017-08-18 12:00:00.0'", default=None)
    parser.add_argument('--API_TOKEN', dest='API_TOKEN', type=str,
                        required=True,
                        help="LCO API Token")
    parser.add_argument('--PROPOSAL_ID', dest='PROPOSAL_ID', type=str,
                        required=True,
                        help="Proposal ID")
    parser.add_argument("--doSubmission",  action="store_true", default=False)

    args = parser.parse_args()

    # Parse command line
    exposure_time = args.exposure_time
    # API token obtained from https://observe.lco.global/accounts/profile/
    API_TOKEN = args.API_TOKEN
    # Proposal IDs may be found here: https://observe.lco.global/proposals/
    PROPOSAL_ID = args.PROPOSAL_ID 
    filters = opts.filt.split(",")
    
    if not opts.tstart is None:
        tstart = Time(args.tstart)
        tend = Time(args.tend)
    else:
        tstart = Time.now()
        tend = Time.now() + TimeDelta(14*u.day)    
    
    tstart = str(tstart.isot).replace("T"," ")
    tend = str(tend.isot).replace("T"," ")

    submit_observation(name, ra, dec,
                       PROPOSAL_ID, API_TOKEN, 
                       tstart=tstart, tend=tend,
                       exposure_time = exposure_time,
                       doSubmission=args.doSubmission)
