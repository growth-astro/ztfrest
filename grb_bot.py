
# Import the relevant packages
import psycopg2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table, Column
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
from datetime import datetime

import healpy as hp

from ligo.skymap import plot, postprocess
import requests
from bs4 import BeautifulSoup

from penquins import Kowalski

from functions_misc import make_triplet, plot_triplet, get_cutouts
from functions_misc import get_dust_info, plot_lc
from functions_misc import get_xmatch_clu_glade, get_bgal_ebv

slack_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt")
with open(slack_token, "r") as f:
    access_token = f.read()
web_client = WebClient(token=access_token)

def GRB_unnoculted(path_GBM, path_out):
    """ Create a new skymap with the earth-oculted portion of the map masked
    
    new map in path.replace('.fit','_unocculted.fit')
    
    haversine and find_greedy_credible_levels from GBM_tools 
    https://fermi.gsfc.nasa.gov/ssc/data/analysis/rmfit/gbm_data_tools/gdt-docs/index.html
    """
    
    def haversine(lon1, lat1, lon2, lat2, deg=True):
        """Calculates the angular separation between two points using the
        haversine equation. If degrees are passed, degrees are returned. else
        the input/output is assumed to be radians.
        lon -> azimuth
        lat -> zenith

        Args:
            lon1 (float): lon/az of first point
            lat1 (float): lat/zen of first point
            lon2 (float): lon/az of second point
            lat2 (float): lat/zen of second point
            deg (bool, optional): True if input/output in degrees.

        Returns:
            float: Angular separation between points
        """
        if deg:
            lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
        d_lat = 0.5 * (lat2 - lat1)
        d_lon = 0.5 * (lon2 - lon1)

        a = np.sin(d_lat) ** 2 + (np.sin(d_lon) ** 2 * np.cos(lat1) * np.cos(lat2))
        alpha = 2. * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

        if deg:
            alpha = np.rad2deg(alpha)

        return alpha

    def find_greedy_credible_levels(p):
        """Calculate the credible values of a probability array using a greedy
        algorithm.

        Args:
            p (np.array): The probability array

        Returns:    
             np.array: The credible values
        """
        p = np.asarray(p)
        pflat = p.ravel()
        i = np.argsort(pflat)[::-1]
        cs = np.cumsum(pflat[i])
        cls = np.empty_like(pflat)
        cls[i] = cs
        return cls.reshape(p.shape)

    # getting header information
    hdu = fits.open(path_GBM)
    geo_ra = hdu['HEALPIX'].header['GEO_RA']
    geo_dec = hdu['HEALPIX'].header['GEO_DEC']
    geo_location = [geo_ra,geo_dec]

    try:
        geo_radius = hdu['HEALPIX'].header['GEO_RAD']
    except:
        geo_radius = 67.5

    prob, sig= hp.read_map(path_GBM, field=(0, 1), memmap=False,
                                    verbose=False)
    # prob_mask = prob > 0.0
    npix = len(prob)
    nside = hp.npix2nside(npix)

    theta, phi = hp.pix2ang(nside, np.arange(npix))
    ra = np.rad2deg(phi)
    dec = np.rad2deg(np.pi / 2.0 - theta)
    ang = haversine(*geo_location, ra, dec)

    # masking earth
    geo_mask = (ang <= geo_radius)
    new_prob = np.copy(prob)
    new_prob[geo_mask] = 0.0
    # renormalize
    new_prob /= np.sum(new_prob)
    # have to redo the significance
    new_sig = 1.0 - find_greedy_credible_levels(new_prob)

    unocculted_path = path_GBM.replace('.fit','_unocculted.fit')
    # get arrays in proper order, and write the healpix data to disk
    prob_arr = hp.reorder(new_prob, r2n=True)
    sig_arr = hp.reorder(new_sig, r2n=True)
    columns = ['PROBABILITY', 'SIGNIFICANCE']
    hp.write_map(unocculted_path, (prob_arr, sig_arr), nest=True, coord='C',
                 overwrite=True, \
                 column_names=columns,
                 extra_header=hdu[1].header.cards)
    hdulist = fits.open(unocculted_path)
    hdulist[0].header.extend(hdu[0].header.cards)
    hdulist[1].name = 'HEALPIX'
    hdulist[1].header['TTYPE1'] = (
    'PROBABILITY', 'Differential probability per pixel')
    hdulist[1].header['TTYPE2'] = (
    'SIGNIFICANCE', 'Integrated probability')
    hdulist.writeto(path_out, overwrite=True, checksum=True)

def in_out(map_path,ra_obj,dec_obj, top_fraction = 0.95 ):
    
    """
    Calculate if a list of objects are inside 'top_fraction' localisation region of a given skymap
    """
    # Read skymap, calculate top pixels
    #     top_fraction = 0.95 # limit skymap top 90% region
    skymap = hp.read_map(map_path)
    npix = len(skymap)
    nside = hp.npix2nside(npix)

    # Convert to astropy Table, easier for manipulation
    indices = np.arange(len(skymap))
    tm = Table(data=(indices, skymap), names=('id', 'prob'))
    tm.sort('prob')
    cs = np.cumsum(tm['prob'])
    cs.name='cumsum'
    tm.add_column(cs)

    top_pix = (tm['cumsum'] > 1 - top_fraction)
    tp = Column(data=top_pix, name="top")
    tm.add_column(tp)

    # Cast as a set for easier comparison below
    top_subset = set(tm['id'][tm['top']])

    inside = False
    pix = hp.ang2pix(nside,ra_obj,dec_obj, lonlat=True)
    if pix in top_subset:
        inside = True 
    return inside

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

def pretty_plot(map_path,ra,dec,cand_name='name'):

    # loading relevant info
    prob = hp.read_map(map_path)
    ipix = np.argmax(prob)
    npix = len(prob)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, ipix)

    # getting the center
    ra_max = np.rad2deg(phi)
    dec_max = np.rad2deg(0.5 * np.pi - theta)-10
    center = str(int(SkyCoord(ra_max,dec_max,unit='deg').ra.hms.h))+'h '+ str(int(SkyCoord(ra_max,dec_max,unit='deg').dec.dms.d))+'d'
    sky_area_per_pix = 4 * 180**2 / np.pi /npix
    area_95 = np.sum(np.cumsum(-np.sort(-prob))<=0.95) * sky_area_per_pix
    area_50 = np.sum(np.cumsum(-np.sort(-prob))<=0.5) * sky_area_per_pix
    print('AREA 50% ',area_50)
    print('AREA 95% ',area_95)

    #plotting
    fig = plt.figure()
    ax = plt.axes(projection='astro globe', center=center)
    ax.imshow_hpx(prob, cmap='cylon')
    cl = postprocess.find_greedy_credible_levels(prob) * 100
    ax.contour_hpx(cl, levels=[50, 95], colors=['black', 'black'], linewidths=[0.5, 0.5])
    ax.grid()
    ax.scatter(ra,dec,marker='*',color='blue',transform=ax.get_transform('world'),zorder=10)
    ax.set_title(fits.open(map_path)[0].header['OBJECT']+' / '+str(fits.open(map_path)[0].header['TRIGTIME'])+' \n'+cand_name)

    return fig

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
   
        doTrigger, trigger_action = False, 'check'
        for mess in payload["messages"]:
            print(mess)
            message_ts = float(mess["ts"])
            if np.abs(message_ts - thread_ts) > delay_thresh:
                continue
            txt = mess['text']
            txtsplit = list(filter(None,txt.split(" ")))
            if len(txtsplit) == 0: continue
            if txtsplit[0] == "grb":
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
        name, trigger_action = 'ZTF20aakqxsq', 'check'
        #name, trigger_action = 'ZTF20aamsouh', 'check'

    message = []
    message.append("Hi <@{0}>! You are interested in ztfrest GRBs, right? Let me get right on that for you.".format(user))
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

    jd = np.min(alerts['jd'])

    # get GRB list
    htmlpage = "https://gcn.gsfc.nasa.gov/fermi_grbs.html"
    page = requests.get(htmlpage)
    soup = BeautifulSoup(page.content, 'html.parser')
    rows = soup.find_all('td')

    trignums = []
    for ii, row in enumerate(rows):
        line = row.get_text()
        if "POS_MAP_URL" in line:
            trignum = int(rows[ii-7].get_text())
            trigdate = rows[ii-6].get_text()
            trigtime = rows[ii-5].get_text()
            skymap = list(filter(None,line.split(" ")))[-1]

            if trignum in trignums: continue

            year = int("20" + trigdate[:2])
            month = int(trigdate[3:5])
            day = int(trigdate[6:8])

            hour = int(trigtime[:2])
            minute = int(trigtime[3:5])
            second = int(trigtime[6:8])

            t = Time(datetime(year, month, day, hour, minute, second), scale='utc')
            dt = jd - t.jd
            if dt < 2 and dt > 0:

                fermipage = "https://gcn.gsfc.nasa.gov/other/%d.fermi" % trignum
                page = requests.get(fermipage)
                soup2 = BeautifulSoup(page.content, 'html.parser')
                rows_fermi = str(soup2).split("\n")
                for row_fermi in rows_fermi:
                    if 'LOC_URL' in row_fermi:
                        skymap = list(filter(None,row_fermi.split(" ")))[-1]
                        skymap = skymap.replace('http://', 'https://')
                        skymap = skymap.replace('_locplot_', '_healpix_')
                        skymap = skymap.replace('.png', '.fit')
                outdir = "skymaps/%d" % trignum
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                filename = os.path.join(outdir,skymap.split("/")[-1])
                if not os.path.isfile(filename):
                    wget_command = "wget %s" % skymap
                    os.system(wget_command)
                    mv_command = "mv %s %s" % (skymap.split("/")[-1], filename)
                    os.system(mv_command)

                unocculted_path = filename.replace('.fit','_unocculted.fit')
                GRB_unnoculted(filename, unocculted_path)
                inout = in_out(unocculted_path, ra, dec, top_fraction = 0.95 )
                fig = pretty_plot(unocculted_path, ra, dec, cand_name=name)
                upload_fig(fig, user, "skymap.png", channel_id)

                message = []
                message.append("Fermi Trigger number: %d" % (trignum))
                message.append("Delay to first detection: %.2f" % (dt))
                message.append("Within 95 percentile: %s" % (inout))

                web_client.chat_postMessage(
                    channel=channel_id,
                    text="\n".join(message)
                )

            trignums.append(trignum)

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
        print('Looking for some GRBing to do!')
        run_on_event(channel)
        #except:
        #    pass
        time.sleep(5)
