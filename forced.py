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

from ForcePhotZTF.force_lc import get_cutout_data, download_images_diffpsf_refdiff
from ForcePhotZTF.force_mcmc import get_forced_phot_mcmc
from ForcePhotZTF.force_maxlike import get_forced_phot_maaxlike
from ForcePhotZTF.refine_lc import get_recerence_jds

def do_forcephot(targetdir, name, ra, dec, start_jd, end_jd,
                 ncpu=1, limit=100, do_mcmc=False,
                 program_ids = [1,2], verbose=False):
    """Perform forced photometry"""

    lightcurvesdir = "%s/%s" % (targetdir, "lightcurves")
    if not os.path.isdir(lightcurvesdir):
        os.makedirs(lightcurvesdir)

    targetdir = targetdir + "/"
    forcedfile = os.path.join(lightcurvesdir, f"{name}_force_phot_nob.csv")

    if not os.path.isfile(forcedfile):

        print(f">>> Forced photometry for target {name}")
        print(">>> Downloading images..")
        under_limit = download_images_diffpsf_refdiff(targetdir, ra, dec,
                                                      start_jd=start_jd,
                                                      end_jd=end_jd)
        if under_limit is True or under_limit is None:
            pass
        else:
            print(f"Exiting forced photometry for {name}")
            return
        print(">>> Images Downloaded")
        r_psf = 3
        r_bkg_in = 10
        r_bkg_out = 15
        print(">>> Getting cutouts...")
        get_cutout_data(name, targetdir, ra, dec, r_psf=r_psf,
                        r_bkg_in=r_bkg_in, r_bkg_out=r_bkg_out, verbose=verbose)
        print(">>> Got cutouts")
        get_recerence_jds(name, targetdir, only_partnership=False, retain_iband=True,
                          oldsuffix='_info.fits', newsuffix='_info_ref.fits', verbose=True)
        print(">>> Running forced phot")
        if do_mcmc is True:
            get_forced_phot_mcmc(name, targetdir, ncpu)
        else:
            get_forced_phot_maaxlike(name, targetdir, ra, dec)
        print(f">>> Done {name}")

        # Clean after yourself
        os.system(f'rm -f {targetdir}/images_*/*fits')

    figsdir = "%s/%s" % (targetdir, "figures")
    if not os.path.isdir(figsdir):
        os.makedirs(figsdir)

    data = pd.read_csv(forcedfile) 

    idxs = []
    for ii, programid in  enumerate(data['programid']):
        if programid in program_ids:
            idxs.append(ii)
    data = data.iloc[idxs]    

    # plotting code
    times = Time(data['jdobs'], format='jd')
    daysago = TimeDelta(times - Time.now(), scale='tt')

    filts = []
    for ii, filt in enumerate(data['filter']):
        filts.append(filt.replace("'","").replace("b",""))
    filts = np.array(filts)
    g = np.where(filts == 'g')[0]
    r = np.where(filts == 'r')[0]
    i = np.where(filts == 'i')[0]

    fig = plt.figure(figsize=(8,6))
    plt.errorbar(daysago[r].value, data["Fmcmc"].iloc[r], yerr=data['Fmcmc_unc'].iloc[r],
                 fmt='.', color='r', label='r-band')
    plt.errorbar(daysago[i].value, data['Fmcmc'].iloc[i], yerr=data['Fmcmc_unc'].iloc[i],
                 fmt='.', color='goldenrod', label='i-band')
    plt.errorbar(daysago[g].value, data['Fmcmc'].iloc[g], yerr=data['Fmcmc_unc'].iloc[g],
                 fmt='.', color='g', label='g-band')
    #plt.ylim(23.0, 18.0)
    #plt.xlim(200, 0)
    plt.xlabel('Days Ago')
    plt.ylabel('Flux')
    plotname = os.path.join(figsdir, '%s.png' % name)
    plt.savefig(plotname)
    plt.grid()
    plt.show()
    plt.close()

    return fig

def main_maxlike(name, t0, targetdir_base,
                 daydelta_before, daydelta_after):
    """Main function to launch maxlike forced phot"""

    t = t0[t0['name'] == name]

    # Increase the precision of the localization
    ra, dec = np.median(np.array(t['ra'])), np.median(np.array(t['dec']))

    # Start and end dates of the light curve
    start_jd = np.max(t['jd']) - daydelta_before
    end_jd = np.max(t['jd']) + daydelta_after

    targetdir = f"{targetdir_base}{name}/"
    # len_maxlike = len(glob.glob(f"{targetdir}/lightcurves/force_phot_{name}_maxlikelihood_lc.fits")) 
    try:
        do_forcephot(targetdir, name, ra, dec, start_jd, end_jd, ncpu=1, limit=2000)
        print(f"Completed forced photometry for {name}")
        gc.collect()
        status = 'success'
    except:
        status = 'problematic'
        print(f"Problematic: {name} failed forced phot")

    # Clean after yourself
    all_images = glob.glob(f"{targetdir}/images_*/*fits")
    subprocess.call(['rm', "-f"] + all_images)

    return (name, status)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Submit LCO.')
    parser.add_argument('--name', dest='name', type=str,
                        required=True, help='Object Name')
    parser.add_argument('--ra', dest='ra', type=float,
                        required=True, help='Right Ascension')
    parser.add_argument('--declination', dest='declination', type=float,
                        required=True, help='Declination')
    parser.add_argument('--date-start', dest='date_start', type=str,
                        required=False,
                        help="Start date of the query, in ISO format. \
                        Example: '2017-08-17 12:41:04.4'", default=None)
    parser.add_argument('--date-end', dest='date_end', type=str,
                        required=False,
                        help="End date of the query, in ISO format. \
                        Example: '2017-08-18 12:00:00.0'", default=None)
    parser.add_argument("-c", "--channel", type=str, default="partnership")

    parser.add_argument("-o", "--outdir", type=str, default="forcedphot")

    args = parser.parse_args()

    # Parse command line
    if not args.date_start is None:
        tstart = Time(args.date_start)
        tend = Time(args.date_end)
    else:
        tstart = Time.now() - TimeDelta(100*u.day)
        tend = Time.now()

    if args.channel == 'partnership':
        program_ids = [1,2]
    elif args.channel == 'caltech':
        program_ids = [1,2,3]
    else:
        print('Sorry, I do not know that channel...')
        exit(0)

    do_forcephot(args.outdir, args.name, args.ra, args.declination,
                 tstart.jd, tend.jd,
                 ncpu=1, limit=100, do_mcmc=True,
                 program_ids=program_ids, verbose=False)
