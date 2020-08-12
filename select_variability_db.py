"""
author: Igor Andreoni
email: andreoni@caltech.edu

Select significantly variable light curves
obtained with forced photometry.
"""

import os
import pdb
import glob

from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy as np
import pandas as pd
from scipy import optimize
from ztfquery import query
from astroquery.vizier import Vizier

from functions_db import connect_database

Vizier.ROW_LIMIT = 999999999  # Remove the limit of 50 rows



def stack_lc(tbl, days_stack=1., snt_det=3, snt_ul=5):
    """Given a dataframe with a maxlike light curve,
    stack the flux """

    if 'jdobs' in list(tbl.colnames):
        key_jd = 'jdobs'
    elif 'jd' in list(tbl.colnames):
        key_jd = 'jd'
    else:
        print("What is the column for the JD??")
        pdb.set_trace()
    t_out = Table([[],[],[],[],[],[],[],[],[]],
                  names=(key_jd, 'flux', 'flux_unc', 'zp', 'ezp',
                         'mag', 'mag_unc', 'limmag', 'filter'),
                  dtype=('double', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'S'))
    # Bin separately by filter
    filters = list(set(tbl['filter']))
    for f in filters:

        t = tbl[tbl['filter'] == f]

        bins = np.arange(int(np.max(t[key_jd]) - np.min(t[key_jd]))+2)
        dt0 = np.min(t[key_jd]) - int(np.min(t[key_jd]))
        if dt0 <= 0.4:
            start = int(np.min(t[key_jd])) - 0.6
        else:
            start = int(np.min(t[key_jd])) + 0.4
        bins = bins + start
        for b in bins:
            temp = t[(t[key_jd] > b) & (t[key_jd] < b+1)]
            if len(temp) == 0:
                continue
            new_jd = np.mean(np.array(temp[key_jd]))

            if len(set(temp['zp'])) == 1:
                zp = temp['zp'][0]
                flux = np.array(temp['Flux_maxlike'])
                flux_unc = np.array(temp['Flux_maxlike_unc'])
                flux[np.isnan(flux)] = 0
                # Use weights only if there are only detections
                if np.min(flux/flux_unc) >= snt_det:
                    weights = (flux/flux_unc)**2
                    new_flux = np.sum(np.array(temp['Flux_maxlike'])*weights)/np.sum(weights)
                else:
                    new_flux = np.mean(np.array(temp['Flux_maxlike']))
                new_flux_unc = np.sqrt(np.sum(np.array(temp['Flux_maxlike_unc'])**2))/len(temp)
            else:
                zp = temp['zp'][0]
                flux1 = np.array(temp['Flux_maxlike'])
                flux1_unc = np.array(temp['Flux_maxlike_unc'])
                zp1 = np.array(temp['zp'])
                flux = 10**((2.5*np.log10(flux1) - zp1 + zp) / 2.5)
                flux_unc = 10**((2.5*np.log10(flux1_unc) - zp1 + zp) / 2.5)
                flux[np.isnan(flux)] = 0
                # Use weights only if there are only detections
                if np.min(flux/flux_unc) >= snt_det:
                    weights = (flux/flux_unc)**2
                    new_flux = np.sum(flux*weights)/np.sum(weights)
                else:
                    new_flux = np.mean(flux)
                new_flux_unc = np.sqrt(np.sum(flux_unc**2))/len(temp)
            if new_flux/new_flux_unc > snt_det:
                mag_stack = -2.5*np.log10(new_flux) + zp
                mag_unc_stack = np.abs(-2.5*np.log10(new_flux-new_flux_unc) + 2.5*np.log10(new_flux))
                maglim_stack = 99.
            else:
                mag_stack = 99.
                mag_unc_stack = 99.
                maglim_stack = -2.5 * np.log10(snt_ul * new_flux_unc) + zp
            ezp = np.sum(temp['ezp']**2)/len(temp)
            t_out.add_row([new_jd, new_flux, new_flux_unc, zp, ezp, mag_stack,
                           mag_unc_stack, maglim_stack, f])

    return t_out


def query_metadata(ra, dec, username, password, zquery,
                   start_jd=None, end_jd=None,
                   out_csv=None):
    """Use ZTFquery to get more reliable upper limits"""

    if start_jd is None and end_jd is None:
        zquery.load_metadata(kind = 'sci', radec = [str(ra.deg), str(dec.deg)], size = 0.003,
                             auth=[username, password])
    else:
        if start_jd is not None and end_jd is None:
            sql_query='obsjd>'+repr(start_jd)
        elif start_jd is None and end_jd is not None:
            sql_query='obsjd<'+repr(end_jd)
        elif start_jd is not None and end_jd is not None:
            sql_query='obsjd<'+repr(end_jd)+'+AND+'+'obsjd>'+repr(start_jd)
        zquery.load_metadata(kind = 'sci', radec = [str(ra.deg), str(dec.deg)], size = 0.003,
                             sql_query=sql_query,
                             auth=[username, password])
    out = zquery.metatable
    final_out = out.sort_values(by=['obsjd'])
    if out_csv is not None:
        final_out.to_csv(out_csv)

    return final_out


def do_fit(errfunc, pinit, time, mag, magerr):
    out = optimize.leastsq(errfunc, pinit,
                           args=(time, mag, magerr), full_output=1)

    pfinal = out[0]
    covar = out[1]

    index = pfinal[1]
    amp = pfinal[0]

    try:
        indexErr = np.sqrt(covar[1][1])
    except TypeError:
        indexErr = 99.9
    try:
        ampErr = np.sqrt(covar[0][0]) * amp
    except TypeError:
        ampErr = 99.9

    return pfinal, covar, index, amp, indexErr, ampErr


def select_variability(tbl, hard_reject=[], update_database=False,
                       read_database=True,
                       use_forced_phot=False, stacked=False,
                       baseline=0.02, var_baseline={'g': 6, 'r': 8, 'i': 10},
                       max_duration_tot=30., max_days_g=1e5, snr=4,
                       index_rise=-0.0, index_decay=0.0,
                       path_secrets_db='db_access.csv',
                       save_plot=False, path_plot='./',
                       show_plot=False, use_metadata=False,
                       path_secrets_meta='../kowalski/secrets.csv',
                       save_csv=False, path_csv='./',
                       path_forced='./forced_phot/'):

    """
    Select candidates based on their duration and on their evolution rate.

    ---
    Parameters

    tbl astropy table
        table with the photometry of the candidates
        from the AVRO packets. It can be created using get_lc_kowalski.py
        and it usually includes:
        name,ra,dec,jd,magpsf,sigmapsf,filter,magzpsci,magzpsciunc,programid,
        field,rcid,pid,sgscore1,sgscore2,sgscore3,distpsnr1,distpsnr2,distpsnr3

    hard_reject list of str
        list of candidates that have to be ignored

    update_database bool
        if True, it updates the psql database with the results

    read_database bool
        if True, it reads light curves from the psql database

    use_forced_phot bool
        if True, forced ForcePhotZTF photometry will be used;
        if False, only alerts will be considered.

    stacked bool
        if True, the light curve will be stacked nightly in flux

    baseline float
        min time (days) between first and last detection for the fit
        to be performed (e.g.: baseline=2. means that if a transient
        is found only on a single night, it will not be fit to find
        its rise and decay rates)

    var_baseline dict
        if no evolution (1-sigma) is measured beyond these thresholds,
        the candidate is "rejected". Different filters have different 
        thresholds, the default being {'g': 6, 'r': 8, 'i': 10}

    max_duration_tot float
        max time (days) allowed between the first and last detections
        in any band

    max_days_g float
        max time (days) allowed between the first and last detections
        in g band

    snr float
        min signal-to-noise ratio for forced photometry data points;
        the maxlike method of ForcePhotZTF tends to underestimate the
        errors, so snr=4 is preferred over snr=3.
 
    index_rise float
        negative number, minimum rise rate (mag/day) for a candidate to
        be selected, if the rise index can be measured at all
        (e.g.: index_rise=-0.5 will allow you to select only
        those candidates rising faster than 0.5 mag/day).
 
    index_decay float
        positive number, minimum decay rate (mag/day) for a candidate to
        be selected, if the decay index can be measured at all
        (e.g.: index_decay=0.3 will allow you to select only
        those candidates decaying faster than 0.3 mag/day).

    path_secrets_db str
        path to the CSV secrets file to access the psql db.
        The file will need to have: 
        db,host,dbname,port,user,password

    path_forced str
        path to the directory where forced photometry light curves are stored

    save_plot bool
        save the plot to a png file

    path_plot str
         path to the folder where to save the plot

    show_plot
         display the plot while the code is running.
         The code will resume running when the plot is manually closed.

    use_metadata bool
        if True, use ztfquery to fetch upper limits based onZTF pipeline
        non-detection of the transient (for plotting only)

    path_secrets_meta str
        path to the CSV secrets file to access ztfquery.
        The file will need to have:
        ztfquery_user, ztfquery_pwd

    save_csv bool
        if True, the light curve is saved in a CSV file

    path_csv str
        path to the folder where the light curve will be saved

       
    ---
    Returns

    if update_database=True, it automatically updates the database with the
    following information:
    duration_tot, duration per band, rise or decay rate (index) per band

    selected list of str
        list of candidates that meet the input selection criteria

    rejected
        list of candidates that don't meet the input selection criteria

    cantsay
        list of candidates without enough information to tell
        if they meet the input selection criteria
    """

    # Useful definitions, specific for ZTF
    candidates = set(tbl["name"])
    filters = ['g', 'r', 'i']
    filters_id = {'1': 'g', '2': 'r', '3': 'i'}
    colors = {'g': 'g', 'r': 'r', 'i': 'y'}

    if update_database is True or read_database is True:
        # Connect to psql db
        con, cur = connect_database(update_database=update_database,
                                    path_secrets_db=path_secrets_db)

    if save_plot is True:
        if not os.path.isdir(path_plot):
            os.makedirs(path_plot) 
    if save_csv is True:
        if not os.path.isdir(path_csv):
            os.makedirs(path_csv)

    names_select = []
    names_reject = []
    empty_lc = []

    # Get forced phot for all the candidates
    candidates = list(candidates)[0:100]
    if read_database is True and use_forced_phot is True:
        str_names = "'" + "','".join(candidates) + "'"
        # table name
        if stacked is True:
            table_name = "lightcurve_stacked"
            column_names = "name, jd, flux, flux_unc, mag, mag_unc, \
limmag, filter, zp, ezp, programid, field, ccdid, qid"
        else:
            table_name = "lightcurve_forced"
            column_names = "name, jd, filter, programid, \
field, mag, mag_unc, limmag, zp, ezp, flux_maxlike, flux_maxlike_unc"
        # Read the light curve from the database
        t_pd = pd.read_sql_query(f"SELECT {column_names} \
                                 from {table_name} \
                                 where name IN ({str_names})", con)
        # If the table is empty, return
        if t_pd.empty:
            print("There is no forced photometry in the database \
for any of the given candidates!")

            return None, None, None

        t_forced = Table.from_pandas(t_pd)

    for name in candidates:
        # Is the candidate to be ignored?
        if name in hard_reject:
            continue
        # Check if the forced photometry light curve is available
        if use_forced_phot is True:
            with_forced_phot = True
            if read_database is True:
                t = t_forced[t_forced['name'] == name]
            else:
                # Read the light curve from a file
                files = glob.glob(f"{path_forced}/*{name}*maxlike*fits")
                if len(files) == 0:
                    print(f"No forced photometry available for {name}: skipping")
                    continue
                elif len(files) > 1:
                    print(f"WARNING: more than one light curve found for {name}")
                    print(f"Using {files[0]}")
                    filename = files[0]
                else:
                    filename = files[0]
                t = Table(fits.open(filename)[1].data)
            empty = False
            if len(t) == 0:
                empty = True
                empty_lc.append(name)
                print(f"Empty forced photometry light curve for {name}: skipping")
                continue
                # Keep going with only the alerts?
                #t = tbl[tbl['name'] == name]
                #t_ul = t[:0].copy()
            else:
                t_ul = t[t["mag"] > 50]
                t = t[t["mag"] < 50]
                # Fix the column names
                t.rename_column('mag', 'magpsf')
                t.rename_column('mag_unc', 'sigmapsf')
                if read_database is False:
                    t.rename_column('jdobs', 'jd')
                    t_ul.rename_column('jdobs', 'jd')

                # Add missing epochs packets
                for l in tbl[tbl['name'] == name]:
                    if len(t_ul) > 0:
                        min_delta_ul = np.min(np.abs(t_ul['jd'] - l['jd']))
                    else:
                        min_delta_ul = np.inf
                    if len(t) > 0:
                        min_delta_det = np.min(np.abs(t['jd'] - l['jd']))
                    else:
                        min_delta_det = np.inf
                    # Min time from forced photometry for using alerts
                    if stacked is True:
                        # 12 hours for stacked images
                        min_delta = 12./24.
                    else:
                        # 15 seconds for non-stacked light curves
                        min_delta = 15./24./60./60.
                    if np.min([min_delta_ul, min_delta_det]) > min_delta:
                        if stacked is False and read_database is False:
                            new_row = [l['jd'], l['filter'], np.nan,
                                       np.nan, np.nan, np.nan, 1, l['field'],
                                       l['rcid'], np.nan, np.nan, np.nan,
                                       np.nan, np.nan, np.nan, np.nan, np.nan,
                                       np.nan, np.nan, np.nan, '', '', np.nan,
                                       np.nan, np.nan, np.nan, np.nan,
                                       l['magpsf'], l['sigmapsf'], np.nan]
                        elif stacked is False and read_database is True:
                            new_row = [name, l['jd'], l['filter'], 1,
                                       l['field'], l['magpsf'], l['sigmapsf'],
                                       99., 0., 0., 0., 0.]
                        elif stacked is True and read_database is True:
                            new_row = [name, l['jd'], np.nan, np.nan,
                                       l['magpsf'], l['sigmapsf'], 99.,
                                       l['filter'], 0., 0., l['programid'], l['field'], 0., 0.]
                        else:
                            new_row = [l['jd'], np.nan, np.nan, np.nan, np.nan,
                                       l['magpsf'], l['sigmapsf'], np.nan,
                                       l['filter'], 1]
                        t.add_row(new_row)
        else:
            with_forced_phot = False
            empty = False
            t = tbl[tbl['name'] == name]
            t_ul = t[:0].copy()

        # Reject those with only upper limits
        if len(t) == 0 and empty is False:
            names_reject.append(name)
            continue

        # Determine the light curve starting time
        t0 = min(t["jd"])

        # Reject if the overall duration is longer than ??  days
        # or if there is only 1 detection
        try:
            if update_database is True:
                # FIXME do we want to report only the max duration_tot between alerts, forcephot, and stack?
                cur.execute(f"UPDATE candidate SET \
                            duration_tot = {np.max(t['jd']) - np.min(t['jd'])}\
                            where name = '{name}'")
            if np.max(t['jd']) - np.min(t['jd']) > max_duration_tot or np.max(t['jd']) - np.min(t['jd']) == 0:
                names_reject.append(name)
                continue
        except ValueError:
            print("Failed calculating max(t['jd']) - min(t['jd']) > 10.")
            pdb.set_trace()
        try:
            plt.close()
        except:
            pass
        #plt.clf()
        plt.figure(figsize=(8,6))
        plt.subplot(1, 1, 1)

        plotted = False
        #print(f"-------- {name}")

        for f in filters:
            tf = t[t['filter'] == f]
            if len(tf) == 0:
                continue
            if use_metadata is False:
                tf_ul = t_ul[t_ul['filter'] == f]
                if len(tf_ul) > 0:
                    tf_ul["jd"] = tf_ul["jd"] - t0
                    plt.plot(np.array(tf_ul["jd"]), np.array(tf_ul["limmag"]),
                             colors[f]+'v', markeredgecolor=colors[f],
                             markerfacecolor='w')
                    plt.plot([],[], 'kv', label='UL')
            # Correct the start time
            tf["jd"] = tf["jd"] - t0

            #brightest, faintest detections
            bright = np.min(tf["magpsf"])
            try:
                brighterr = tf["sigmapsf"][tf["magpsf"] == bright][0]
            except:
                print("problems with",
                      "brighterr = tf['sigmapsf'][tf['magpsf'] == bright][0]")
                print(tf)
                pdb.set_trace()
                continue
            bright_jd = tf["jd"][tf["magpsf"] == bright][0]
            faint = np.max(tf["magpsf"])
            fainterr = tf["sigmapsf"][tf["magpsf"] == faint][0]
            faint_jd = tf["jd"][tf["magpsf"] == faint][0]

            # First and last detections
            first = np.min(tf["jd"])
            last = np.max(tf["jd"])

            # Add the information regarding the duration in the db
            # duration_g is the max number of days between the first
            # detection and the last detection in g band.
            if update_database is True:
                cur.execute(f"UPDATE candidate SET \
                            duration_{f} = {last-first} \
                            where name = '{name}'")

            time = np.array(tf["jd"])
            mag = np.array(tf["magpsf"])
            magerr = np.array(tf["sigmapsf"])
            plt.errorbar(np.array(tf["jd"][tf['programid']!=1]),
                         np.array(tf["magpsf"][tf['programid']!=1]),
                         fmt=colors[f]+'s',
                         yerr=np.array(tf["sigmapsf"][tf['programid']!=1]),
                         markeredgecolor='k', markersize=8)
            plt.errorbar(np.array(tf["jd"][tf['programid']==1]),
                         np.array(tf["magpsf"][tf['programid']==1]),
                         fmt=colors[f]+'o',
                         yerr=np.array(tf["sigmapsf"][tf['programid']==1]),
                         markeredgecolor='k', markersize=8)

            plt.plot([],[], 'ks', label='programid=2,3')
            plt.plot([],[], 'ko', label='programid=1')

            # max_days_g is the max number of days between the first
            # detection in any band and the last g-band detection
            if update_database is True:
                cur.execute(f"UPDATE candidate SET \
                            max_days_{f} = {np.max(tf['jd'])} \
                            where name = '{name}'")

            # SELECT: not enough baseline - no action taken
            if np.abs(last-first) < baseline:
                continue

            # SELECT: no variability between the first and last detection - rejection
            if bright+brighterr > faint-fainterr and np.abs(bright_jd - faint_jd) >= var_baseline[f]:
                names_reject.append(name)
                continue

            # SELECT: if a g-band detection is present xx days
	    # after the first detection, reject
            if f == 'g' and np.max(tf['jd']) > max_days_g:
                names_reject.append(name)
                continue

            onlyrise = False
            onlyfade = False

            if bright_jd < first + baseline:
                onlyfade = True
                riseorfade = 'fade'
            elif bright_jd > last - baseline:
                onlyrise = True
                riseorfade = 'rise'

            # Fit
            fitfunc = lambda p, x: p[0] + p[1] * x
            errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
            pinit = [1.0, -1.0]

            if onlyrise or onlyfade:
                pfinal, covar, index, amp, indexErr, ampErr = do_fit(errfunc, pinit, time, mag, magerr)
                plt.plot(time, fitfunc(pfinal, time), color=colors[f],
                         label=f"{f}, index= {'{:.2f}'.format(index)}+-{'{:.2f}'.format(indexErr)}")

                # Add info to the database
                if update_database is True and with_forced_phot is False:
                    cur.execute(f"UPDATE candidate SET \
                                index_{riseorfade}_{f} = {index} \
                                where name = '{name}'")
                elif update_database is True and with_forced_phot is True:
                    if stacked is True:
                        column = f"index_{riseorfade}_stack_{f}"
                    else:
                        column = f"index_{riseorfade}_forced_{f}"
                    cur.execute(f"UPDATE candidate SET \
                                {column} = {index} \
                                where name = '{name}'")

                # SELECT: slow evolution over a time longer than the baseline
                if ((index > 0 and index <= index_decay) and
                    (last-first) > var_baseline[f]):
                    names_reject.append(name)
                else:
                    plotted = True
            else:
                indexrise = np.where(time <= bright_jd)
                indexfade = np.where(time >= bright_jd)
                for i, riseorfade in zip([indexrise, indexfade], ['rise', 'fade']):
                    time_new = time[i[0]]
                    mag_new = mag[i[0]]
                    magerr_new = magerr[i[0]]

                    faint_new = np.max(mag_new)
                    fainterr_new = magerr_new[np.where(mag_new == faint_new)[0]]
                    faint_jd_new = time_new[np.where(mag_new == faint_new)[0]]
                    # Format check
                    if type(magerr_new) is list or type(magerr_new) is np.ndarray:
                        magerr_new = magerr_new[0]
                    if type(fainterr_new) is list or type(fainterr_new) is np.ndarray:
                        fainterr_new = fainterr_new[0]
                    if type(faint_jd_new) is list or type(faint_jd_new) is np.ndarray:
                        faint_jd_new = faint_jd_new[0]

                    # SELECT: no evolution
                    try:
                        if (bright+brighterr > faint_new - fainterr_new and
                            np.abs(faint_jd_new - bright_jd) >= var_baseline[f]):
                            names_reject.append(name)
                            plt.errorbar(time, mag, yerr=magerr,
                                         fmt=colors[f]+'.',
                                         markeredgecolor='k',
                                         markersize=8)
                            continue
                    except:
                        print(bright, brighterr, faint_new, fainterr_new)
                        pdb.set_trace()
                    pfinal, covar, index, amp, indexErr, ampErr = do_fit(errfunc, pinit, time_new, mag_new, magerr_new)
                    plt.plot(time_new, fitfunc(pfinal, time_new),
                             color=colors[f],
                             label=f"{f}, index= {'{:.2f}'.format(index)}+-{'{:.2f}'.format(indexErr)}")

                    if update_database is True and with_forced_phot is False:
                        cur.execute(f"UPDATE candidate SET \
                                    index_{riseorfade}_{f} = {index} \
                                    where name = '{name}'")
                    elif update_database is True and with_forced_phot is True:
                        cur.execute(f"UPDATE candidate SET \
                                    index_{riseorfade}_forced_{f} = {index} \
                                    where name = '{name}'")

                    # SELECT: slow evolution
                    if ((index > 0 and index <= index_decay) and
                        (last-first) > var_baseline[f]):
                        names_reject.append(name)
                    else:
                        plotted = True

        if name in names_reject:
            plt.close()
            continue

        if plotted is True:
            # The candidate was selected!
            names_select.append(name)

        if plotted is True and (show_plot is True or save_plot is True):
            # The following is for the file naming
            if use_forced_phot is True:
                forcedbool = 1
            else:
                forcedbool = 0
            if stacked is True:
                stackbool = 1
            else:
                stackbool = 0

            if use_metadata:
                # Fetch metadata on all the available ZTF images
                #start_jd = Time('2018-03-01 00:00:00', format='iso').jd
                #end_jd =  Time('2030-02-01 04:50:59.998', format='iso').jd
                start_jd, end_jd = None, None

                # Read the secrets
                secrets = ascii.read(path_secrets_meta, format='csv')
                username = secrets['ztfquery_user'][0]
                password = secrets['ztfquery_pwd'][0]

                zquery = query.ZTFQuery()
                coords = SkyCoord(ra=np.mean(tbl[tbl['name'] == name]['ra']*u.deg),
                                  dec=np.mean(tbl[tbl['name'] == name]['dec']*u.deg))
                metadata = query_metadata(coords.ra, coords.dec, username, password,
                                          zquery, start_jd=start_jd, end_jd=end_jd,
                                          out_csv=None)
                t_ul = Table([[],[],[],[],[],[],[],[]],
                             names=('jd', 'magpsf', 'sigmapsf', 'filter',
                                   'snr', 'limmag', 'seeing', 'programid'),
                             dtype=('double','f','f','S','f','f','f','int'))
                for j, ml, fid, s, pid in zip(metadata['obsjd'],
                                              metadata['maglimit'],
                                              metadata['fid'],
                                              metadata['seeing'],
                                              metadata['pid']):
                    if not (j in t['jd']):
                        new_row = [j, 99.9, 99.9, filters_id[str(fid)],
                                   np.nan, ml, s, 0]
                        t_ul.add_row(new_row)
                for f in filters:
                    tf_ul = t_ul[t_ul['filter'] == f]
                    if len(tf_ul) > 0:
                        tf_ul["jd"] = tf_ul["jd"] - t0
                        plt.plot(np.array(tf_ul["jd"]), np.array(tf_ul["limmag"]),
                                 colors[f]+'v', markeredgecolor=colors[f],
                                 markerfacecolor='w')
                        plt.plot([],[], 'kv', label='UL')
 
            plt.title(f"{name}")
            plt.xlabel('Time [days]', fontsize=18)
            plt.ylabel('mag', fontsize=18)

            plt.tick_params(axis='both',    # changes apply to the x-axis
                            which='both',   # both major and minor ticks are affected
                            labelsize=16)

            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), fontsize=16)

            plt.gca().invert_yaxis()

            if save_plot is True:
                # Save the image
                plt.savefig(f"{path_plot}/lc_{name}_forced{forcedbool}_stacked{stackbool}.png")
            if save_csv is True:
                # Save the light curve
                t_union = vstack([t, t_ul])
                ascii.write(t_union,
                            f"{path_csv}/lc_{name}_forced{forcedbool}_stacked{stackbool}.csv",
                            format='csv', overwrite=True)
            if show_plot is True:
                plt.show()
        else:
            plt.close()

    if update_database is True:
        # Update the hard_reject flag
        names_reject_str = "','".join(list(names_reject))
        cur.execute(f"UPDATE candidate SET \
                    hard_reject = 1 \
                    where name IN ('{names_reject_str}')")
        # Commit the changes
        con.commit()

        # Close the connection
        cur.close()
        con.close()

    print(f"{len(set(empty_lc))} empty light curves")
    print(f"Select {set(names_select)}")
    print(f"Reject {set(names_reject)}")
    cantsay = list(n for n in candidates if
                   not (n in names_select) and not (n in names_reject))
    print(f"Cannot say {set(cantsay)}")
    print(f"{len(set(names_reject))}/{len(candidates)} objects rejected")
    print(f"{len(set(names_select))}/{len(candidates)} objects selected")
    print(f"{len(set(cantsay))}/{len(candidates)} objects cannot say")

    return names_select, names_reject, cantsay


if __name__ == "__main__":
    filename1 = '../paper_kn_ZTF/lc_results_2y_min0015_max12_ndethist3.csv'
    tbl1 = ascii.read(filename1, format='csv')
    tbl = tbl1

    filename_reject = '../paper_kn_ZTF/hard_rejects_candidates2.csv'
    hard_reject = ascii.read(filename_reject, format='csv')['name']

    use_metadata = False  # Get upper limits with ZTFquery

    candidates = set(tbl["name"])
    filters = ['g', 'r', 'i']
    filters_id = {'1': 'g', '2': 'r', '3': 'i'}
    colors = {'g': 'g', 'r': 'r', 'i': 'y'}

    path_forced = '../paper_kn_ZTF/lc_all_ebv03_forced/'

    # Min time between for the fit to be performed
    baseline = 0.125
    # Min time over which you want to see any significant evolution
    var_baseline = {'g': 6, 'r': 8, 'i': 10}
    # Max number of days from the first detection for a g-band point to be present
    max_days_g = 15

    snr = 4

    index_rise = -2.0
    index_decay = 1.0

    # Do you want to update the db?
    update_database = False

    # Use forced photometry, if available:
    use_forced_phot = True

    # Do you want to stack the photometry?
    stacked = True

    selected, rejected, cantsay = select_variability(tbl,
                       hard_reject=[], update_database=False,
                       read_database=True,
                       use_forced_phot=True, stacked=False,
                       baseline=1.0, var_baseline={'g': 6, 'r': 8, 'i': 10},
                       max_duration_tot=15., max_days_g=10, snr=4,
                       index_rise=index_rise, index_decay=index_decay,
                       path_secrets_db='/Users/igor/data/paper_kn_ZTF/database/db_access.csv',
                       save_plot=True, path_plot='./plots/',
                       show_plot=False, use_metadata=False,
                       path_secrets_meta='../kowalski/secrets.csv',
                       save_csv=True, path_csv='./lc_csv',
                       path_forced=path_forced)
