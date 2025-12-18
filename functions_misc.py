"""
Suite of functions for plotting cutouts and light curves
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from astropy.stats import sigma_clipped_stats, mad_std
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord 
import astropy.units as u


def make_triplet(alert, normalize: bool = False):
    """
        Feed in alert packet
    """
    from bson.json_util import loads, dumps
    import gzip
    import io
    from astropy.io import fits
    from matplotlib.colors import LogNorm

    cutout_dict = dict()

    for cutout in ('science', 'template', 'difference'):
        cutout_data = loads(dumps([alert[f'cutout{cutout.capitalize()}\
']['stampData']]))[0]

        # unzip
        with gzip.open(io.BytesIO(cutout_data), 'rb') as f:
            with fits.open(io.BytesIO(f.read())) as hdu:
                data = hdu[0].data
                # replace nans with zeros
                cutout_dict[cutout] = np.nan_to_num(data)
                # normalize
                if normalize:
                    cutout_dict[cutout] /= np.linalg.norm(cutout_dict[cutout])

        # pad to 63x63 if smaller
        shape = cutout_dict[cutout].shape
        if shape != (63, 63):
            cutout_dict[cutout] = np.pad(cutout_dict[cutout],
                                         [(0, 63 - shape[0]),
                                          (0, 63 - shape[1])],
                                         mode='constant', constant_values=1e-9)

    triplet = np.zeros((63, 63, 3))
    triplet[:, :, 0] = cutout_dict['science']
    triplet[:, :, 1] = cutout_dict['template']
    triplet[:, :, 2] = cutout_dict['difference']

    return triplet


def plot_triplet(tr, show_fig=True):
    fig = plt.figure(figsize=(8, 2), dpi=120)
    ax1 = fig.add_subplot(131)
    ax1.axis('off')
    mean, median, std = sigma_clipped_stats(tr[:, :, 0])
    ax1.imshow(tr[:, :, 0], vmin = median - 2*std, vmax = median + 3*std)
    #ax1.imshow(tr[:, :, 0], origin='upper', cmap=plt.cm.bone, norm=LogNorm())
    ax1.title.set_text('Science')
    ax2 = fig.add_subplot(132)
    ax2.axis('off')
    mean, median, std = sigma_clipped_stats(tr[:, :, 1])
    ax2.imshow(tr[:, :, 1], vmin = median - 2*std, vmax = median + 3*std)
    #ax2.imshow(tr[:, :, 1], origin='upper', cmap=plt.cm.bone, norm=LogNorm())
    ax2.title.set_text('Reference')
    ax3 = fig.add_subplot(133)
    ax3.axis('off')
    mean, median, std = sigma_clipped_stats(tr[:, :, 2])
    ax3.imshow(tr[:, :, 2], vmin = median - 2*std, vmax = median + 3*std)
    #ax3.imshow(tr[:, :, 2], origin='upper', cmap=plt.cm.bone)
    ax3.title.set_text('Difference')

    if show_fig:
        plt.show()
    return fig

def get_cutouts(name, username, password):
    """Query kowalski to get the candidate stamps"""
    from penquins import Kowalski

    k = Kowalski(username=username, password=password, verbose=False)

    if type(name) == str:
        list_names = [name]
    elif type(name) == list:
        list_names = name
    else:
        print(f"{name} must be a list or a string")
        return None

    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts",
                   "filter": {
                              'objectId': {'$in': list(list_names)}
                              },
                   "projection": {
                                  "objectId": 1,
                                  "candidate.jd": 1,
                                  "candidate.ra": 1,
                                  "candidate.dec": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.fid": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.programid": 1,
                                  "candidate.field": 1,
                                  "candidate.rcid": 1,
                                  "cutoutScience": 1,
                                  "cutoutTemplate": 1,
                                  "cutoutDifference": 1,
                                  }
                       },
         "kwargs": {"hint": "objectId_1"}
         }

    r = k.query(query=q)

    if r['default']['data'] == []:
        print("No candidates to be checked?")
        return None
    else:
        alerts = r['default']['data']

    return alerts

def get_dust_info(coords):
    """
    Get coordinates as a SkyCoord object
    and returns E(B-V) 
    """
    from dustmaps.planck import PlanckQuery
    from dustmaps.sfd import SFDQuery
    
    planck = PlanckQuery()
    ebv = planck(coords)
    #print('E(B-V) = {:.3f} mag'.format(ebv))

    return ebv


def plot_lc(name, con, cur, forced=True, stack=False,
            plot_alerts=True, save=False, reddening=False,
            plot_cow=True, plot_gfo=True, plot_bulla=True, filtermatch = 'g',
            plot_gw=False, inset=False, tr=None, writecsv=False,
            show_fig=True, program_ids=[1,2,3], grb_time=None,
            plot_alertfp=False):
    '''Plot the light curve of a candidate'''

    color_dict = {'g': 'green', 'r': 'red', 'i': 'y'}
    forced_dict = {'1': 's', '0': 'o'}

    if forced is False:
        plot_alerts = True
        lc = pd.DataFrame(columns=['jd', 'mag', 'mag_unc', 'filter',
                                   'limmag', 'forced', 'programid'])
    # FPZTF
    elif forced is True and plot_alertfp is False:
        if stack is True:
            table = 'lightcurve_stacked'
        elif forced is True:
            table = 'lightcurve_forced'
        lc = pd.read_sql_query(f"SELECT jd, mag, mag_unc, filter, limmag, programid FROM {table} WHERE name = '{name}'", con)
        lc["forced"] = np.ones(len(lc))
        lc['origin'] = ['fpztf' for i in range(0, len(lc))]
    # Alert FP
    else:
        table = 'lightcurve_alertfp'
        lc = pd.read_sql_query(f"SELECT jd, magpsf, sigmapsf, filter, limmag, programid, origin FROM {table} WHERE name = '{name}'", con) 
        forced_idx = np.where(lc['origin']=='alertfp')[0]
        lc["forced"] = [1 if i in forced_idx else 0 for i in range(0, len(lc))]
        #lc_nondet = lc[lc['origin'] == 'magul']
        #lc = lc[lc['origin'] != 'magul']
        lc = lc.rename(columns={'magpsf':'mag', 'sigmapsf':'mag_unc'})
    
    if plot_alerts is True:
        alerts = pd.read_sql_query(f"SELECT jd, magpsf, sigmapsf, filter, programid FROM lightcurve WHERE name = '{name}'", con)

        # Remove the alerts if they are way too many
        ##if len(alerts) > 20:
        ##    print("TOO MANY ALERTS!!!! Cutting to 20...")
        ##    idx = np.random.choice(np.arange(len(alerts)), size=20, replace=False)
        ##    alerts = alerts.iloc[idx]

        d_alerts = {'jd': alerts['jd'],
                    'mag': alerts['magpsf'],
                    'mag_unc': alerts['sigmapsf'],
                    'filter': alerts['filter'],
                    'limmag': np.ones(len(alerts))*99.0,
                    'forced': np.zeros(len(alerts)),
                    'programid': alerts['programid'],
                    'origin': ['alert' for i in range(0, len(alerts))]}
        lc_alerts = pd.DataFrame(data=d_alerts)

        # Remove alerts where forced photometry is present
        if forced is True or plot_alertfp is True:
            jd_forced = list(lc["jd"])
            for j in set(alerts['jd'].values):
                # If the time difference between the alert and any 
                # forced phot is <5s, remove the alert
                if not (lc_alerts.empty) and not (lc.empty):
                    if np.min(np.abs(lc['jd'].values - j)) < 5./60/60/24.:
                        lc_alerts.drop(lc_alerts[lc_alerts['jd'] == j].index,
                                       inplace=True)
        
        lc = pd.concat([lc, lc_alerts], ignore_index=True) # alerts + forced
    
    ## Separate out program ID = 3 if not in progam_ids to plot
    idx = np.where(np.in1d(lc['programid'],program_ids))[0]
    if len(idx) < 2:
        print('Fewer than 2 detections ... continuing.')
        return None

    if len(idx) != len(lc):
        print('Warning: plotting %d/%d of available points' % (len(idx), len(lc)))
    
    lc = lc.iloc[idx]

    # Plot
    if plot_gw is True:
        gw_info = {'ZTF19acbqtue': {'gw_name': 'S190930t',
                                    'gw_time': Time('2019-09-30 14:34:07',
                                    format='iso')}
                   }
        if name in gw_info.keys():
            t0 = gw_info[name]['gw_time'].jd
            xlabel = gw_info[name]['gw_name']
        else:
            t0 = np.min(lc['jd'].values)
            xlabel = Time(t0, format='jd').iso[0:1]
    else:
        if lc.empty:
            print(f"Empty light curve for {name} with forced={forced}, \
stack={stack}, plot_alerts={plot_alerts}")
            return
        if len(np.where(lc['mag'] < 50)[0]) == 0:
            print(f"No detections in light curve for {name} with forced={forced}, stack={stack}, plot_alerts={plot_alerts}")
            return
        t0 = np.min(lc[lc['mag'] < 50]['jd'].values)
        xlabel = Time(t0, format='jd').iso[0:19]

    # Correct for Galactic extinction
    if reddening is True:
        # Get the coordinates
        coords_tbl = pd.read_sql_query(f"SELECT ra, dec FROM candidate \
                                WHERE name = '{name}'", con)
        coords = SkyCoord(ra=coords_tbl["ra"][0].value*u.deg,
                          dec=coords_tbl["dec"][0].value*u.deg)
        ebv = get_dust_info(coords)
        
        """
        # Dict with correction factor to go from ebv to A_lambda
        corr_dict = ["g": 3.3, "r", "i"]
        """

    fig, ax1 = plt.subplots(1, 1, figsize=(9,6))
    ymin, ymax = np.inf, -np.inf

    for f in set(lc['filter']):
        tf = lc[lc['filter'] == f]

        """
        # Correct for the extinction
        if reddening is True:
            corr = ebv * corr_dict[f]
            tf["mag"] = tf["mag"] - 
            tf["limmag"] = tf["limmag"] - 
        """
        tf_det = tf[(tf['mag'] < 50.) & (tf['origin'] != 'magul')]
        tf_ul = tf[(tf['mag'] > 50.) | (tf['origin'] == 'magul')]
        
        for isforced in [0,1]:
            if isforced == 0:
                label = f"{f} alerts"
            else:
                if stack is True:
                    label = f"{f} forced phot stacked"
                elif forced is True and plot_alertfp is False:
                    label = f"{f} forced phot"
                elif plot_alertfp is True:
                    label = f"{f} alert fp"
            tf_det2 = tf_det[tf_det['forced'] == isforced]
            if len(tf_det2) == 0:
                continue
            ax1.errorbar(tf_det2['jd'].values - t0,
                         tf_det2['mag'], yerr=tf_det2['mag_unc'],
                         color=color_dict[f], markeredgecolor='k',
                         fmt=forced_dict[str(int(isforced))], label=label)
            if np.min(tf_det2['mag']) < ymin:
                ymin = np.min(tf_det2['mag'])
            if np.max(tf_det2['mag']) > ymax:
                ymax = np.max(tf_det2['mag'])

        if len(tf_ul) != 0:
            ax1.errorbar(tf_ul['jd'].values - t0, tf_ul['limmag'],
                         markeredgecolor=color_dict[f],
                         markerfacecolor='w', fmt='v')
            plt.plot([],[], 'kv', markeredgecolor='k', markerfacecolor='w',
                     label='upper limits')

            if np.min(tf_det2['limmag']) < ymin:
                ymin = np.min(tf_det2['limmag'])
            if np.max(tf_det2['limmag']) > ymax:
                ymax = np.max(tf_det2['limmag'])

    # Determine the row at which the filter we want
    # to match (filtermatch), peaks
    if plot_cow or plot_gfo:
        peak_row = lc.iloc[lc[lc['mag'] < 50.][lc['filter'] == filtermatch]['mag'].idxmin()]
    
    # Overplot 2018cow in the filters for which there is
    # candidate photometry (limits)
    if plot_cow:
        AT2018cow = pd.read_fwf('../comparison_photometry/AT2018cow_photometry_table.dat',
                                sep = ' ', comment = '#', header = None,\
            names = ['mjd', 'telescope', 'filter', 'mag', 'ABmag',
                     'magerr', 'source'])
        peak_row_cow = AT2018cow.iloc[AT2018cow[AT2018cow["filter"] == peak_row['filter']]['mag'].idxmin()]
        peak_offset = peak_row['mag'] - peak_row_cow['mag']
        mjd_offset = t0 - peak_row['jd'] + peak_row_cow['mjd']
        for f in set(lc['filter']):
            cow_filt = AT2018cow[AT2018cow['filter'] == f]
            ax1.plot(cow_filt['mjd'] - mjd_offset,
                     cow_filt['ABmag'] + peak_offset,
                     color=color_dict[f], linestyle = ':')
            plt.plot([],[], color = 'black', linestyle = ':', label='AT 2018cow')
    # Overplot 2017gfo in the filters for which there is candidate photometry (limits)
    if plot_gfo:
        AT2017gfo = pd.read_csv('../comparison_photometry/AT2017gfo_optical_photometry_smartt+17.txt',
                                sep = ' ', comment = '#',na_values='-')
        peak_row_gfo = AT2017gfo.iloc[[AT2017gfo[peak_row['filter']].idxmin()]]
        peak_offset = peak_row['mag'] - peak_row_gfo[peak_row['filter']]
        mjd_offset = t0 - peak_row['jd'] + peak_row_gfo['Phase']
        for f in set(lc['filter']):
            nanmask = np.isfinite(AT2017gfo[f])
            ax1.plot(AT2017gfo['Phase'][nanmask] - mjd_offset.values,
                     AT2017gfo[f][nanmask] + peak_offset.values,
                     color=color_dict[f], linestyle = '-')
            plt.plot([],[], color = 'black', linestyle = '-',
                     label='AT 2017gfo')
    # Overplot decline rates in g and r from Bulla models    
    if plot_bulla:
        peak_r = lc.iloc[lc.query("mag < 50. & filter == 'r'")['mag'].idxmin()]
        peak_g = lc.iloc[lc.query("mag < 50. & filter == 'g'")['mag'].idxmin()]
        x = np.arange(0,7,0.01)
        ax1.plot(x + peak_r['jd'] - t0, x * 0.3 + peak_r['mag'],
                 linestyle = '--', color = 'red', alpha=0.5)
        ax1.plot(x + peak_g['jd'] - t0, x * 0.5 + peak_g['mag'],
                 linestyle = '--', color = 'green', alpha=0.5)
        plt.plot([],[], color = 'black', linestyle = '--',
                 label='Bulla upper limits')
    if not grb_time is None:
        plt.plot([grb_time - t0, grb_time - t0], [ymin, ymax], 'k--')
    
    plt.gca().invert_yaxis()

    ax1.set_xlabel(f"Days from {xlabel}", fontsize=18)
    ax1.set_ylabel("Magnitude (AB)", fontsize=18)

    # Legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=16)

    ax1.tick_params(axis='both',       # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    labelsize=16)

    # Add the discovery image as inset
    if inset is True:
        # Define the position
        # Reference
        cm = plt.cm.cubehelix
        left, bottom, width, height = [0.51, 0.6, 0.23, 0.23]
        ax2 = fig.add_axes([left, bottom, width, height])
        ax2.axis('off')

        # From Tomas 
        s_max, s_min = [1.001,1.001], [1.001,1.001]
        img_data = tr[:, :, 0]
        imgd_scaled = np.log10(img_data)
        vmax = np.sort(imgd_scaled[28:36,28:36].flatten())[-3]
        npixel = (len(img_data)+1)**2
        imgd_flat = img_data.flatten()
        imgd_scaled[imgd_scaled<0]=np.nanmedian(img_data)
        v_onesig = np.log10(np.nanmedian(img_data) - mad_std(imgd_flat[np.isnan(imgd_flat)==False])*1.4826)
        vmin= max(v_onesig, np.nanmin(imgd_scaled))
        ax2.imshow(imgd_scaled, cmap=cm, vmax=s_max[0]*vmax, vmin=vmin*s_min[0])
        
        '''
        # Subtraction
        left, bottom, width, height = [0.68, 0.6, 0.23, 0.23]
        ax3 = fig.add_axes([left, bottom, width, height])
        ax3.axis('off')
        mean, median, std = sigma_clipped_stats(tr[:, :, 2])
        ax3.imshow(tr[:, :, 2], vmin = median - 2*std, vmax = median + 3*std)
        '''

        # Reference
        left, bottom, width, height = [0.68, 0.6, 0.23, 0.23]
        ax3 = fig.add_axes([left, bottom, width, height])
        ax3.axis('off')
        s_max, s_min = [1.001,1.001], [1.001,1.001]
        img_data = tr[:, :, 1]
        imgd_scaled = np.log10(img_data)
        vmax = np.sort(imgd_scaled[28:36,28:36].flatten())[-3]
        npixel = (len(img_data)+1)**2
        imgd_flat = img_data.flatten()
        imgd_scaled[imgd_scaled<0]=np.nanmedian(img_data)
        v_onesig = np.log10(np.nanmedian(img_data) - mad_std(imgd_flat[np.isnan(imgd_flat)==False])*1.4826)
        vmin= max(v_onesig, np.nanmin(imgd_scaled))
        ax3.imshow(imgd_scaled, cmap=cm, vmax=s_max[0]*vmax, vmin=vmin*s_min[0])     

        # Add candidate name
        ax2.text(0.8, 0.57, name, color='k', fontsize=16, ha='center', va='center', transform=ax1.transAxes)
    
    if save is True or writecsv is True:
        if forced is True:
            forcedbool = 1
        else:
            forcedbool = 0
        if stack is True:
            stackbool = 1
        else:
            stackbool = 0
        plt.savefig(f"lc_{name}_forced{forcedbool}_stacked{stackbool}.png")
        if writecsv is True:
            lc.to_csv(f"lightcurves/lc_{name}_forced{forcedbool}_stacked{stackbool}.csv")
    if show_fig:
        plt.show()
    return fig

def get_xmatch_clu_glade(list_names, con, cur):
    # Filter for match in either CLU or GLADE 
    str_names = "'" + "', '".join(list_names) + "'"
    galaxy_match_clu = pd.read_sql_query(f"SELECT name, \
    clu_ra, clu_dec, clu_distmpc, clu_z, clu_dist_kpc, clu_sep_arcsec \
    FROM crossmatch \
    WHERE name in ({str_names}) and clu_ra is not NULL", con)

    galaxy_match_glade = pd.read_sql_query(f"SELECT name, \
    glade_ra, glade_dec, glade_dist_mpc, glade_z, glade_dist_kpc, glade_sep_arcsec \
    FROM crossmatch \
    WHERE name in ({str_names}) and glade_ra is not NULL", con)

    return galaxy_match_clu, galaxy_match_glade


def get_bgal_ebv(list_names, con, cur):
    """
    Get information from the db about the Galactic latitude and extinction
    """

    str_names = "'" + "', '".join(list_names) + "'"
    bgal_ebv = pd.read_sql_query(f"SELECT name, \
    ra, dec, b_gal, ebv \
    FROM candidate \
    WHERE name in ({str_names})", con)
    
    return bgal_ebv

def agn_b_scores(name,username,password,colors=False):
    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
             "catalog": 'ZTF_alerts',
             "filter": {"objectId": name},
             "projection": {"_id": 0, "cutoutScience": 0, "cutoutTemplate": 0, "cutoutDifference": 0},
         }
         }
    r = k.query(query=q)
    alerts = r['default']['data']
    ra,dec = alerts[0]['candidate']['ra'],alerts[0]['candidate']['dec']
    
    cc = SkyCoord(ra,dec,unit=(u.deg,u.deg))
    table = Irsa.query_region(coordinates=cc, catalog="allwise_p3as_psd", spatial="Cone",radius=2 * u.arcsec)

    # AGN WISE
    if len(table['w1mpro'])==0:
            agn=False
            temp_points = 6
    else:
            w1,w1_err,w2,w2_err,w3,w3_err = table['w1mpro'],table['w1sigmpro'],table['w2mpro'],table['w2sigmpro'],table['w3mpro'],table['w3sigmpro']
            if w1-w2>0.8+0.1 and w2_err<0.5 and w1_err<0.5:
                agn = True
                temp_points = -2
            elif w2-w3>2.5+0.1 and w2_err<0.5 and w3_err<0.5:
                agn = True
                temp_points = -2
            elif w1-w2>0.8 and w2_err<0.5 and w1_err<0.5:
                agn = True
                temp_points = 0
            elif w2-w3>2.5 and w2_err<0.5 and w3_err<0.5:
                agn = True
                temp_points = 0
            elif w1-w2>0.8-0.2 and w2_err<0.5 and w1_err<0.5:
                agn = False
                temp_points = 2
            elif w2-w3>2.5-0.3 and w2_err<0.5 and w3_err<0.5:
                agn = False
                temp_points = 2
            elif w1-w2>0.8-0.5 and w2_err<0.5 and w1_err<0.5:
                agn = False
                temp_points = 4
            elif w2-w3>2.5-0.5 and w2_err<0.5 and w3_err<0.5:
                agn = False
                temp_points = 4

            else:
                agn=False
                temp_points = 6
    # low b
    if np.abs(cc.galactic.b.value)<15:
        b_temp_points =  -10 
    else:
        b_temp_points = 0
    
    if colors:
        return temp_points,agn,[w1-w2,w2-w3]
    else:
        return temp_points,b_temp_points  
