#!/usr/bin/env python

'''
Query Kowalski searching for transients
given a set of constraints.

'''
import os
import numpy as np
import json
import pdb
import pickle

from astropy.time import Time
from astropy.io import ascii
#from astropy.io import fits
from astropy.table import Table
import pandas as pd
import requests
import psycopg2
import corner

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
import healpy as hp

import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
plt.rcParams['xtick.labelsize']=30
plt.rcParams['ytick.labelsize']=30

import pymultinest

from gwemlightcurves.sampler import *
from gwemlightcurves.KNModels import KNTable
from gwemlightcurves.sampler import run
from gwemlightcurves import __version__
from gwemlightcurves import lightcurve_utils, ztf_utils, Global

def do_knfit(df, outputDir = './knfit',
             errorbudget=1.0,
             n_live_points = 100,
             evidence_tolerance = 0.5):

    if len(df.index) < 2:
        print('Not enough data for candidate %s... continuing.', df.candname.values[0])
        return

    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)

    candname = df.name.values[0].decode()
    filters = list(set(df.filtname))

    mint = 0.0
    maxt = 7.0
    dt = 0.05
    tt = np.arange(mint,maxt,dt) 
    doWaveformExtrapolate = True

    ZPRange = 5.0
    T0Range = 2.0

    ModelPath = "../gwemlightcurves/output/svdmodels/"

    T0 = np.inf
    mag_min = np.inf

    data_out = {}
    for index, row in df.iterrows():
        filt = row.filtname
        mag = row.magpsf
        dmag = row.sigmapsf
        mjd = Time(row.jd, format='jd').mjd
        magzp = row.magzpsci

        if 99.0 == mag:
            mag = magzp
            dmag = np.inf
        else:
            T0 = np.min([T0, mjd])
            mag_min = np.min([mag_min, mag])

        if not filt in data_out:
            data_out[filt] = np.empty((0,3), float)
        data_out[filt] = np.append(data_out[filt],np.array([[mjd,mag,dmag]]),axis=0)    

    distmod = mag_min - -16
    distance = 10**((distmod/5.0) + 1.0) / 1e6

    for ii,key in enumerate(data_out.keys()):
        if key == "t":
            continue
        else:
            data_out[key][:,0] = data_out[key][:,0] - T0
            data_out[key][:,1] = data_out[key][:,1] - 5*(np.log10(distance*1e6) - 1)

    for ii,key in enumerate(data_out.keys()):
        if key == "t":
            continue
        else:
            idxs = np.intersect1d(np.where(data_out[key][:,0]>=mint)[0],np.where(data_out[key][:,0]<=maxt)[0])
            data_out[key] = data_out[key][idxs,:]

    for ii,key in enumerate(data_out.keys()):
        idxs = np.where(~np.isnan(data_out[key][:,2]))[0]
        if key == "t":
            continue
        else:
            data_out[key] = data_out[key][idxs,:]

    for ii,key in enumerate(data_out.keys()):
        if not key in filters:
            del data_out[key]

    for ii,key in enumerate(data_out.keys()):
        if ii == 0:
            samples = data_out[key].copy()
        else:
            samples = np.vstack((samples,data_out[key].copy()))

    idx = np.argmin(samples[:,0])
    samples = samples[idx,:]

    Global.data_out = data_out
    Global.errorbudget = errorbudget
    Global.ZPRange = ZPRange
    Global.T0Range = T0Range
    Global.doLightcurves = 1
    Global.filters = filters
    Global.doWaveformExtrapolate = doWaveformExtrapolate

    modelfile = os.path.join(ModelPath,'Bu2019inc_mag.pkl')
    with open(modelfile, 'rb') as handle:
        svd_mag_model = pickle.load(handle)
    Global.svd_mag_model = svd_mag_model

    modelfile = os.path.join(ModelPath,'Bu2019inc_lbol.pkl')
    with open(modelfile, 'rb') as handle:
        svd_lbol_model = pickle.load(handle)
    Global.svd_lbol_model = svd_lbol_model

    plotDir = os.path.join(outputDir,candname)
    if not os.path.isdir(plotDir):
        os.makedirs(plotDir)

    max_iter = -1
    best = []

    parameters = ["t0","mej","phi","theta","zp"]
    labels = [r"$T_0$",r"${\rm log}_{10} (M_{\rm ej})$",r"$\Phi$",r"$\Theta$","ZP"]
    n_params = len(parameters)
    pymultinest.run(myloglike_Bu2019inc_ejecta, myprior_Bu2019inc_ejecta, n_params, importance_nested_sampling = False, resume = True, verbose = True, sampling_efficiency = 'parameter', n_live_points = n_live_points, outputfiles_basename='%s/2-'%plotDir, evidence_tolerance = evidence_tolerance, multimodal = False, max_iter = max_iter)

    multifile = lightcurve_utils.get_post_file(plotDir)
    data = np.loadtxt(multifile)

    t0, mej, phi, theta, zp, loglikelihood = data[:,0], 10**data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]
    idx = np.argmax(loglikelihood)
    t0_best, mej_best, phi_best, theta_best, zp_best = data[idx,0], 10**data[idx,1], data[idx,2], data[idx,3], data[idx,4]
    zp_mu, zp_std = 0.0, Global.ZPRange
    zp_best = scipy.stats.norm(zp_mu, zp_std).ppf(zp_best)
    tmag, lbol, mag = Bu2019inc_model_ejecta(mej_best,phi_best,theta_best)

    pcklFile = os.path.join(plotDir,"data.pkl")
    f = open(pcklFile, 'wb')
    pickle.dump((data_out, data, tmag, lbol, mag,
                 t0_best, zp_best, n_params, labels, best), f)
    f.close()

    title_fontsize = 30
    label_fontsize = 30

    plotName = "%s/corner.pdf"%(plotDir)
    figure = corner.corner(data[:,:-1], labels=labels,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True, title_kwargs={"fontsize": title_fontsize},
                           label_kwargs={"fontsize": label_fontsize}, title_fmt=".2f",
                           smooth=3,
                           color="coral")
    figure.set_size_inches(14.0,14.0)
    plt.savefig(plotName)
    plt.close()

    tmag = tmag + t0_best
    #colors=cm.rainbow(np.linspace(0,1,len(filters)))
    colors=cm.Spectral(np.linspace(0,1,len(filters)))[::-1]

    color2 = 'coral'
    color1 = 'cornflowerblue'
    
    plotName = "%s/models_panels.pdf"%(plotDir)
    #plt.figure(figsize=(20,18))
    plt.figure(figsize=(20,28))
    
    cnt = 0
    for filt, color in zip(filters,colors):
        cnt = cnt+1
        if cnt == 1:
            ax1 = plt.subplot(len(filters),1,cnt)
        else:
            ax2 = plt.subplot(len(filters),1,cnt,sharex=ax1,sharey=ax1)
    
        if not filt in data_out: continue
        samples = data_out[filt]
        t, y, sigma_y = samples[:,0], samples[:,1], samples[:,2]
        idx = np.where(~np.isnan(y))[0]
        t, y, sigma_y = t[idx], y[idx], sigma_y[idx]
        if len(t) == 0: continue
    
        idx = np.where(np.isfinite(sigma_y))[0]
        plt.errorbar(t[idx],y[idx],sigma_y[idx],fmt='o',c=color, markersize=16, label='%s-band'%filt)
    
        idx = np.where(~np.isfinite(sigma_y))[0]
        plt.errorbar(t[idx],y[idx],sigma_y[idx],fmt='v',c=color, markersize=16)
    
        magave = lightcurve_utils.get_mag(mag,filt)
        ii = np.where(~np.isnan(magave))[0]
        f = interp.interp1d(tmag[ii], magave[ii], fill_value='extrapolate')
        maginterp = f(tt)
        #plt.plot(tt,maginterp+zp_best,'--',c=color,linewidth=3)
        #plt.fill_between(tt,maginterp+zp_best-errorbudget,maginterp+zp_best+errorbudget,facecolor=color,alpha=0.2)
    
        plt.plot(tt,maginterp+zp_best,'--',c=color2,linewidth=3)
        plt.fill_between(tt,maginterp+zp_best-errorbudget,maginterp+zp_best+errorbudget,facecolor=color2,alpha=0.2)
    
        plt.ylabel('%s'%filt,fontsize=48,rotation=0,labelpad=40)
    
        plt.xlim([0.0, 10.0])
        plt.ylim([-22.0,-8.0])
        plt.gca().invert_yaxis()
        plt.grid()
    
        if cnt == 1:
            ax1.set_yticks([-22,-18,-14,-10])
            plt.setp(ax1.get_xticklabels(), visible=False)
            #l = plt.legend(loc="upper right",prop={'size':36},numpoints=1,shadow=True, fancybox=True)
        elif not cnt == len(filters):
            plt.setp(ax2.get_xticklabels(), visible=False)
        plt.xticks(fontsize=36)
        plt.yticks(fontsize=36)
    
    ax1.set_zorder(1)
    plt.xlabel('Time [days]',fontsize=48)
    plt.tight_layout()
    plt.savefig(plotName)
    plt.close()

    return

