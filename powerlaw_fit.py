
import os, sys
from astropy.io import ascii
import numpy as np
from astropy.time import Time
import scipy.stats
import optparse
import pymultinest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# Useful functions for the power-law fit
from scipy.optimize import curve_fit
from astropy.table import Table, Row, Column

import corner

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()
    parser.add_option("--doPlots",  action="store_true", default=False)
    parser.add_option("-p","--plotDir",default="plots")

    parser.add_option("-l","--lightcurve",default="data/ZTF21aahifke_lc_all_frank.csv")

    opts, args = parser.parse_args()

    return opts

def plaw(t, f0, alpha=1.0, t0=None):
    """
    Power law function:
    f / f0 = (t - t0) ** -alpha
    """
    #global t0
    return f0 * (t - t0) ** -alpha

def pmag(t, f0, alpha=1.0, t0=None):
    """
    Use indices from plaw fit to calculate magnitudes
    -2.5 log(f/f0) = 2.5 * alpha * log(dt)
    m = -2.5 log(f)
    m + 2.5 log(f0) = 2.5 * alpha * log(dt)
    """
    #global t0
    return 2.5 * alpha * np.log10(t - t0) - 2.5 * np.log10(f0)

def eflux(emag, flux):
    """
    Error propoagation:
        m - m0 = dm = 2.5 log(1 + df / f0)
        10.0 ** (dm / 2.5) = 1 + df / f0
        df = f0 * (10.0 ** (dm / 2.5) - 1)
    """
    return flux * (10.0 ** (emag / 2.5) - 1)

def myprior(cube, ndim, nparams):

        cube[0] = cube[0]*(tmax - tmin) + tmin
        cube[1] = cube[1]*3.0
        cube[2] = cube[2]*20.0 - 10.0

def myloglike(cube, ndim, nparams):
    t0fit = cube[0]
    alpha = cube[1]
    f0 = 10**cube[2]

    mod = pmag(t, f0, alpha=alpha, t0=t0fit)
    idx1 = np.where(~np.isnan(y))[0]
    idx2 = np.where(np.isnan(y))[0]

    chisq = -(1/2)*np.sum(((y[idx1] - mod[idx1])**2.)/(dy[idx1]**2.))/(len(t[idx1]) - len(parameters))

    #gaussprobvals = np.sum(np.log(1-scipy.stats.norm.cdf(dy[idx2], mod[idx2], 0.1)))

    return chisq
    
# Parse command line
opts = parse_commandline()
baseplotDir = opts.plotDir
if not os.path.isdir(baseplotDir):
    os.makedirs(baseplotDir)

lc = ascii.read(opts.lightcurve, format='csv')
parameters = ["t0", "alpha", "f0"]
n_params = len(parameters)

n_live_points = 1000
evidence_tolerance = 0.1
max_iter = -1

#for f in ["g", "r", "i"]:
for f in ["r"]:
    tab = lc[lc['fit'] == 1]
    tab = tab[tab['filter'] == f]
    #tab = tab[tab['mag'] < 50]

    mjd = Time(tab["jd"], format='jd').mjd
    tab.rename_column("mag_unc", "e_mag")
    tab['mjd'] = mjd
    # Select data
    tab.sort('mjd')

    idx = np.where(tab['mag'] > 50)[0]
    tab['mag'][idx] = np.nan
    tab['e_mag'][idx] = tab['limmag'][idx]

    # T0
    tmax = np.min(tab[tab['mag'] < 50]['mjd'])
    tmin =  np.max(tab[np.isnan(tab['mag'])]['mjd'])

    # Some checks
    print("Data that will be used for the fit:")
    print("MJD, mag, instrument")
    for l in tab:
        print(l['mjd'], l['mag'], l['instrument'])
        print("---")

    plotDir = os.path.join(baseplotDir, f)
    if not os.path.isdir(plotDir):
        os.makedirs(plotDir)

    t = tab['mjd']
    y = tab['mag']
    dy = tab["e_mag"]

    pymultinest.run(myloglike, myprior, n_params, importance_nested_sampling = False, resume = True, verbose = True, sampling_efficiency = 'parameter', n_live_points = n_live_points, outputfiles_basename='%s/2-'%plotDir, evidence_tolerance = evidence_tolerance, multimodal = False, max_iter = max_iter)
   
    multifile = os.path.join(plotDir,'2-post_equal_weights.dat')
    data = np.loadtxt(multifile)

    labels = [f"$t_0$", f"$\\alpha$", f"$\\log_{10} f_0$"]
    plotName = "%s/corner.pdf"%(plotDir)
    figure = corner.corner(data[:,:-1], labels=labels,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 24},
                       label_kwargs={"fontsize": 24}, title_fmt=".2f",
                       smooth=3,
                       color="coral")
    figure.set_size_inches(14.0,14.0)
    plt.savefig(plotName)
    plt.close()

    idx = np.argmax(data[:,-1])
    t_0, alpha, f0 = data[idx,:-1]
    tt = np.linspace(np.min(t), np.max(t), 1000)
    mod = pmag(tt, 10**f0, alpha=alpha, t0=t_0)

    plotName = "%s/lightcurve.pdf"%(plotDir)
    plt.figure()
    plt.errorbar(t-tmin,y,dy,fmt='o',c='k')
    plt.plot(tt-tmin, mod,'k--',linewidth=2)
    plt.xlabel('Time [days] [t0 = %.5f]' % tmin,fontsize=24)
    plt.ylabel('Magnitude',fontsize=24)
    plt.grid()
    plt.gca().invert_yaxis()
    plt.savefig(plotName)
    plt.close()

