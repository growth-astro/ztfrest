import gc
import glob
import subprocess

import numpy as np
from astropy.io import ascii

from ForcePhotZTF.force_lc import get_cutout_data, download_images_diffpsf_refdiff
from ForcePhotZTF.force_mcmc import get_forced_phot_mcmc
from ForcePhotZTF.force_maxlike import get_forced_phot_maaxlike
from ForcePhotZTF.refine_lc import get_recerence_jds


import multiprocessing as mp
from multiprocessing import cpu_count
from tqdm import tqdm


def do_forcephot(targetdir, name, ra, dec, start_jd, end_jd,
                 ncpu=1, limit=100, do_mcmc=False,
                 verbose=False):
    """Perform forced photometry"""

    print(f">>> Forced photometry for target {name}")
    print(">>> Downloading images..")
    under_limit = download_images_diffpsf_refdiff(targetdir, ra, dec,
                                                  start_jd=start_jd,
                                                  end_jd=end_jd, limit=limit)
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


def trigger_forced_photometry(t0, targetdir_base,
                              daydelta_before=7., daydelta_after=14.):
    """
    Trigger forced photometry using ForcePhotZTF (Yao et al., 2019)
    for given ZTF transients

    ---
    Parameters

    t0 astropy.table
        table with alerts photometry for the transient candidates
        the table must contain information such as the name of
        the transient and its equatorial coordinates (in degrees)

    targetdir_base str
        directory hosting the forced photometry data products

    daydelta_before float
        number of days before the last detection that will
        be considered for the forced photometry

    daydelta_after float
        number of days after the last detection that will
        be considered for the forced photometry

    ---
    Returns

    The function generates forced photometry light curves

    success list
        list of targets with successful forced photometry generation
    problematic list
        list of targets with problematic forced photometry generation
    """

    print(f"There are a total of {len(set(t0['name']))} candidates considered \
for forced photometry")

    # How many CPUs are available and will be used?
    ncpu_avail = cpu_count()
    print(f"{ncpu_avail} available CPUs")
    if ncpu_avail == 0:
        print("ERROR: No CPU available for forced photometry!")
        print("Exiting...")
        exit()
    elif ncpu_avail <= 4:
        print(f"WARNING: Only {ncpu_avail} CPUs available, using all of them")
        ncpu = ncpu_avail
    else:
        ncpu = ncpu_avail - 4
    print(f"{ncpu} CPUs will be used")
    
    # List of candidate names
    list_names = set(t0['name'])

    # Parallelize the light curve generation
    pool = mp.Pool(processes=ncpu)
    process_list =[]

    for name in list_names:
        process_list.append(pool.apply_async(main_maxlike,
                                             args=(name, t0,
                                                   targetdir_base,
                                                   daydelta_before,
                                                   daydelta_after,)))
    # Add a progress bar
    process_list = tqdm(process_list)
    results = [p.get() for p in process_list]
    pool.close()

    success = [n[0] for n in results if n[1] == 'success']
    problematic = [n[0] for n in results if n[1] == 'problematic']

    print(">>> FORCED PHOTOMETRY: ALL FINISHED")
    print("Total targets:", len(set(t0['name'])))
    print(f"Success {len(success)}:", success)
    print(f"problematic {len(problematic)}:", problematic)

    return success, problematic


if __name__ == "__main__":

    # Define the target directory
    targetdir_base = '/data/ia/paper_kn_ZTF/lc_2019_May_2020_Jan_ebv01_forced/'
    daydelta = 1.

    t0 = ascii.read('lc_results_2y_min0015_max6_ndethist2_CLUnew.csv', format='csv')

    print(f"There are a total of {len(set(t0['name']))} candidates considered")

    ncpu_avail = cpu_count()
    print("{0} available CPUs".format(ncpu_avail))
    ncpu = ncpu_avail - 4
    print("{0} used CPUs".format(ncpu))
    
    list_names = set(t0['name'])

    # Parallelize the light curve generation
    pool = mp.Pool(processes = ncpu)
    process_list =[]

    for name in list_names:
        process_list.append(pool.apply_async(main_maxlike, args=(name,t0, 
                                                                 targetdir_base,
                                                                 daydelta_before,
                                                                 daydelta_after,)))
    # Add a progress bar
    process_list = tqdm(process_list)
    results = [p.get() for p in process_list]
    pool.close()

    success = [n[0] for n in results if n[1] == 'success']
    problematic = [n[0] for n in results if n[1] == 'problematic']

    print(">>> ALL FINISHED")
    print("Total targets:", len(set(t0['name'])))
    print(f"Success {len(success)}:", success)
    print(f"problematic {len(problematic)}:", problematic)

