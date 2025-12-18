"""
author: sarah
Create tables with replicate schema in 2025 database with additional columns as needed.
created 08/08/2025
This is a record of the tables and columns added.

candidate
    - add columns: index_rise/fade_g/r/i
                   index_nondet_rise_g/r/i
lightcurve
    - add columns: 'origin'  ('prv' or 'alert')
lightcurve_forced
lightcurve_stacked -- need to choose either alert forced photometry or ForcePhotZTF*
crossmatch - update if needed?

------ NEW: ------
*lightcurve_stacked_alertfp -- if decide to do both
lightcurve_alertfp -- alert forced photometry data

UPDATED data tables 10/31/2025. Dropped all tables (empty) and re-made with functions below.

"""
import glob
import pdb
from socket import gethostname

import numpy as np
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, vstack, unique
from astroquery.vizier import Vizier
from astropy.cosmology import Planck15 as cosmo
from astropy.io.misc.hdf5 import read_table_hdf5
import pandas as pd
import psycopg2


def connect_database(update_database=False, path_secrets_db='db_access.csv',
                     dbname='db_kn_2025_admin'):
    """
    Establish a connection to the psql database

    ----
    Parameters

    update_database bool
        if True, the connection will be established as admin
        for w+r privileges
    path_secrets_db str
        path to the CSV file with the db access information    
    dbname
        specifies which database and privledges (user or admin)
    ----
    Returns

    con, cur
        connection and cursor for the interaction with the db
    """
    if update_database is True:
        confirm = str(input(f"You are accessing {dbname} with admin privledges. \nContinue (y/n)?"))
        if confirm != 'y':
            print("Exiting.")
            quit()
    # Read the secrets
    info = ascii.read(path_secrets_db, format='csv')
    # Admin access only if writing is required
    if update_database is True and dbname is None:
        info_db = info[info['db'] == 'db_kn_rt_admin']
    elif update_database is False and dbname is None:
        info_db = info[info['db'] == 'db_kn_rt_user']
    elif update_database is True and dbname is not None:
        info_db = info[info['db'] == dbname]
    elif update_database is False and dbname is not None:
        info_db = info[info['db'] == dbname]
    if gethostname() == 'usnik':
        host = 'localhost'
    else:
        host = info_db['host'][0]
    db_kn = f"host={host} dbname={info_db['dbname'][0]} \
port={info_db['port'][0]} user={info_db['user'][0]} \
password={info_db['password'][0]}"
    con = psycopg2.connect(db_kn)
    cur = con.cursor()

    return con, cur

def alter_columns(con,cur):
    """
    fixing my mistakes :/
    """
    #cur.execute("ALTER TABLE candidate DROP COLUMN index_nondet_rise")
    #cur.execute("ALTER TABLE candidate ADD index_nondet_rise_g FLOAT")
    #cur.execute("ALTER TABLE candidate ADD index_nondet_rise_r FLOAT")
    #cur.execute("ALTER TABLE candidate ADD index_nondet_rise_i FLOAT")
    #cur.execute("ALTER TABLE candidate ALTER COLUMN ra TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE candidate ALTER COLUMN dec TYPE DECIMAL(10,7)")
    #fixlist = ['clu_ra', 'clu_dec', 'glade_ra', 'glade_dec', 'gaia_ra', 'gaia_dec', 'ls_ra', 'ls_dec']
    #for col in fixlist:
    #    cur.execute(f"ALTER TABLE crossmatch ALTER COLUMN {col} TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve ALTER COLUMN ra TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve ALTER COLUMN dec TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve_forced ALTER COLUMN ra TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve_forced ALTER COLUMN dec TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve_alertfp ALTER COLUMN ra TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve_alertfp ALTER COLUMN dec TYPE DECIMAL(10,7)")
    #cur.execute("ALTER TABLE lightcurve ALTER COLUMN pid TYPE BIGINT")
    #cur.execute("ALTER TABLE lightcurve_forced ALTER COLUMN pid TYPE BIGINT")
    #cur.execute("ALTER TABLE lightcurve_alertfp ALTER COLUMN pid TYPE BIGINT")
    #cur.execute ("DROP TABLE lightcurve_forced") # deleted and re-made with correct columns
    
    #con.commit()
    return
# ------------- CANDIDATE --------------------------------------------------
def create_table_candidate(con, cur):
    # Table for the candidates
    cur.execute("CREATE TABLE candidate(name TEXT NOT NULL PRIMARY KEY, \
                ra NUMERIC(9,6), dec NUMERIC(9,6), \
                sgscore1 FLOAT, sgscore2 FLOAT, sgscore3 FLOAT, \
                distpsnr1 FLOAT, distpsnr2 FLOAT, distpsnr3 FLOAT, \
                index_rise FLOAT, index_fade FLOAT, \
                gaia_match INT, gaia_stellar INT, \
                clu_match INT, glade_match INT, comment TEXT)")
    cur.execute("ALTER TABLE candidate ADD index_fade_g FLOAT") 
    cur.execute("ALTER TABLE candidate ADD index_fade_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_i FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_i FLOAT")
    cur.execute("ALTER TABLE candidate ADD mindet3 INT") # min number of detections
    cur.execute("ALTER TABLE candidate ADD mindet2 INT")
    
    cur.execute("ALTER TABLE candidate ADD index_fade_stack_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_stack_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_stack_i FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_stack_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_stack_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_stack_i FLOAT")
    
    cur.execute("ALTER TABLE candidate ADD duration_g FLOAT") 
    cur.execute("ALTER TABLE candidate ADD duration_r FLOAT") 
    cur.execute("ALTER TABLE candidate ADD duration_i FLOAT") 
    cur.execute("ALTER TABLE candidate ADD duration_tot FLOAT")
    
    cur.execute("ALTER TABLE candidate ADD index_fade_forced_g FLOAT") 
    cur.execute("ALTER TABLE candidate ADD index_fade_forced_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_forced_i FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_forced_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_forced_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_forced_i FLOAT")
    
    cur.execute("ALTER TABLE candidate ADD index_fade_alertfp_g FLOAT") 
    cur.execute("ALTER TABLE candidate ADD index_fade_alertfp_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_alertfp_i FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_alertfp_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_alertfp_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_alertfp_i FLOAT")
    
    # keep track of candidates with deep upper limits:
    cur.execute("ALTER TABLE candidate ADD index_nondet_rise_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_nondet_rise_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_nondet_rise_i FLOAT")
    
    cur.execute("ALTER TABLE candidate ADD hard_reject INT")
    
    cur.execute("ALTER TABLE candidate ADD max_days_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD max_days_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD max_days_i FLOAT")
    
    cur.execute("ALTER TABLE candidate ADD o3 INT")
    
    cur.execute("ALTER TABLE candidate ADD before_trigger INT")
    cur.execute("ALTER TABLE candidate ADD ebv FLOAT")
    cur.execute("ALTER TABLE candidate ADD b_gal FLOAT")
    
    # commit the changes
    con.commit()
    
# ---------- CROSSMATCH ---------------------------------------------------
def create_table_crossmatch(con, cur):
    cur.execute("CREATE TABLE crossmatch(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, clu_id INT, clu_ra NUMERIC(9,6), clu_dec NUMERIC(9,6), \
                clu_z FLOAT, clu_zerr FLOAT, clu_distmpc FLOAT, \
                clu_mstar FLOAT, clu_sfr_fuv FLOAT, clu_sfr_ha FLOAT, \
                clu_w1mpro FLOAT, clu_w1sigmpro FLOAT, clu_w2mpro FLOAT, \
                clu_w2sigmpro FLOAT, clu_w3mpro FLOAT, clu_w3sigmpro FLOAT, \
                clu_w4mpro FLOAT, clu_w4sigmpro FLOAT, clu_type_ned TEXT, \
                clu_a FLOAT, clu_b2a FLOAT, clu_dist_kpc FLOAT, clu_sep_arcsec FLOAT, \
                glade_sep_arcsec FLOAT, glade_name_gwgc TEXT, \
                glade_name_sdss TEXT, glade_name_hl TEXT,\
                glade_name_2mass TEXT, glade_ra NUMERIC(9,6), glade_dec NUMERIC(9,6), \
                glade_z FLOAT, glade_dist_mpc FLOAT, glade_dist_mpc_err FLOAT,\
                glade_dist_kpc FLOAT, glade_bmag FLOAT, glade_bmagerr FLOAT, \
                gaia_plx FLOAT, gaia_bmag FLOAT, gaia_bmagerr FLOAT, \
                gaia_ra NUMERIC(9,6), gaia_dec NUMERIC(9,6), gaia_sep_arcsec FLOAT,\
                ls_ra NUMERIC(9,6), ls_dec NUMERIC(9,6), ls_sep_arcsec FLOAT, ls_z_spec FLOAT, \
                ls_z_phot_median FLOAT, ls_z_phot_std FLOAT, \
                ls_photoz_checked INT,  ls_type TEXT, \
                ls_z_phot_l95 FLOAT, ls_z_phot_u95 FLOAT, ls_fluxz FLOAT)")
                
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_spec FLOAT")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_phot_median FLOAT")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_phot_u95 FLOAT")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_phot_l95 FLOAT")
    
    # commit the changes
    con.commit()


# ---------- LIGHTCURVE ---------------------------------------------------
def create_table_lightcurve(con, cur):
    """
    alert+prv photometry
    """
    cur.execute("CREATE TABLE lightcurve(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, ra NUMERIC(9,6), dec NUMERIC(9,6), jd FLOAT, \
                magpsf FLOAT, sigmapsf FLOAT, filter TEXT, \
                magzpsci FLOAT, magzpsciunc FLOAT, programid INT, \
                field INT, rcid INT, pid BIGINT)")
    # Add column: origin ('alert' or 'prv')        
    cur.execute("ALTER TABLE lightcurve ADD origin TEXT")
    # commit the changes
    con.commit()
    
# ---------- LIGHTCURVE_FORCED --------------------------------------------
def create_table_lightcurve_forced(con, cur):

    # ForcePhotZTF forced photometry
    #cur.execute("CREATE TABLE lightcurve_forced(id INTEGER NOT NULL PRIMARY KEY, \
                #name TEXT, ra NUMERIC(9,6), dec NUMERIC(9,6), jd FLOAT, \
                #magpsf FLOAT, sigmapsf FLOAT, filter TEXT, \
                #magzpsci FLOAT, magzpsciunc FLOAT, programid INT, \
                #field INT, rcid BIGINT, pid BIGINT)")
    cur.execute("CREATE TABLE lightcurve_forced(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, jd FLOAT, filter TEXT, seeing FLOAT, gain FLOAT, zp FLOAT, \
                ezp FLOAT, programid INT, field INT, ccdid INT, qid INT, filterid INT, \
                moonra NUMERIC(9,6), moondec NUMERIC(9,6), moonillf FLOAT, moonphase FLOAT, \
                airmass FLOAT, nbad FLOAT, nbadbkg FLOAT, bkgstd FLOAT, bkgmed FLOAT, diffimgname TEXT, \
                psfimgname TEXT, flux_maxlike FLOAT, flux_maxlike_unc FLOAT, chi2_nu FLOAT, \
                fratio FLOAT, fratio_unc FLOAT, mag FLOAT, mag_unc FLOAT, limmag FLOAT);")
    # commit the changes
    con.commit()
    
# --------- LIGHTCURVE_STACKED ---------------------------------------------
def create_table_lightcurve_stacked(con, cur):
    cur.execute("CREATE TABLE lightcurve_stacked(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, jd FLOAT, \
                flux FLOAT, flux_unc FLOAT, \
                mag FLOAT, mag_unc FLOAT, limmag FLOAT, filter TEXT, \
                zp FLOAT, ezp FLOAT, programid INT, \
                field INT, ccdid INT, qid INT)")
    # commit the changes
    con.commit()

# --------- LIGHTCURVE_ALERTFP ----------------------------------------------
def create_table_lightcurve_alertfp(con,cur):
    # Forced photometry from the ZTF alerts
    cur.execute("CREATE TABLE lightcurve_alertfp(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, ra NUMERIC(9,6), dec NUMERIC(9,6), jd FLOAT, \
                magpsf FLOAT, sigmapsf FLOAT, filter TEXT, \
                magzpsci FLOAT, magzpsciunc FLOAT, programid INT, \
                field INT, rcid INT, pid BIGINT, origin TEXT, \
                limmag FLOAT, snr FLOAT, forcediffimflux FLOAT, forcediffimfluxunc FLOAT)")
    # commit changes
    con.commit()