'''
Functions to connect to, manage, and write to kn database.
note: all create_<tablename> functions are for old database, see new_db_2025.py
for current create_<tablename>, column names, data types, etc.
'''

import glob
import pdb
from socket import gethostname

import numpy as np
np.set_printoptions(legacy="1.25") # for compatibility with psycopg2

from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, vstack, unique
from astroquery.vizier import Vizier
from astropy.cosmology import Planck15 as cosmo
from astropy.io.misc.hdf5 import read_table_hdf5
import pandas as pd
import psycopg2
#from sqlalchemy import create_engine

# the default connection is the 2025 database
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

    ----
    Returns

    con, cur
        connection and cursor for the interaction with the db
    """

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

def connect_database_pd(update_database=False, path_secrets_db='db_access.csv',
                     dbname=None):
    """ Connect to psql database via sqlalchemy for use with pandas """
    
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

    user = info_db['user'][0]
    password = info_db['password'][0]
    port = info_db['port'][0]
    dbname = info_db['dbname'][0]
    url = "postgresql://{0}:{1}@{2}:{3}/{4}".format(user,password,host,port,dbname)
    
    con = create_engine(url)

    return con

def add_cv(con, cur, filename):
    """Check which names are in the CV list"""

    cv_list = ascii.read(filename)['name']

    # Get all names
    cur.execute("select name from candidate")
    r = cur.fetchall()
    names = list(l[0] for l in r)
    names_match = list(n for n in names if n in cv_list)
    print(names_match)
    print("Number of matches:", len(names_match))

    # Repeat for those with hard rejection
    cur.execute("select name from candidate where hard_reject = 0")
    r = cur.fetchall()
    names = list(l[0] for l in r)
    names_match = list(n for n in names if n in cv_list)
    print(names_match)
    print("Number of matches not hard rejected:", len(names_match))


def add_fbot(con, cur, filename):
    """Check which names are in the FBOT candidate list"""

    fbot_list = ascii.read(filename, format='csv')['name']

    # Get all names
    cur.execute("select name from candidate")
    r = cur.fetchall()
    names = list(l[0] for l in r)
    names_match = list(n for n in names if n in fbot_list)
    print(names_match)
    print("Number of matches:", len(names_match))

    # Repeat for those with hard rejection
    cur.execute("select name from candidate where hard_reject = 0")
    r = cur.fetchall()
    names = list(l[0] for l in r)
    names_match = list(n for n in names if n in fbot_list)
    print(names_match)
    print("Number of matches not hard rejected:", len(names_match))


def check_glade(ra, dec, rad=1, dist_min=4):
    """Check if the sources have possible match in GLADE"""

    coord = SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg))

    glade = Vizier.query_region(coord, radius=rad*u.deg,
                                   catalog="VII/281/glade2")

    if len(glade) == 0:
        return None, None, None

    # Check if any of the galaxies found are within an
    # acceptable distance
    glade_select = []
    dist_kpc_select = []
    sep_select = []

    for g in glade[0]:
        # Ignore galaxies too nearby
        if g['Dist'] < dist_min:
            continue
        sep = SkyCoord(ra=g['RAJ2000']*u.deg,
                       dec=g['DEJ2000']*u.deg).separation(coord)
        dist_kpc = g['Dist']*(10**3)*np.sin(sep)/np.cos(sep)
        if dist_kpc >= 0 and dist_kpc < 120:
            glade_select.append(g)
            dist_kpc_select.append(dist_kpc)
            sep_select.append(sep)

    if len(glade_select) == 0:
        return None, None, None
    else:
        return glade_select, sep_select, dist_kpc_select


def ls_coords_datalab(ra, dec, radius_nuclear=1.5):
    '''Query legacy survey via DataLab'''

    from dl import queryClient as qc
    result = qc.query(sql='SELECT ra,dec from smash_dr1.object LIMIT 10')
    print(result)
    exit()

    return


def populate_table_ls(tbl, con, cur, catalog='dr8_south'):
    '''Query the legacy survey table and populate it'''

    #Read the secrets file and make the connection to the database
    info = ascii.read('./db_access.csv', format='csv')

    info_photoz = info[info['db'] == 'photoz']
    db_photoz = f"host={info_photoz['host'][0]} port={info_photoz['port'][0]} dbname={info_photoz['dbname'][0]} user={info_photoz['user'][0]} password={info_photoz['password'][0]}"
    connection_photoz = psycopg2.connect(db_photoz)
    cursor_photoz = connection_photoz.cursor()

    # remove those candidates that already have a match
    cur.execute("select name from crossmatch where ls_sep_arcsec is not NULL")
    r = cur.fetchall()
    names_skip = list(l[0] for l in r)
    print(names_skip, len(set(names_skip)))
    names = set(list(n for n in list(tbl['name'])))
    for name in list(names):
        if name in names_skip:
            continue
        ra = np.mean(tbl['ra'][tbl['name'] == name])
        dec = np.mean(tbl['dec'][tbl['name'] == name])
        #####ls_coords_datalab(ra, dec, radius_nuclear=1.5)
        ######
        ####exit()
        ######
        ls_info = query_coords_nuclear(ra, dec, cursor_photoz,
                                       radius_nuclear=1.5,
                                       datalab=False,
                                       check_quality=False,
                                       catalog=catalog)
        if ls_info is None:
            continue
        marks = '?,?,?,?,?,?,?,?,?,?,?,?,?'
        if use_sqlite is False:
            marks.replace("?", "%s")
        for l in ls_info:
            cur.execute(f"INSERT INTO crossmatch (id, name, \
                        ls_ra, ls_dec, ls_sep_arcsec, ls_z_spec, \
                        ls_z_phot_median, ls_z_phot_std, ls_type, \
                        ls_z_phot_l95, ls_z_phot_u95, ls_fluxz, \
                        ls_photoz_checked) \
                        VALUES ({marks})",
                        (None, name, l['ls_ra'],
                         l['ls_dec'], l['ls_sep_arcsec'], l['ls_z_spec'],
                         l['ls_z_phot_median'], l['ls_z_phot_std'],
                         l['ls_type'], l['ls_z_phot_l95'],
                         l['ls_z_phot_u95'], l['ls_fluxz'],
                         l['ls_photoz_checked']))

        con.commit()


def get_sigma_limit(cursor_photoz, ra_field, dec_field, flux, dmag=0.1, query_radius=2, confidence=0.95, catalog='dr8_north'):
    """For a given position in the sky, create a distribution of z_phot_std
    for a certain mag_z"""

    dflux = abs(flux * 10**(-dmag/2.5) - flux)
    flux_max = flux + dflux
    flux_min = flux - dflux
    query = f'''SELECT z_phot_median, z_phot_std, z_spec, "PARALLAX", "FLUX_Z", q3c_dist({ra_field}, {dec_field}, "RA", "DEC") * 3600 as sep_arcsec FROM {catalog} where q3c_radial_query("RA", "DEC", {ra_field}, {dec_field},{query_radius}) and "FLUX_Z" > {flux_min} and "FLUX_Z" < {flux_max} and "PARALLAX" = 0 and z_phot_median > 0 ORDER BY sep_arcsec '''
    cursor_photoz.execute(query)
    result = cursor_photoz.fetchall()
    z_phot_std_list = list(r[1] for r in result)
    std = np.std(z_phot_std_list)
    m = np.median(z_phot_std_list)
    std_median = np.sqrt(np.median(abs(z_phot_std_list - m)**2))

    return m,std_median


def query_coords_nuclear(ra,dec,cursor_photoz,radius=0.1, radius_photoz='0.5',
                         radius_nuclear=1., catalog='dr8_north', datalab=True,
                         check_quality=True):
    '''Query the database to search for matches at the given coordinates'''
    print("Starting query")

    #Crossmatch with photoz database
    if datalab is True:
        from dl import queryClient as qc
        '''
        query = "SELECT pz.z_phot_median, pz.z_phot_std, pz.z_phot_l95, pz.z_phot_u95, \
                pz.z_spec, "PARALLAX", "FLUX_Z", \
                q3c_dist({ra}, {dec}, "RA", "DEC") * 3600  as sep_arcsec, \
                "RA", "DEC", "TYPE" from ls_dr8.photoz LIMIT 10"
        result = qc.query(sql=query)
        '''
    else:
        query = f'''SELECT z_phot_median, z_phot_std, z_phot_l95, z_phot_u95, \
                z_spec, "PARALLAX", "FLUX_Z", \
                q3c_dist({ra}, {dec}, "RA", "DEC") * 3600  as sep_arcsec, \
                "RA", "DEC", "TYPE" \
                FROM {catalog} where \
                q3c_radial_query("RA", "DEC", {ra}, {dec}, 20 * 0.0002777) \
                ORDER BY sep_arcsec'''

        cursor_photoz.execute(query)
        result = cursor_photoz.fetchall()

    if len(result) == 0:
        print("No match")
        return None
    else:
        print("Found a match")
    nuclear = False
    check_int = None
    info_out = []
    if result[0][7] <= radius_nuclear:
        print("NUCLEAR")
        nuclear = True
    else:
        print("OFF-NUCLEAR")
    checked_z_std = False
    for r in result:
        #If nuclear, skip everything further than 5 arcsec
        if nuclear and r[7] > 5.:
            continue
        #Select only good mags
        flux = r[6]
        mag_z = -2.5 * np.log10(flux)+22.5
        photoz, photoz_err = None, None

        if mag_z < 21. and check_quality is True:
            #Get information about the distribution of errors 
            #in the std of the photoz 
            if checked_z_std is False:
                median_std_magz, std_magz = get_sigma_limit(cursor_photoz, ra,
                                                            dec, flux, dmag=0.1,
                                                            query_radius=2,
                                                            confidence=0.95)
                checked_z_std = True

            if r[1] < median_std_magz + 2*std_magz:
                photoz, photoz_err = r[0], r[1]
                check_int = 1
        elif mag_z < 21. and check_quality is False:
            check_int = 0
            photoz, photoz_err = r[0], r[1]
        info = {'ls_ra': r[8], 'ls_dec': r[9], 'ls_sep_arcsec': r[7],
                'ls_z_spec': r[4], 'ls_z_phot_median': photoz,
                'ls_z_phot_std': photoz_err, 'ls_type': r[10],
                'ls_z_phot_l95': r[2], 'ls_z_phot_u95': r[3],
                'ls_fluxz': r[6],'ls_photoz_checked': check_int}
        info_out.append(info)

    return info_out


def set_dist_kpc_ls(con, cur):
    """Calculate the separation in kpc for LS photoz"""

    query = "select id, ls_z_spec, ls_z_phot_median, ls_z_phot_l95, \
             ls_z_phot_u95, ls_sep_arcsec from crossmatch where \
             (ls_z_spec > 0) \
             or (ls_z_phot_median > 0)"
    r = pd.read_sql_query(query, con)
    info = Table.from_pandas(r)
    for l in info:
        if l['ls_z_spec'] > 0 and np.ma.is_masked(l['ls_z_spec']) is False and np.isnan(l['ls_z_spec']) == False:
            sep_kpc_spec = np.tan(l['ls_sep_arcsec']*u.arcsec)*cosmo.luminosity_distance(l['ls_z_spec']).to(u.kpc).value
            cur.execute(f"UPDATE crossmatch SET \
                        ls_dist_kpc_spec = {sep_kpc_spec} \
                        where id = {l['id']}")
        else:
            sep_kpc_spec = np.nan
        if l['ls_z_phot_median'] > 0 and np.ma.is_masked(l['ls_z_phot_median']) == False and np.isnan(l['ls_z_phot_median']) == False:
            sep_kpc_phot_median = np.tan(l['ls_sep_arcsec']*u.arcsec)*cosmo.luminosity_distance(l['ls_z_phot_median']).to(u.kpc).value
            cur.execute(f"UPDATE crossmatch SET \
                        ls_dist_kpc_phot_median = {sep_kpc_phot_median} \
                        where id = {l['id']}")
        else:
            sep_kpc_phot_median = np.nan
        if l['ls_z_phot_l95'] > 0 and np.ma.is_masked(l['ls_z_phot_l95']) == False and np.isnan(l['ls_z_phot_l95']) == False:
            sep_kpc_phot_l95 = np.tan(l['ls_sep_arcsec']*u.arcsec)*cosmo.luminosity_distance(l['ls_z_phot_l95']).to(u.kpc).value
            cur.execute(f"UPDATE crossmatch SET \
                        ls_dist_kpc_phot_u95 = {sep_kpc_phot_l95} \
                        where id = {l['id']}")
        else:
            sep_kpc_phot_l95 = np.nan
        if l['ls_z_phot_u95'] > 0 and np.ma.is_masked(l['ls_z_phot_u95']) == False and np.isnan(l['ls_z_phot_u95']) == False:
            sep_kpc_phot_u95 = np.tan(l['ls_sep_arcsec']*u.arcsec)*cosmo.luminosity_distance(l['ls_z_phot_u95']).to(u.kpc).value
            cur.execute(f"UPDATE crossmatch SET \
                        ls_dist_kpc_phot_l95 = {sep_kpc_phot_u95} \
                        where id = {l['id']}")
        else:
            sep_kpc_phot_u95 = np.nan
    # Commit the changes
    con.commit()


def populate_table_glade(tbl, con, cur):
    '''Retrieve GLADE information and populate
    the crossmatch table'''

    names = set(list(n for n in list(tbl['name'])))

    # remove those candidates that already have a match
    cur.execute("select name from candidate where glade_match is not NULL")
    r = cur.fetchall()
    names_skip = list(l[0] for l in r)

    # Marks for the ingestion
    marks = '?,?,?,?,?,?,?,?,?,?,?,?,?'
    if use_sqlite is False:
        marks = marks.replace("?", "%s")
        cur.execute("SELECT MAX(id) from crossmatch")
        maxid = cur.fetchall()[0][0]
    else:
        maxid = None
    for name in list(names):
        if name in names_skip:
            continue
        ra = np.mean(tbl['ra'][tbl['name'] == name])
        dec = np.mean(tbl['dec'][tbl['name'] == name])
        glade_match = 0
        try:
            glade_list, sep_list, dist_kpc_list = check_glade(ra, dec)
        except timeout:
            print("Time out on {name}?")
            print("Trying again")
            try:
                glade_list, sep_list, dist_kpc_list = check_glade(ra, dec)
            except timeout:
                print("Problem persisting for {name}")
                continue

        if glade_list is not None:
            glade_match = 1
            for g, s, dk in zip(glade_list, sep_list, dist_kpc_list):
                for k in g.colnames:
                    if np.ma.is_masked(g[k]):
                        try:
                            g[k] = None
                        except TypeError:
                            if k == 'PGC':
                                g[k] = 1
                if use_sqlite is False:
                    maxid += 1
                cur.execute(f"INSERT INTO crossmatch (id, name, \
                            glade_sep_arcsec, glade_name_gwgc, \
                            glade_name_sdss, glade_name_hl,\
                            glade_name_2mass, glade_ra, glade_dec, \
                            glade_z, glade_dist_mpc, \
                            glade_dist_kpc, glade_bmag) \
                            VALUES ({marks})",
                            (maxid, name, float(s.arcsec), g['GWGC'],
                             g['SDSS-DR12'], g['HyperLEDA'],
                             g['_2MASS'], g['RAJ2000'], g['DEJ2000'],
                             g['z'], g['Dist'], float(dk),
                             float(g['Bmag'])))
                print(f"Match found for {name}, sep {s.arcsec}, dist_kpc {dk}")
        # Add info to the candidate table
        cur.execute(f"UPDATE candidate SET \
                    glade_match = {glade_match} \
                    where name = '{name}'")

        con.commit()


def check_gaia(ra, dec, rad=1.5):
    """Check if the sources have a counterpart in Gaia and return its
    parallax."""

    coord = SkyCoord(ra=ra, dec=dec,unit=(u.deg, u.deg))
    gaia = Vizier.query_region(coord, radius=rad*u.arcsec,
                               catalog="I/345/gaia2")
    if len(gaia) == 0:
        return None, None, None

    else:
        plx = gaia[0]['Plx'][gaia[0]['Plx'] == np.max(gaia[0]['Plx'])][0]
        if np.ma.is_masked(plx):
            plx = None
        sep = coord.separation(SkyCoord(ra=gaia[0]['RA_ICRS'],
                                        dec=gaia[0]['DE_ICRS']))

    # Note that plx is a scalar!
    return gaia[0], sep, plx


def populate_table_gaia(tbl, con, cur):
    '''Retrieve Gaia information and populate
    the crossmatch table'''
    
    names = set(list(n for n in list(tbl['name'])))
    # remove those candidates that already have a match
    cur.execute("select name from candidate where mindet3 is not NULL")
    r = cur.fetchall()
    names_skip = list(l[0] for l in r)

    # Marks for the ingestion
    marks = '?,?,?,?,?,?,?,?'
    if use_sqlite is False:
        marks = marks.replace("?","%s")
        cur.execute("SELECT MAX(id) from crossmatch")
        maxid = cur.fetchall()[0][0]
    else:
        maxid = None

    for name in list(names):
        if name in names_skip:
            continue
        ra = np.mean(tbl['ra'][tbl['name'] == name])
        dec = np.mean(tbl['dec'][tbl['name'] == name])
        gaia, sep, plx = check_gaia(ra, dec)
        gaia_match = 0
        gaia_stellar = 0
        if gaia is None:
            continue
        else:
            gaia_match = 1
            for g, s in zip(gaia, sep):
                if use_sqlite is False:
                    maxid += 1
                for k in g.colnames:
                    if np.ma.is_masked(g[k]):
                        g[k] = None
                cur.execute(f"INSERT INTO crossmatch (id, name, gaia_plx, \
                            gaia_bmag, gaia_bmagerr, gaia_ra, gaia_dec, \
                            gaia_sep_arcsec) VALUES ({marks})",
                            (maxid, name,  g['Plx'], g['BPmag'], g['e_BPmag'],
                             g['RA_ICRS'], g['DE_ICRS'], s.arcsec))

            # Check if the plx is such that the source can be
            # classified stella. For the threshold see Luri+2018
            if plx is not None and np.abs(plx) > 1.081:
                gaia_stellar = 1
                print(f"Found stellar source {name}")
            # Add info to the candidate table
            cur.execute(f"UPDATE candidate SET \
                         gaia_match = {gaia_match}, \
                         gaia_stellar = {gaia_stellar} \
                         where name = '{name}'")
        con.commit()


def populate_table_lightcurve(tbl, con, cur, programids=[1,2,3]):
    '''Add the lightcurve information for each candidate'''

    # remove those candidates that already have an entry
    str_names = "'"+"','".join(list(tbl['name']))+"'"
    cur.execute(f"select name, jd from lightcurve \
where name in ({str_names})")
    r = cur.fetchall()
    names_skip = list((l[0], l[1]) for l in r)

    # Marks for the ingestion
    #cur.execute("SELECT * FROM lightcurve LIMIT 0;")
    #num_cols = len([desc[0] for desc in cur.description])
    num_cols = 15
    marks = ",".join(["%s"]*num_cols)
    cur.execute("SELECT MAX(id) from lightcurve")
    maxid = cur.fetchall()[0][0]

    if maxid is None:
        maxid = 0

    for l in tbl:
        # Skip if the combination name+jd is already present
        if (l['name'], l['jd']) in names_skip:
            continue
        if not int(l['programid']) in programids:
            continue
        maxid += 1
        cur.execute(f"INSERT INTO lightcurve (id, name, ra, dec, \
                    jd, magpsf, sigmapsf, filter, \
                    magzpsci, magzpsciunc, programid, \
                    field, rcid, pid, origin) VALUES ({marks})",
                    (maxid, l['name'],l['ra'], l['dec'], l['jd'],
                     np.float64(l['magpsf']), np.float64(l['sigmapsf']),
                     l['filter'], np.float64(l['magzpsci']),
                     np.float64(l['magzpsciunc']), int(l['programid']),
                     int(l['field']), int(l['rcid']), int(l['pid']), l['origin']
                     ))
    con.commit()


def populate_table_candidate(tbl, con, cur):
    '''Populate the table with candidate information'''

    names = set(list(n for n in list(tbl['name'])))

    # remove those candidates that already have an entry
    str_names = "'"+"','".join(list(names))+"'"
    cur.execute(f"select name from candidate \
where name in ({str_names})")
    r = cur.fetchall()
    names_skip = list(l[0] for l in r)

    # Marks for the ingestion
    #cur.execute("SELECT * FROM candidate LIMIT 0;")
    #num_cols = len([desc[0] for desc in cur.description])
    marks = ",".join(["%s"]*9)

    for name in list(names):
        if name in names_skip:
            continue
        ra = np.mean(tbl['ra'][tbl['name'] == name])
        dec = np.mean(tbl['dec'][tbl['name'] == name])
        sgscore1 = list(s['sgscore1'] for s in tbl[tbl['name']==name] if s['sgscore1']>=0 and s['sgscore1'] <= 1.)
        sgscore2 = list(s['sgscore2'] for s in tbl[tbl['name']==name] if s['sgscore2']>=0 and s['sgscore2'] <= 1.)
        sgscore3 = list(s['sgscore3'] for s in tbl[tbl['name']==name] if s['sgscore3']>=0 and s['sgscore3'] <= 1.)
        distpsnr1 = list(s['distpsnr1'] for s in tbl[tbl['name']==name] if s['distpsnr1']>=0)
        distpsnr2 = list(s['distpsnr2'] for s in tbl[tbl['name']==name] if s['distpsnr2']>=0)
        distpsnr3 = list(s['distpsnr3'] for s in tbl[tbl['name']==name] if s['distpsnr3']>=0)

        if len(sgscore1) >= 1:
            sgscore1 = sgscore1[0]
            distpsnr1 = distpsnr1[0]
        else:
            sgscore1 = None
            distpsnr1 = None
        if len(sgscore2) >= 1:
            sgscore2 = sgscore2[0]
            distpsnr2 = distpsnr2[0]
        else:
            sgscore2 = None
            distpsnr2 = None
        if len(sgscore3) >= 1:
            sgscore3 = sgscore3[0]
            distpsnr3 = distpsnr3[0]
        else:
            sgscore3 = None
            distpsnr3 = None
        cur.execute(f"INSERT INTO candidate (name, ra, dec, \
                    sgscore1, sgscore2, sgscore3, \
                    distpsnr1, distpsnr2, distpsnr3) \
                    VALUES ({marks})",
                    (name, np.float64(ra), np.float64(dec),
                    np.float64(sgscore1), np.float64(sgscore2),
                    np.float64(sgscore3),
                    np.float64(distpsnr1), np.float64(distpsnr2),
                    np.float64(distpsnr3)))

    con.commit()


def populate_table_clu(con, cur, tbl=None,
                       max_dist=100.,
                       path_clu='CLU_20190708_marshalFormat.hdf5'):
    '''
    Crossmatch the candidates with the CLU galaxy catalog.

    ---
    Parameters

    con, cur
        connection and cursor for the psql database

    tbl astropy.table
        photometry table, if provided all the candidates in the table
        will be crossmatched. If None, all the candidates in the db
        with clu_match NULL will be crossmatched.

    max_sep float
        largest projected distance (in kpc) from CLU galaxies

    path_clu str
        path to the CLU catalog filename

    ---
    Returns

    It populates the crossmatch table with CLU galaxies
    and it updates clu_match with a boolean value in the
    candidate table.
    '''
    import numpy as np
    if tbl is None:
        # Get candidates that do not have been matched already
        cur.execute("select name, ra, dec from candidate \
where clu_match is NULL")
        r = cur.fetchall()
        names = list(l[0] for l in r)
        ra = list(l[1] for l in r)
        dec = list(l[2] for l in r)
    else:
        names = list(n for n in set(tbl['name']))
        ra = list(np.mean(tbl['ra'][tbl['name'] == name])
                  for name in list(names))
        dec = list(np.mean(tbl['dec'][tbl['name'] == name])
                   for name in list(names))
    coords = SkyCoord(ra=np.array(ra)*u.deg, dec=np.array(dec)*u.deg)

    # Marks for the ingestion
    #cur.execute("SELECT * FROM crossmatch LIMIT 0;")
    #num_cols = len([desc[0] for desc in cur.description])
    marks = ",".join(["%s"]*24)
    cur.execute("SELECT MAX(id) from crossmatch")
    maxid = cur.fetchall()[0][0]
    if maxid is None:
        maxid = 0

    # Read CLU
    clu = read_table_hdf5(path_clu)
    clu = clu[clu['distmpc'] > 4]
    clu_coords = SkyCoord(ra=clu['ra']*u.deg, dec=clu['dec']*u.deg)

    names_match = []
    names_no_match = []
    for name, coord in zip(names, coords):
        sep = clu_coords.separation(coord)
        dist_kpc = clu['distmpc']*(10**3)*np.sin(sep)/np.cos(sep)
        condition0 = dist_kpc >= 0
        clu_match = clu[condition0]
        sep_match = sep[condition0]
        dist_kpc = dist_kpc[condition0]
        condition = dist_kpc < 120
        clu_match = clu_match[condition]
        sep_match = sep_match[condition]
        dist_kpc = dist_kpc[condition]

        if len(clu_match) > 0:
            names_match.append(name)
            for c, d, s in zip(clu_match, dist_kpc, sep_match):
                maxid += 1
                cur.execute(f"INSERT INTO crossmatch (id, name, clu_id, \
                clu_ra, clu_dec, clu_z, clu_zerr, clu_distmpc, \
                clu_mstar, clu_sfr_fuv, clu_sfr_ha, \
                clu_w1mpro, clu_w1sigmpro, clu_w2mpro, \
                clu_w2sigmpro, clu_w3mpro, clu_w3sigmpro, \
                clu_w4mpro, clu_w4sigmpro, clu_type_ned, \
                clu_a, clu_b2a, clu_dist_kpc, clu_sep_arcsec) \
                VALUES ({marks})",
                (maxid, name, int(c['cluid']),
                 c['ra'], c['dec'], c['z'], c['zerr'],
                 c['distmpc'], c['mstar'], c['sfr_fuv'], c['sfr_ha'],
                 c['w1mpro'], c['w1sigmpro'],
                 c['w2mpro'], c['w2sigmpro'],
                 c['w3mpro'], c['w3sigmpro'],
                 c['w4mpro'], c['w4sigmpro'],
                 c['type_ned'], c['a'], c['b2a'], float(d), s.arcsec))
            con.commit()
        else:
            names_no_match.append(name)

        # Update the candidate table
        if len(names_match) > 0:
            names_match_str = "'" + "','".join(names_match) + "'"
            cur.execute(f"UPDATE candidate SET \
                        clu_match = 1 \
                        where name in ({names_match_str})")
        if len(names_no_match) > 0:
            names_no_match_str = "'" + "','".join(names_no_match) + "'"
            cur.execute(f"UPDATE candidate SET \
                        clu_match = 0 \
                        where name in ({names_no_match_str})")
    # Commit the changes
    con.commit()


def add_columns(con, cur):
    """Add additional necessary columns"""

    print("Adding new columns")
    cur.execute("ALTER TABLE candidate ADD index_fade_stack_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_stack_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_fade_stack_i FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_stack_g FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_stack_r FLOAT")
    cur.execute("ALTER TABLE candidate ADD index_rise_stack_i FLOAT")
    print("Added 0")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_spec FLOAT")
    print("Added 1")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_phot_median FLOAT")
    print("Added 2")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_phot_u95 FLOAT")
    print("Added 4")
    cur.execute("ALTER TABLE crossmatch ADD ls_dist_kpc_phot_l95 FLOAT")
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
    cur.execute("ALTER TABLE candidate ADD hard_reject INT")

    # first detection in a given band
    cur.execute("ALTER TABLE candidate ADD max_days_g FLOAT")
    print("Added all.")
    con.commit()


def set_hard_reject(con, cur, reject_filename, append=True):
    """Read the file with hard-rejected candidates and update the db"""

    t = ascii.read(reject_filename, format='csv')
    names_reject = set(list(n for n in list(t['name'])))

    # Get all the names of the candidates in the db
    if append is True:
        cur.execute("select name from candidate where hard_reject is NULL")
    else:
        cur.execute("select name from candidate")
    r = cur.fetchall()
    names_db = list(l[0] for l in r) 

    for name in names_db:
        if name in names_reject:
            cur.execute(f"UPDATE candidate SET \
                         hard_reject = 1 where name = '{name}'")
        else:
            cur.execute(f"UPDATE candidate SET \
                         hard_reject = 0 where name = '{name}'")
    con.commit()

    # Check:
    cur.execute("select name from candidate WHERE hard_reject = 1")
    r = cur.fetchall()
    names_check = list(l[0] for l in r)
    print(f"There are {len(names_check)}/{len(names_db)} candidates hard rejected")
    print(f"{len(names_check)}/{len(t)} names in the hard_reject list.")


def create_table_candidate(con, cur):
    # Table for the candidates
    cur.execute("CREATE TABLE candidate(name TEXT NOT NULL PRIMARY KEY, \
                ra DOUBLE, dec DOUBLE, \
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
    # commit the changes
    con.commit()


def create_table_crossmatch(con, cur):
    cur.execute("CREATE TABLE crossmatch(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, clu_id INT, clu_ra DOUBLE, clu_dec DOUBLE, \
                clu_z FLOAT, clu_zerr FLOAT, clu_distmpc FLOAT, \
                clu_mstar FLOAT, clu_sfr_fuv FLOAT, clu_sfr_ha FLOAT, \
                clu_w1mpro FLOAT, clu_w1sigmpro FLOAT, clu_w2mpro FLOAT, \
                clu_w2sigmpro FLOAT, clu_w3mpro FLOAT, clu_w3sigmpro FLOAT, \
                clu_w4mpro FLOAT, clu_w4sigmpro FLOAT, clu_type_ned BLOB, \
                clu_a FLOAT, clu_b2a FLOAT, clu_dist_kpc FLOAT, clu_sep_arcsec FLOAT, \
                glade_sep_arcsec FLOAT, glade_name_gwgc TEXT, \
                glade_name_sdss TEXT, glade_name_hl TEXT,\
                glade_name_2mass TEXT, glade_ra DOUBLE, glade_dec DOUBLE, \
                glade_z FLOAT, glade_dist_mpc FLOAT, glade_dist_mpc_err FLOAT,\
                glade_dist_kpc FLOAT, glade_bmag FLOAT, glade_bmagerr FLOAT, \
                gaia_plx FLOAT, gaia_bmag FLOAT, gaia_bmagerr FLOAT, \
                gaia_ra DOUBLE, gaia_dec DOUBLE, gaia_sep_arcsec FLOAT,\
                ls_ra DOUBLE, ls_dec DOUBLE, ls_sep_arcsec FLOAT, ls_z_spec FLOAT, \
                ls_z_phot_median FLOAT, ls_z_phot_std FLOAT, \
                ls_photoz_checked INT,  ls_type TEXT, \
                ls_z_phot_l95 FLOAT, ls_z_phot_u95 FLOAT, ls_fluxz FLOAT)")
    # commit the changes
    con.commit()


def create_table_lightcurve(con, cur):
    cur.execute("CREATE TABLE lightcurve(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, ra DOUBLE, dec DOUBLE, jd FLOAT, \
                magpsf FLOAT, sigmapsf FLOAT, filter TEXT, \
                magzpsci FLOAT, magzpsciunc FLOAT, programid INT, \
                field INT, rcid BIGINT, pid INT)")
    # commit the changes
    con.commit()


def create_table_lightcurve_forced(con, cur):
    cur.execute("CREATE TABLE lightcurve_forced(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, ra DOUBLE, dec DOUBLE, jd FLOAT, \
                magpsf FLOAT, sigmapsf FLOAT, filter TEXT, \
                magzpsci FLOAT, magzpsciunc FLOAT, programid INT, \
                field INT, rcid BIGINT, pid INT)")
    # commit the changes
    con.commit()

 
def create_table_lightcurve_stacked(con, cur):
    cur.execute("CREATE TABLE lightcurve_stacked(id INTEGER NOT NULL PRIMARY KEY, \
                name TEXT, jd FLOAT, \
                flux FLOAT, flux_unc FLOAT, \
                mag FLOAT, mag_unc FLOAT, limmag FLOAT, filter TEXT, \
                zp FLOAT, ezp FLOAT, programid INT, \
                field INT, ccdid INT, qid INT)")
    # commit the changes
    con.commit()

def populate_table_lightcurve_alertfp(tbl, con, cur, programids=[1,2,3]):
    """
    Populate table with forced photometry from ZTF alerts
    """
     
    str_names = "'"+"','".join(list(tbl['name']))+"'"
    cur.execute(f"select name, jd from lightcurve_alertfp \
where name in ({str_names})")
    r = cur.fetchall()
    names_skip = list((l[0], l[1]) for l in r) # list of (name,jd) already present in db
    #print(names_skip)
    # Marks for the ingestion
    cur.execute("SELECT * FROM lightcurve_alertfp LIMIT 0;")
    num_cols = len([desc[0] for desc in cur.description])
    marks = ",".join(["%s"]*num_cols) # 19 cols
    cur.execute("SELECT MAX(id) from lightcurve_alertfp")
    maxid = cur.fetchall()[0][0]
    if maxid is None:
        maxid = 0
    
    #tbl_to_write = tbl[tbl['origin'] == 'alertfp']
    tbl_to_write = tbl[:]
    
    for l in tbl_to_write:
        # Skip if the combination name+jd is already present
        if (l['name'], l['jd']) in names_skip:
            continue
        if not int(l['programid']) in programids:
            continue
        maxid += 1
        cur.execute(f"INSERT INTO lightcurve_alertfp (id, name, ra, dec, \
                    jd, magpsf, sigmapsf, limmag, snr, filter, \
                    magzpsci, magzpsciunc, programid, \
                    field, rcid, pid, forcediffimflux, forcediffimfluxunc, origin) VALUES ({marks})",
                    (maxid, l['name'], l['ra'], l['dec'], l['jd'],
                     np.float64(l['magpsf']), np.float64(l['sigmapsf']), 
                     np.float64(l['limmag']), np.float64(l['snr']),
                     l['filter'], np.float64(l['magzpsci']),
                     np.float64(l['magzpsciunc']), int(l['programid']),
                     int(l['field']), int(l['rcid']), int(l['pid']),
                     np.float64(l['forcediffimflux']), np.float64(l['forcediffimfluxunc']), 
                     l['origin']
                     ))
    con.commit()
    

def populate_table_lightcurve_forced(con, cur, tbl, targetdir_base,
                                     programids=[1,2,3]):
    """Populate the table with forced photometry measurements
    from the maxlike output FITS files."""

    # Combinations of name and jd already present, to skip
    cur.execute("select name, jd from lightcurve_forced")
    r = cur.fetchall()
    names_skip = list((l[0], l[1]) for l in r)

    # Get max ID
    cur.execute("SELECT MAX(id) from lightcurve_forced")
    maxid = cur.fetchall()[0][0]
    if maxid is None:
        maxid = 0

    for name in set(tbl['name']):
        # File names like force_phot_ZTF17aacbbmv_maxlikelihood_lc.fits
        filename = glob.glob(f"{targetdir_base}/{name}/l*/for*maxl*fits")
        if filename == []:
            continue
        forced = Table(fits.open(filename[0])[1].data)
        forced.rename_column('jdobs', 'jd')
        forced.rename_column('fieldid', 'field')
        # Upload the results in the database
        for l in forced:
            if (name, l['jd']) in names_skip:
                continue
            if not int(l['programid']) in programids:
                # FIXME remove, only as a check for tomorrow
                print("SKIPPED DATA POINT")
                continue
            keys = l.colnames
            keys_string = ", ".join(keys)
            marks = ", ".join(['%s']*(2+len(keys)))
            maxid += 1
            values = tuple([maxid, name] + [l[k] for k in keys])
 
            cur.execute(f"INSERT INTO lightcurve_forced (id, name, \
                        {keys_string}) \
                        VALUES ({marks})",
                        values)

    con.commit() 


def populate_table_lightcurve_stacked(con, cur, names):
    """Populate the table with stacking of forced photometry measurements
    from the forced photometry light curve in the psql db"""
    from select_variability_db import stack_lc
    
    # Files to skip
    names_str = "'" + "','".join(list(names)) + "'"
    cur.execute(f"select name, jd from lightcurve_stacked \
where name in ({names_str})")
    r = cur.fetchall()
    names_skip = list((l[0], l[1]) for l in r)
    if len(names_skip) > 0:
        list_jds = list(j[1] for j in names_skip)

    # Get the forced phot light curves
    query = f"select * from lightcurve_forced \
where name in ({names_str})"
    r = pd.read_sql_query(query, con)
    forced_all = Table.from_pandas(r)
    forced_all.rename_column('flux_maxlike', 'Flux_maxlike')
    forced_all.rename_column('flux_maxlike_unc', 'Flux_maxlike_unc')

    # Marks for the ingestion
    marks = ",".join(["%s"]*12)

    # Find the max id
    cur.execute("SELECT MAX(id) from lightcurve_stacked")
    maxid = cur.fetchall()[0][0]
    if maxid is None:
        maxid = 0

    for name in names:
        # Forced phot table for the candidate
        try:
            forced = forced_all[forced_all['name'] == name]
        except IndexError:
            continue
        # Stack the data points
        forced_stack = stack_lc(forced, snt_det=4, days_stack=1.)
        # Upload the results in the database
        for l in forced_stack:
            # Skip data point already present in the db
            if len(names_skip) > 0:
                if (name, l['jd']) in names_skip:
                    continue
                # Check necessary for stacks!
                elif np.min(np.abs(list_jds - l['jd'])) < 1./60/60/24.:
                    continue
                     
            # Increment the ID number
            maxid += 1
            # Upload the photometry in the database
            cur.execute(f"INSERT INTO lightcurve_stacked (id, name, \
                        jd, flux, flux_unc, \
                        mag, mag_unc, limmag, filter, \
                        zp, ezp, programid) \
                        VALUES ({marks})",
                        (maxid, name, float(l['jd']), float(l['flux']),
                         float(l['flux_unc']), float(l['mag']),
                         float(l['mag_unc']), float(l['limmag']),
                         l['filter'], float(l['zp']), float(l['ezp']),
                         int(l['programid'])))
        con.commit()


def get_dust_info(coords):
    """
    Get coordinates as a SkyCoord object
    and returns E(B-V) 
    """
    from dustmaps.planck import PlanckQuery
    #from dustmaps.config import config
    #if config['datadir'] == None:
        # create data directory
    #    config['datadir']= data directory name
    try:
        planck = PlanckQuery()
    except FileNotFoundError:
        # Download the dustmap
        import dustmaps.planck
        dustmaps.planck.fetch()
        planck = PlanckQuery()

    ebv = planck(coords)

    return ebv


def populate_extinction(con, cur):
    """
    Find the dust extinction for all those candidates
    still lacking of ebv information
    """

    coords_tbl = pd.read_sql_query(f"SELECT name, ra, dec FROM candidate \
                                   WHERE ebv is NULL", con)
    if coords_tbl.empty:
        return
    names = coords_tbl["name"]
    ras = np.array(coords_tbl["ra"])
    decs = np.array(coords_tbl["dec"])
    coords = SkyCoord(ra=ras*u.deg, dec=decs*u.deg)
    # Get E(B-V)
    ebv = get_dust_info(coords)

    # Update the database
    for n, e in zip(names, ebv):
        cur.execute(f"UPDATE candidate SET \
                    ebv = {e} \
                    where name = '{n}'")
    # Commit the changes
    con.commit()


def populate_gal_lat(con, cur):
    """
    Find the Galactic latitude for all those candidates
    still lacking of b_gal information
    """
    # Find candidates to update
    coords_tbl = pd.read_sql_query(f"SELECT name, ra, dec FROM candidate \
                                   WHERE b_gal is NULL", con)
    if coords_tbl.empty:
        return
    names = coords_tbl["name"]
    ras = np.array(coords_tbl["ra"])
    decs = np.array(coords_tbl["dec"])
    coords = SkyCoord(ra=ras*u.deg, dec=decs*u.deg)
    # Convert the coordinates from equatorial to Galactic
    bgal = coords.galactic.b

    # Update the database
    for n, b in zip(names, bgal):
        cur.execute(f"UPDATE candidate SET \
                    b_gal = {b.deg} \
                    where name = '{n}'")
    # Commit the changes
    con.commit()


def get_table_names(con,cur):
    cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema NOT IN \
    ('pg_catalog','information_schema') AND table_type='BASE TABLE';")
    tables = [table[0] for table in cur.fetchall()]
    return tables
def get_col_names(tblname, con, cur):
    cur.execute(f"SELECT * FROM {tblname} LIMIT 0")
    colnames = [desc[0] for desc in cur.description]
    return colnames


if __name__ == '__main__':

    # Define the input filename
    tbl = ascii.read('test_lightcurves.csv')

    # Connect to psql db
    # Read the secrets
    info = ascii.read('./db_access.csv', format='csv')
    info_db = info[info['db'] == 'db_kn_rt_admin']
    db_kn = f"host={info_db['host'][0]} dbname={info_db['dbname'][0]} \
port={info_db['port'][0]} user={info_db['user'][0]} \
password={info_db['password'][0]}"

    con, cur = connect_database(update_database=False,
                                path_secrets_db='db_access.csv')
 
    # Create the tables
    #create_table_candidate(con, cur)
    #create_table_crossmatch(con, cur)
    #create_table_lightcurve(con, cur)

    #cur.execute("UPDATE crossmatch SET \
    #             ls_photoz_checked = 1 \
    #             where ls_z_phot_median is not NULL")
    #cur.execute("select name, ls_z_phot_median from crossmatch where ls_photoz_checked is not NULL")
    #a = cur.fetchall()
    #for l in a:
    #    print(l)
    #pdb.set_trace()
    #con.commit()
    #exit()

    #The following commands are to retrieve information from the database
    #for tab in ['candidate', 'crossmatch', 'lightcurve', 'lightcurve_forced', 'lightcurve_stacked']:
    #    cur.execute(f"PRAGMA table_info({tab})")
    #    r = cur.fetchall()
    #    for l in r:
    #        print(l)

    #exit()
    # Populate the tables
    # Write in the db if a clu match is present
    #update_clu_match(con, cur)
    #exit()

    #print("Populating table candidate..")
    #populate_table_candidate(tbl, con, cur)
    #print("populating table lightcurve...")
    #populate_table_lightcurve(tbl, con, cur)
    #print("populating table CLU...")
    #populate_table_clu(tbl, con, cur)
    #print("CLU populated!")
    #print('Populating GLADE..')
    #populate_table_glade(tbl, con, cur)
    #print('GLADE populated!')
    #print('Populating gaia..')
    #populate_table_gaia(tbl, con, cur)
    #print('gaia populated!')
    #print('Populating legacy survey')
    #populate_table_ls(tbl, con, cur, catalog='dr8_south')
    #print('legacy survey populated!')

    # Set boolean values for hard rejection
    #add_columns(con, cur)
    #set_hard_reject(con, cur, '../hard_rejects_candidates2.csv', append=True)
    # Create table for forced photometry
    #create_table_lightcurve_stacked(con, cur)

    # Set the distance in kpc for ls photoz
    #set_dist_kpc_ls(con, cur)
    #print("Populating the forced photometry table")
    #populate_table_lightcurve_forced(con, cur, tbl, clobber_all=False)
    #print("Populating table for stacked forced photometry")
    #populate_table_lightcurve_stacked(con, cur, tbl, clobber_all=False)

    # Match with CV list
    #cv_filename = '../CV_list_23Mar2020.txt'
    #add_cv(con, cur, cv_filename)

    # Match with FBOT cand list
    #fbot_filename = '../list_fbot_anna.csv'
    #add_fbot(con, cur, fbot_filename)

    cur.close()
    con.close()

