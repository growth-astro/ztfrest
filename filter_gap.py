'''
Query Kowalski searching for transients
given a set of constraints.
'''

import json
import requests
import datetime
import pdb
import os

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time, TimeDelta

from penquins import Kowalski
from functions_db import connect_database
import psycopg2
from astropy.io.misc.hdf5 import read_table_hdf5

def xmatch_clu(tbl=None,
               max_dist=200.,
               min_dist=20.,
               max_dist_kpc=130,
               min_dist_kpc=10,
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

    if tbl is None:
        return
    else:
        names = list(n for n in set(tbl['name']))
        ra = list(np.mean(tbl['ra'][tbl['name'] == name])
                  for name in list(names))
        dec = list(np.mean(tbl['dec'][tbl['name'] == name])
                   for name in list(names))
    coords = SkyCoord(ra=np.array(ra)*u.deg, dec=np.array(dec)*u.deg)

    # Read CLU
    clu = read_table_hdf5(path_clu)
    clu = clu[clu['distmpc'] > min_dist]
    clu = clu[clu['distmpc'] < max_dist]
    clu_coords = SkyCoord(ra=clu['ra']*u.deg, dec=clu['dec']*u.deg)

    names_match = []
    names_no_match = []
    for name, coord in zip(names, coords):
        sep = clu_coords.separation(coord)
        dist_kpc = clu['distmpc']*(10**3)*np.sin(sep)/np.cos(sep)
        condition0 = dist_kpc >= min_dist_kpc
        clu_match = clu[condition0]
        sep_match = sep[condition0]
        dist_kpc = dist_kpc[condition0]
        condition = dist_kpc <= max_dist_kpc
        clu_match = clu_match[condition]
        sep_match = sep_match[condition]
        dist_kpc = dist_kpc[condition]

        if len(clu_match) > 0:
            names_match.append(name)
            # Name of the output file
            outfile_clu = f"crossmatch_clu_dist{min_dist}-{max_dist}Mpc_projdist{min_dist_kpc}-{max_dist_kpc}.csv"
            # create the file with the right header if it doesn't exisit
            if os.path.isfile(outfile_clu) is False:
                with open(outfile_clu, "w") as ofclu:
                    ofclu.write("name, clu_id, clu_ra, clu_dec, clu_z, clu_zerr, clu_distmpc, clu_dist_kpc, clu_sep_arcsec\n")
            # Append the information
            with open(outfile_clu, "a") as ofclu:
                for c, d, s in zip(clu_match, dist_kpc, sep_match):
                    ofclu.write(f"{name}, {int(c['cluid'])}, \
{c['ra']}, {c['dec']}, {c['z']}, {c['zerr']}, {c['distmpc']}, \
{d}, {s.arcsec}\n")
        else:
            names_no_match.append(name)

    return names_match, names_no_match


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def print_query_params(args, ra_center, dec_center):
    '''Print a summary of the query parameters'''

    print("#-----")
    print("Cone search parameters:")
    print(f"A list of {len(ra_center)} coordinate pairs will be explored")
    print(f"Search radius {args.radius} arcmin")
    if args.after_trigger or args.jd_trigger > 0:
        print(f"Only sources detected for the first time \
after {Time(args.jd_trigger, format='jd').iso} will be considered")
    print(f"Minimum time between the first and last alert {args.min_days} days")
    print(f"Maximum time between the first and last alert {args.max_days} days")    
    print(f"Query divided in {args.slices} slices")
    print("#-----")
    print(" ")

    return


def get_programidx(program_name, username, password):
    ''' Given a marshal science program name, it returns its programidx'''

    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi',
                      auth=(username, password))
    programs=json.loads(r.text)
    program_dict={p['name']:p['programidx'] for i,p in enumerate(programs)}

    try:
        return program_dict[program_name]
    except KeyError:
        print(f'The user {username} does not have access to \
the program {program_name}')
        return None


def get_candidates_growth_marshal(program_name, username, password):
    ''' Query the GROWTH db for the science programs '''

    programidx=get_programidx(program_name, username, password)
    if programidx==None:
        return None
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', \
        auth=(username, password), data={'programidx':str(programidx)})
    sources=json.loads(r.text)
    sources_out=[]
    for s in sources:
            coords=SkyCoord(ra=s['ra']*u.deg, dec=s['dec']*u.deg, frame='icrs')
            sources_out.append({"name":s['name'],
                                "ra":coords.ra, "dec":coords.dec,
	                        "classification":s['classification'],
                                "redshift":s['redshift'],
                                "creation_date":s['creationdate']})

    return sources_out


def check_clu_transients(sources_kowalski, clu_sources):
    '''Check if the selected sources are present in the 
    CLU science program.  If so, print out the relevant information.'''

    sources_in_clu = []
    sources_not_in_clu = []
    list_clu_sources = list(s['name'] for s in clu_sources)

    for source in sources_kowalski:
        print("-------")
        if source in list_clu_sources:
            clu_source = clu_sources[np.where(np.array(list_clu_sources) == source)[0][0]]
            try:
                for k in clu_source.keys():
                    print(f"{k}: {clu_source[k]}")
                sources_in_clu.append(source)
            except:
                pdb.set_trace()
        else:
            print(f"{source} was not saved in CLU")
            sources_not_in_clu.append(source)
        print("-------")
    print("Summary:")
    print(f"Sources saved in CLU: {sources_in_clu}")
    print(f"Sources not saved in CLU: {sources_not_in_clu}")

    return


def check_lightcurve_alerts(username, password, list_names,
                            min_days, max_days, verbose=True):
    """Re-query light curve info for a list of candidates\
    and check that their full/updated duration is consistent\
    with the time limits provided"""

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts",
                   "filter": {
                              'objectId': {'$in': list(list_names)}
                              },
                   "projection": {
                                  "objectId": 1,
                                  "candidate.jd": 1,
                                  "candidate.ndethist": 1,
                                  "candidate.jdstarthist": 1,
                                  "candidate.jdendhist": 1,
                                  "candidate.jdendhist": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.programid": 1,
                                  }
                       },
            "kwargs": {"hint": "objectId_1"}
             }

    r = k.query(query=q)
    try:
        if r['data'] == []:
            print("No candidates to be checked?")
            return None
    except (KeyError, TypeError) as e:
        if verbose is True:
            print(f"ERROR in getting light curves! attempt {i}" )
            i += 1
        if i > 5:
            print(f"SKIPPING {len(list_names)} light curves, after 5 attempts")
            import pdb
            pdb.set_trace()
            return None

    old = []
    objectid_list = []
    for info in r['data']:
        if info['objectId'] in old:
            continue
        if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) < min_days:
            continue
        if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > max_days:
            old.append(info['objectId'])
        objectid_list.append(info['objectId'])
    clean_set = set(objectid_list)
    #Remove those objects considered old
    for n in set(old):
        try:
            clean_set.remove(n)
        except:
            pass

    return clean_set


def query_kowalski(kow, list_fields, min_days, max_days,
                   ndethist_min, jd, jd_gap=50.,
                   programids=[1,2], verbose=True):
    '''Query kowalski and apply the selection criteria'''

    # Correct the minimum number of detections
    ndethist_min_corrected = int(ndethist_min - 1)

    jd_start = jd
    jd_end = jd + jd_gap

    #Initialize a set for the results
    set_objectId_all = set([])
    for field in list_fields:
        set_objectId_field = set([])
        q = {"query_type": "find",
             "query": {
                       "catalog": "ZTF_alerts",      
                       "filter": {
                                  'candidate.jd': {'$gt': jd_start, '$lt': jd_end},
                                  'candidate.field': int(field),
                                  # FIXME re-introduce rb
                                  ##'candidate.rb': {'$gt': 0.5},
                                  'candidate.drb': {'$gt': 0.65},
                                  'classifications.braai': {'$gt': 0.65},
                                  'candidate.ndethist': {'$gt': ndethist_min_corrected},
                                  'candidate.magpsf': {'$gt': 17},
                                  'candidate.programid': { '$in': programids}
                                  #'candidate.isdiffpos': 't'
                                   },
                       "projection": {
                                      "objectId": 1,
                                      "candidate.rcid": 1,
                                      "candidate.ra": 1,
                                      "candidate.dec": 1,
                                      "candidate.jd": 1,
                                      "candidate.ndethist": 1,
                                      "candidate.jdstarthist": 1,
                                      "candidate.jdendhist": 1,
                                      "candidate.jdendhist": 1,
                                      "candidate.magpsf": 1,
                                      "candidate.sigmapsf": 1,
                                      "candidate.fid": 1,
                                      "candidate.programid": 1,
                                      "candidate.isdiffpos": 1,
                                      "candidate.ndethist": 1,
                                      "candidate.ssdistnr": 1,
                                      "candidate.rb": 1,
                                      "candidate.drb": 1,
                                      "candidate.distpsnr1": 1,   
                                      "candidate.sgscore1": 1,
                                      "candidate.srmag1": 1,
                                      "candidate.distpsnr2": 1,   
                                      "candidate.sgscore2": 1,
                                      "candidate.srmag2": 1,
                                      "candidate.distpsnr3": 1,   
                                      "candidate.sgscore3": 1,
                                      "candidate.srmag3": 1
                                       }
                       },
            "kwargs": {"hint": "jd_field_rb_drb_braai_ndethhist_magpsf_isdiffpos"}
             }

        #Perform the query
        r = kow.query(query=q)
        if verbose is True:
            print(f"Search completed for field {field}, \
{Time(jd, format='jd').iso} + {jd_gap:.1f} days.")


        #Identify 'candid' for all relevant candidates
        objectId_list = []
        with_neg_sub = []
        old = []
        out_of_time_window = []
        stellar_list = []

        # Try to query kowalski up to 5 times
        i = 1
        no_candidates = False
        while i <= 5:
            try:
                if r['data'] == []:
                    no_candidates = True
                break
            except (KeyError, TypeError) as e:
                if verbose is True:
                    print(f"ERROR! jd={jd}, field={field}, attempt {i}" ) 
                i += 1
        if i > 5:
            print(f"SKIPPING jd={jd}, field={field} after 5 attempts")
            continue

        if no_candidates is True:
            if verbose is True:
                print(f"No candidates on jd={jd}, field={field}")
            continue

        for info in r['data']:    
            if info['objectId'] in old:
                continue
            if info['objectId'] in stellar_list:
                continue
            if np.abs(info['candidate']['ssdistnr']) < 10:
                continue
            if info['candidate']['isdiffpos'] in ['f',0]:
                with_neg_sub.append(info['objectId'])
            if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) < min_days:
                continue
            if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > max_days:
                old.append(info['objectId'])
            try:
                if (np.abs(info['candidate']['distpsnr1']) < 1.5 and info['candidate']['sgscore1'] > 0.5):
                    stellar_list.append(info['objectId'])
            except:
                pass
            try:
                if (np.abs(info['candidate']['distpsnr1']) < 15. and
                           info['candidate']['srmag1'] < 15. and
                           info['candidate']['srmag1'] > 0. and
                           info['candidate']['sgscore1'] >= 0.5):
                    continue
            except:
                pass
            try:
                if (np.abs(info['candidate']['distpsnr2']) < 15. and
                           info['candidate']['srmag2'] < 15. and
                           info['candidate']['srmag2'] > 0. and
                           info['candidate']['sgscore2'] >= 0.5):
                    continue
            except:
                pass
            try:
                if (np.abs(info['candidate']['distpsnr3']) < 15. and
                           info['candidate']['srmag3'] < 15. and
                           info['candidate']['srmag3'] > 0. and
                           info['candidate']['sgscore3'] >= 0.5):
                    continue
            except:
                pass

            objectId_list.append(info['objectId'])

        set_objectId = set(objectId_list)

        #Remove those objects with negative subtraction
        for n in set(with_neg_sub):
            try:
                set_objectId.remove(n)
            except:
                pass

        #Remove stellar objects
        for n in set(stellar_list):
            try:
                set_objectId.remove(n)
            except:
                pass

        #Remove those objects considered old
        for n in set(old):
            try:
                set_objectId.remove(n)
            except:
                pass

        #Remove those objects whole alerts go bejond jd_trigger+max_days
        for n in set(out_of_time_window):
            try:
                set_objectId.remove(n)
            except:
                pass
        #print(set_objectId)
        set_objectId_all = set_objectId_all | set_objectId
        #print("Cumulative:", set_objectId_all)

        if verbose is True:
            print("Field", field, len(set_objectId_all))

    return set_objectId_all


def query_and_populate_ls(tbl, con, cur, radius_arcsec=5.,
                         radius_nuclear=1., catalog='both', datalab=True,
                         check_quality=True):
    '''Query the database to search for matches at the given coordinates'''

    # remove those candidates that already have a match
    cur.execute("select name from crossmatch where ls_sep_arcsec is not NULL")
    r = cur.fetchall()
    names_skip = list(l[0] for l in r)
    names = set(list(n for n in list(tbl['name'])))
    names_matched = []
    for name in list(names):
        if name in names_skip:
            continue
        ra = np.mean(tbl['ra'][tbl['name'] == name])
        dec = np.mean(tbl['dec'][tbl['name'] == name])

        # FIXME this is not ncessary with datalab
        if catalog == 'both':
            # First query south
            ls_info = query_coords_ls(ra, dec,
                                       radius_arcsec=radius_arcsec,
                                       radius_nuclear=radius_nuclear,
                                       datalab=datalab,
                                       check_quality=False,
                                       catalog='dr8_south')
            #if ls_info is None:
            #    ls_info = query_coords_ls(ra, dec,
            #                               radius_arcsec=radius_arcsec,
            #                               radius_nuclear=radius_nuclear,
            #                               datalab=datalab,
            #                               check_quality=False,
            #                               catalog='dr8_north')
        else:
            ls_info = query_coords_ls(ra, dec,
                                       radius_arcsec=radius_arcsec,
                                       radius_nuclear=radius_nuclear,
                                       datalab=datalab,
                                       check_quality=False,
                                       catalog=catalog)
        if ls_info is None:
            continue
        # Prepare for inserting into the db
        marks = ",".join(["%s"]*13)
        cur.execute("SELECT MAX(id) from crossmatch")
        maxid = cur.fetchall()[0][0]
        if maxid is None:
            maxid = 0

        for l in ls_info:
            # Remove PSF-shaped underlying sources
            try:
                if l['ls_type'] == 'PSF':
                    continue
            except TypeError:
                continue
            # Remove too nearby measurements
            try:
                if float(l['ls_z_phot_median']) < 0.1:
                    continue
            except TypeError:
                continue
            maxid += 1
            cur.execute(f"INSERT INTO crossmatch (id, name, \
                        ls_ra, ls_dec, ls_sep_arcsec, ls_z_spec, \
                        ls_z_phot_median, ls_z_phot_std, ls_type, \
                        ls_z_phot_l95, ls_z_phot_u95, ls_fluxz, \
                        ls_photoz_checked) \
                        VALUES ({marks})",
                        (maxid, name, l['ls_ra'],
                         l['ls_dec'], l['ls_sep_arcsec'], l['ls_z_spec'],
                         l['ls_z_phot_median'], l['ls_z_phot_std'],
                         l['ls_type'], l['ls_z_phot_l95'],
                         l['ls_z_phot_u95'], l['ls_fluxz'],
                         l['ls_photoz_checked']))
            names_matched.append(name)

        con.commit()

    return list(set(names_matched))


def query_coords_ls(ra,dec,radius_arcsec=5,
                         radius_nuclear=1., catalog='dr8_north', datalab=True,
                         check_quality=True):
    '''Query the database to search for matches at the given coordinates'''

    #Crossmatch with photoz database
    if datalab is True:
        from dl import queryClient as qc
        radius_deg = radius_arcsec / 3600.
        query = qc.query(sql=f"SELECT z_phot_median, z_phot_std, z_phot_l95, ra, dec, \
                             type, flux_z from ls_dr8.photo_z INNER JOIN ls_dr8.tractor \
                             ON ls_dr8.tractor.ls_id = ls_dr8.photo_z.ls_id \
                             where ra > ({ra-radius_deg}) and \
                             ra < ({ra+radius_deg}) and \
                             dec > ({dec-radius_deg}) and \
                             dec < ({dec+radius_deg})")
        result0 = query.split('\n')
        result0 = [r.split(",") for r in result0][1:-1]

        ras = [float(r[3]) for r in result0]
        decs = [float(r[4]) for r in result0]
        result = []
        if len(ras) > 0:
            # Add separation
            gal_coords = SkyCoord(ra=ras*u.deg, dec=decs*u.deg)
            cand_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            sep = cand_coord.separation(gal_coords)
            for i in np.arange(len(ras)):
                result0[i].append(float(sep[i].arcsec))
                # Check that the separation is less than required
                if float(sep[i].arcsec) < radius_arcsec:
                    result.append(result0[i])

        for r in result:
            print(r)
    else:
        #Read the secrets file and make the connection to the database
        info = ascii.read('./db_access.csv', format='csv')

        info_photoz = info[info['db'] == 'photoz']
        db_photoz = f"host={info_photoz['host'][0]} port={info_photoz['port'][0]} dbname={info_photoz['dbname'][0]} user={info_photoz['user'][0]} password={info_photoz['password'][0]}"
        connection_photoz = psycopg2.connect(db_photoz)
        cursor_photoz = connection_photoz.cursor()

        query = f'''SELECT z_phot_median, z_phot_std, z_phot_l95, z_phot_u95, \
                z_spec, "PARALLAX", "FLUX_Z", \
                q3c_dist({ra}, {dec}, "RA", "DEC") * 3600  as sep_arcsec, \
                "RA", "DEC", "TYPE" \
                FROM {catalog} where \
                q3c_radial_query("RA", "DEC", {ra}, {dec}, \
                {radius_arcsec} * 0.0002777)\
                ORDER BY sep_arcsec'''
        cursor_photoz.execute(query)
        result = cursor_photoz.fetchall()

    if len(result) == 0:
        return None

    nuclear = False
    check_int = None
    info_out = []
    if datalab is True:
        separations = [float(s.arcsec) for s in sep]
        if np.min(separations) <= radius_nuclear:
            nuclear = True
    else:
        if result[0][7] <= radius_nuclear:
            nuclear = True
    checked_z_std = False
    for r in result:
        #If nuclear, skip everything further than 5 arcsec
        if nuclear and r[7] > 2*radius_nuclear:
            continue
        #Select only good mags
        flux = float(r[6])
        mag_z = -2.5 * np.log10(flux)+22.5
        photoz, photoz_err = None, None

        # Hard limit for z-band 
        if mag_z > 21:
            continue

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
        if datalab:
            info = {'ls_ra': r[3], 'ls_dec': r[4], 'ls_sep_arcsec': r[7],
                    'ls_z_spec': None, 'ls_z_phot_median': photoz,
                    'ls_z_phot_std': photoz_err, 'ls_type': r[5],
                    'ls_z_phot_l95': r[2], 'ls_z_phot_u95': None,
                    'ls_fluxz': r[6],'ls_photoz_checked': check_int}
        else:
            info = {'ls_ra': r[8], 'ls_dec': r[9], 'ls_sep_arcsec': r[7],
                    'ls_z_spec': r[4], 'ls_z_phot_median': photoz,
                    'ls_z_phot_std': photoz_err, 'ls_type': r[10],
                    'ls_z_phot_l95': r[2], 'ls_z_phot_u95': r[3],
                    'ls_fluxz': r[6],'ls_photoz_checked': check_int}
        info_out.append(info)

    return info_out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski.')

    parser.add_argument('--date-start', dest='date_start', type=str,
                        required=False,
                        help="Start date of the query, in ISO format. \
                        Example: '2017-08-17 12:41:04.4'", default=None)
    parser.add_argument('--date-end', dest='date_end', type=str,
                        required=False,
                        help="End date of the query, in ISO format. \
                        Example: '2017-08-18 12:00:00.0'", default=None)
    parser.add_argument('--min-days', dest='min_days', type=float,
                        required=False, help='Minimum time (days) between the \
                        first and last alert', default=20.)
    parser.add_argument('--max-days', dest='max_days', type=float,
                        required=False, help='Maximum time (days) between the \
                        first and last alert', default=100)
    parser.add_argument('--ndethist', dest='ndethist_min', type=int,
                        required=False,
                        help='Minimum number of detections', default=10)
    parser.add_argument('--out-query', dest='out', type=str, required=False,
                        help='Query output filename, txt',
                        default='results.txt')
    parser.add_argument('--out-lc', dest='out_lc', type=str, required=False,
                        help='Query output light curves (alerts+prv), CSV',
                        default='lightcurves.csv')
    parser.add_argument('--fields', dest='fields', type=str, required=False,
                        help='CSV file with a column of field names',
                        default='selected_fields_ebv01.csv')
    parser.add_argument("--v",  action='store_true',
                        help='Verbose: print out information on kowalski \
                        query status',
                        default=False)
    parser.add_argument("--doForcePhot", action="store_true",
                        default=False)
    parser.add_argument('--targetdir-base', dest='targetdir_base', type=str,
                        required=False,
                        help='Directory for the forced photometry',
                        default='./forced_photometry/')
    parser.add_argument("--doLCOSubmission",  action="store_true",
                        default=False)
    parser.add_argument("--doLCOStatus",  action="store_true",
                        default=False)
    parser.add_argument('--lco-programs', dest='lco_programs',
                        type=str, required=False,
                        default='NOAO2020B-005,TOM2020A-008')
    parser.add_argument("--doKNFit",  action="store_true", default=False)   
    parser.add_argument("--doCheckAlerts",  action="store_true",
                        default=False)
    parser.add_argument("--doWriteDb",  action='store_true',
                        help='Write information to the psql database \
                        (needs admin privileges)',
                        default=False)
    parser.add_argument("--doCLU",  action='store_true',
                        help='Crossmatch with the CLU galaxy catalog',
                        default=False)
    parser.add_argument("--path-CLU", dest='path_clu', type=str,
                        help='Path to the CLU galaxy catalog',
                        default='CLU_20190708_marshalFormat.hdf5')
    parser.add_argument('--path-secrets-db', dest='path_secrets_db', type=str,
                        required=False,
                        help="Path to the CSV file including the credentials \
                        to access the psql database", default='db_access.csv')
 
    args = parser.parse_args()

    # Selected fields
    if args.fields is not None:
        print(f"Using fields in {args.fields}")
        t = ascii.read(args.fields)
        list_fields = list(set(f for f in t['field'] if ((f > 156))))
    else:
        list_fields = np.arange(156,1900)

    # Send a warning if you need to have admin permissions
    if args.doWriteDb:
        print("WARNING! You activated a flag to write information \
into the database. If you are admin, this means that the database \
will be updated with the results of your queries.") 

    # Read the secrets
    secrets = ascii.read('./secrets.csv', format = 'csv')
    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    kow = Kowalski(username=username, password=password)
    ##connection_ok = kow.check_connection()
    ##if not connection_ok:
    ##    raise KowalskiError('not connected to Kowalski DB')
    ##print(f'Connection to Kowalski OK: {connection_ok}')

    # Iterate over a certain date range
    if args.date_start is None:
        date_start = Time.now() - datetime.timedelta(days=1)
    else:
        try:
            date_start = Time(args.date_start, format='iso')
        except ValueError:
            print("Invalid start date. It must be a string in ISO format.")
            print("Example: '2017-08-17 12:41:04.4'")
            exit()

    if args.date_end is None:
        date_end = Time.now()
    else:
        try:
            date_end = Time(args.date_end, format='iso')
        except ValueError:
            print("Invalid end date. It must be a string in ISO format.")
            print("Example: '2018-01-01 12:41:04.4'")
            exit()

    sources_kowalski_all = []
    jd_gap = date_end.jd - date_start.jd

    # If the gap is larger than thresh_days, pass a list of jds
    thresh_days = 30.
    if jd_gap < thresh_days:
        list_jd = [date_start.jd]
    else:
        list_jd = np.linspace(date_start.jd, date_end.jd,
                          int((date_end.jd - date_start.jd)/thresh_days)+1)
        jd_gap = list_jd[1] - list_jd[0] + 1

    print("Querying kowalski...")
    for jd in list_jd:    
        #Query kowalski
        sources_kowalski = query_kowalski(kow, list_fields,
                                          args.min_days, args.max_days,
                                          args.ndethist_min,
                                          jd, jd_gap=jd_gap,
                                          verbose=args.v)

        sources_kowalski_all += list(sources_kowalski)
    sources_kowalski_all = set(sources_kowalski_all)

    # Check full light curve duration (alerts)
    print("Checking durations.....")
    clean_set = check_lightcurve_alerts(username, password,
                                        sources_kowalski_all,
                                        args.min_days, args.max_days)
    print("...Done.")

    if clean_set is None:
        print(f"The Kowalski query did not return any candidate \
between {date_start} and {date_end} \
with the specified criteria.")
    else:
        print("Final set:")
        print(clean_set)
        print(f"Total: {len(clean_set)} candidates between {date_start.iso} \
and {date_end.iso}")

        #Print results to an output text file
        with open(args.out, 'a') as f:
            f.write(f"#{args} \n")
            f.write("name \n")
            for n in clean_set:
                f.write(f"{n} \n")

    # Get the light curves
    print("Getting light curves from the alerts...")
    from get_lc_kowalski import get_lightcurve_alerts, \
                                get_lightcurve_alerts_aux, create_tbl_lc

    # If there are candidates at all..
    if clean_set is not None:
        light_curves_alerts = get_lightcurve_alerts(username,
                                                    password,
                                                    clean_set)

        # Add prv_candidates photometry to the light curve
        print("Getting light curves from the alerts prv...")
        light_curves_aux = get_lightcurve_alerts_aux(username,
                                                     password,
                                                     clean_set)
    else:
        light_curves_alerts, light_curves_aux = None, None

    # Are there any candidates?
    if light_curves_alerts is not None and light_curves_aux is not None:
        light_curves = light_curves_alerts + light_curves_aux
    elif light_curves_alerts is not None:
        light_curves = light_curves_alerts
    elif light_curves_aux is not None:
        light_curves = light_curves_aux
    else:
        light_curves = None

    # Create a table and output CSV file
    if light_curves is not None:
        tbl_lc = create_tbl_lc(light_curves, outfile=args.out_lc)
    else:
        tbl_lc = None

    if tbl_lc is None:
        print("No candidates from the Kowalski query")
        print("Exiting...")
        exit()

    if args.doWriteDb:
        # Connect to the database
        con, cur = connect_database(update_database=args.doWriteDb,
                                    path_secrets_db=args.path_secrets_db,
                                    dbname='db_lens_admin')

    # Crossmatch with CLU
    names_matched, names_no_matched = xmatch_clu(tbl=tbl_lc,
               max_dist=200.,
               min_dist=20.,
               max_dist_kpc=130,
               min_dist_kpc=10,
               path_clu='CLU_20190708_marshalFormat.hdf5')

    # In the table, leave only candidates with a match
    index = []
    if names_matched is None:
        print("No matches found with CLU!")
        exit()
    for n in names_matched:
        index += list(np.where(tbl_lc['name'] == n)[0])
    print(f"There were {len(set(tbl_lc['name']))} candidates before crossmatch")
    tbl_lc = tbl_lc[index]
    print(f"There are {len(set(tbl_lc['name']))} candidates after crossmatch")

    # Write output
    jd_DR17 = Time("2021-11-06 00:00:00").jd
    for name in set(tbl_lc['name']):
        tbl_lc_short = tbl_lc[tbl_lc['name'] == name]
        # Check the data rights, DR17
        if np.max(np.array(tbl_lc_short["jd"])) > jd_DR17:
            tbl_lc_short = tbl_lc_short[tbl_lc_short["programid"] < 3]
        tbl_lc_short.write(f"lc_gap_csv/lc_{name}.csv", format='csv', overwrite=True)

    exit()

    # For the lensed project, no variability selection
    """
    # Select based on the variability criteria
    from select_variability_db import select_variability

    # Alerts
    if tbl_lc is not None:
        selected, rejected, cantsay = select_variability(tbl_lc,
                       hard_reject=[], update_database=args.doWriteDb,
                       read_database=True,
                       use_forced_phot=False, stacked=False,
                       baseline=0.5, var_baseline={'g': 6, 'r': 8, 'i': 10},
                       max_duration_tot=15., max_days_g=7., snr=4,
                       index_rise=-1.0, index_decay=0.3,
                       path_secrets_db=args.path_secrets_db,
                       save_plot=True, path_plot='./plots/',
                       show_plot=False, use_metadata=False,
                       path_secrets_meta='../kowalski/secrets.csv',
                       save_csv=True, path_csv='./lc_lens_csv',
                       path_forced='./')
    else:
        selected, rejected, cantsay = None, None, None

    # Check if the select_variability_db function returned any candidate
    if selected is not None:
        # which objects do we care about
        allids = selected + cantsay
        # select only relevant entries from the light curve table
        indexes = list(i for i, n in
                       zip(np.arange(len(tbl_lc)), tbl_lc['name'])
                       if n in allids)
        tbl_lc = tbl_lc[indexes]

    else:
        allids = []

    if args.doCheckAlerts and tbl_lc is not None:
        print("Checking alerts...")
        from alert_check import alert_check_complete
        ind_check_alerts = []
        for objid in allids:
            index_check = alert_check_complete(kow, objid)
            ind_check_alerts.append(index_check)
        ind_check_alerts = np.array(ind_check_alerts)
        allids = np.asarray(allids)[ind_check_alerts<2]

    # Check the database for candidates to do forced phot with
    # FIXME add argument to the arg parser
    read_database = True
    if read_database:
        # Connect to the database
        con, cur = connect_database(update_database=args.doWriteDb,
			            path_secrets_db=args.path_secrets_db,
                                    dbname='db_lens_admin')
        ####
        # Select from the db which candidates need forced photometry
        # In this offline version, select all and only those that
        # were found during the query
        ok_dur = allids 

        # Check which new candidates were already hard rejected
        names_str = "','".join(list(allids))
        cur.execute(f"select name from candidate \
where hard_reject = 1 and name in ('{names_str}')")
        r = cur.fetchall()
        # Bad ones, already rejected
        ko = list(l[0] for l in r)

        names_ok = list(n for n in ok_dur if
                        (not (n in ko)))
        candidates_for_phot = set(list(n for n in allids if
                                       not n in ko) + names_ok)
        # What if there are no candidates?
        if len(candidates_for_phot) == 0:
            print("There are no candidates do do forced photometry with")
            t_for_phot = None
        else:
            # Get the alerts light curve to improve the location accuracy
            lc_for_phot = get_lightcurve_alerts(username,
                                                password,
                                                candidates_for_phot)
            # Create a table in the right format
            t_for_phot = create_tbl_lc(lc_for_phot, outfile=None)
    """

    # For the lensed SN project, use all candidates found by kowalski
    t_for_phot = tbl_lc

    if args.doForcePhot and t_for_phot is not None:
        print("Triggering forced photometry...")
        from forcephot import trigger_forced_photometry

        # Trigger forced photometry
        success, _ = trigger_forced_photometry(t_for_phot,
                                               args.targetdir_base,
                                               daydelta_before=21.,
                                               daydelta_after=35.,
                                               ncpus=8)
        if args.doWriteDb and len(success) > 0:
            # Connect to the database
            con, cur = connect_database(update_database=args.doWriteDb,
                                        path_secrets_db=args.path_secrets_db,
                                        dbname='db_lens_admin')

            # Update the database with forced photometry
            from functions_db import populate_table_lightcurve_forced
            populate_table_lightcurve_forced(con, cur, t_for_phot,
                                             args.targetdir_base,
                                             programids=[1,2])
            print("POPULATED forced photometry table")

            # Update the database with stacked forced photometry
            from functions_db import populate_table_lightcurve_stacked
            populate_table_lightcurve_stacked(con, cur, success)
            print("POPULATED stacked forced photometry table")

            # Close the connection to the db
            cur.close()
            con.close()
    # Again, for the lensed SN project there is no need of selection
    """
    if t_for_phot is not None:
        # Repeat the selection based on forced photometry
        selected, rejected, cantsay = select_variability(t_for_phot,
                   hard_reject=[], update_database=args.doWriteDb,
                   read_database=True,
                   use_forced_phot=True, stacked=False,
                   baseline=0.5, var_baseline={'g': 6, 'r': 8, 'i': 10},
                   max_duration_tot=15., max_days_g=7., snr=4,
                   index_rise=-1.0, index_decay=0.3,
                   path_secrets_db=args.path_secrets_db,
                   save_plot=True, path_plot='./plots/',
                   show_plot=False, use_metadata=False,
                   path_secrets_meta='../kowalski/secrets.csv',
                   save_csv=True, path_csv='./lc_csv',
                   path_forced='./')

        # Repeat the selection based on stacked forced photometry
        selected, rejected, cantsay = select_variability(t_for_phot,
                   hard_reject=[], update_database=args.doWriteDb,
                   read_database=True,
                   use_forced_phot=True, stacked=True,
                   baseline=0.5, var_baseline={'g': 6, 'r': 8, 'i': 10},
                   max_duration_tot=15., max_days_g=7., snr=4,
                   index_rise=-1.0, index_decay=0.3,
                   path_secrets_db=args.path_secrets_db,
                   save_plot=True, path_plot='./plots/',
                   show_plot=False, use_metadata=False,
                   path_secrets_meta='../kowalski/secrets.csv',
                   save_csv=True, path_csv='./lc_csv',
                   path_forced='./')

    # Populate the database with CLU galaxy catalog crossmatch
    if args.doCLU:
        if args.doWriteDb:
            # Connect to the database
            con, cur = connect_database(update_database=args.doWriteDb,
                                        path_secrets_db=args.path_secrets_db,
                                        dbname='db_lens_admin')

            # Import the relevant function
            from functions_db import populate_table_clu
            populate_table_clu(con, cur, tbl=None,
                               max_dist=100.,
                               path_clu=args.path_clu)
            print("POPULATED CLU crossmatch table")
        else:
            print("WARNING: in order to do the CLU galaxy catalog \
crossmatching, you need to have --doWriteDb active")
        
    if args.doKNFit:
        print('Fitting to kilonova grid...')

        from knfit import do_knfit
        for objid in allids:
            t = tbl_lc[tbl_lc['name'] == objid]
            do_knfit(t.to_pandas().rename(columns={"filter": "filtname"}))

    if args.doLCOStatus:
        print('Checking LCO for existing observations...')

        # LCO sometime over next 2 weeks
        tstart = Time.now()
        tend = Time.now() + TimeDelta(14*u.day)
        tstart = str(tstart.isot).replace("T"," ")
        tend = str(tend.isot).replace("T"," ")

        #Read the secrets
        lco_secrets = ascii.read('../lco/secrets.csv', format = 'csv')
        PROPOSAL_ID = lco_secrets['PROPOSAL_ID'][0]
        API_TOKEN = lco_secrets['API_TOKEN'][0]

        lco_programs = args.lco_programs.split(",")

        from lco import check_observations
        obs = check_observations(API_TOKEN, lco_programs=lco_programs)

    if args.doLCOSubmission: 
        print('Triggering LCO...')

        # LCO sometime over next 2 weeks
        tstart = Time.now() 
        tend = Time.now() + TimeDelta(14*u.day)    
        tstart = str(tstart.isot).replace("T"," ")
        tend = str(tend.isot).replace("T"," ")
    
        #Read the secrets
        lco_secrets = ascii.read('../lco/secrets.csv', format = 'csv')
        PROPOSAL_ID = lco_secrets['PROPOSAL_ID'][0]
        API_TOKEN = lco_secrets['API_TOKEN'][0]
    
        from lco import submit_photometric_observation
        from lco import submit_spectroscopic_observation

        for objid in allids:

            if args.doLCOStatus:
                to_observe = True
                for key in obs:
                    if obs[key]["completed"] == 1: # PENDING
                        to_observe = False
                if not to_observe:
                    continue

            t = tbl_lc[tbl_lc['name'] == objid]
            ra, dec = np.median(t['ra']), np.median(t['dec'])
   
            submit_photometric_observation(objid, ra, dec,
                                           PROPOSAL_ID, API_TOKEN,
                                           tstart=tstart, tend=tend,
                                           exposure_time = 300,
                                           doSubmission=False)

            submit_spectroscopic_observation(objid, ra, dec,
                                             PROPOSAL_ID, API_TOKEN,
                                             tstart=tstart, tend=tend,
                                             exposure_time = 300,
                                             doSubmission=False)
    """
    print("Done.")
