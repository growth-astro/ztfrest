'''
Query Kowalski searching for transients
given a set of constraints.
'''

import json
import requests
import datetime
import pdb

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time, TimeDelta

from penquins import Kowalski
from functions_db import connect_database


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


def check_lightcurve_alerts(username, password, list_names, min_days, max_days):
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
    if r['result_data']['query_result'] == []:
        print("No candidates to be checked?")
        return None

    old = []
    objectid_list = []
    for info in r['result_data']['query_result']:
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
                   ndethist_min, jd, jd_gap=50., verbose=True):
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
                                  'candidate.drb': {'$gt': 0.9},
                                  'classifications.braai': {'$gt': 0.8},
                                  'candidate.ndethist': {'$gt': ndethist_min_corrected},
                                  'candidate.magpsf': {'$gt': 12}
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
                if r['result_data']['query_result'] == []:
                    no_candidates = True
                break
            except KeyError:
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

        for info in r['result_data']['query_result']:    
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
                        first and last alert', default=0.01)
    parser.add_argument('--max-days', dest='max_days', type=float,
                        required=False, help='Maximum time (days) between the \
                        first and last alert', default=14.)
    parser.add_argument('--ndethist', dest='ndethist_min', type=int,
                        required=False,
                        help='Minimum number of detections', default=2)
    parser.add_argument('--out-query', dest='out', type=str, required=False,
                        help='Query output filename, txt',
                        default='results.txt')
    parser.add_argument('--out-lc', dest='out_lc', type=str, required=False,
                        help='Query output light curves (alerts+prv), CSV',
                        default='lightcurves.csv')
    parser.add_argument('--fields', dest='fields', type=str, required=False,
                        help='CSV file with a column of field names',
                        default=None)
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
    connection_ok = kow.check_connection()
    if not connection_ok:
        raise KowalskiError('not connected to Kowalski DB')
    print(f'Connection to Kowalski OK: {connection_ok}')

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

    light_curves_alerts = get_lightcurve_alerts(username,
                                                password,
                                                clean_set)

    # Add prv_candidates photometry to the light curve
    print("Getting light curves from the alerts prv...")
    light_curves_aux = get_lightcurve_alerts_aux(username,
                                                 password,
                                                 clean_set)

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

    if args.doWriteDb and tbl_lc is not None:
        # Connect to the database
        con, cur = connect_database(update_database=args.doWriteDb,
                                    path_secrets_db=args.path_secrets_db)

        # Add the candidates to the db
        from functions_db import populate_table_candidate
        populate_table_candidate(tbl_lc, con, cur)
        print("POPULATED candidate table")

        # Upload the light curves to the database
        from functions_db import populate_table_lightcurve
        populate_table_lightcurve(tbl_lc, con, cur)
        print("POPULATED alert lightcurves")

        # Extinction information
        from functions_db import populate_extinction
        populate_extinction(con, cur)

        # Galactic latitude
        from functions_db import populate_gal_lat
        populate_gal_lat(con, cur)

        con.close()
        cur.close()

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
                       save_csv=True, path_csv='./lc_csv',
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
			            path_secrets_db=args.path_secrets_db)
        ####
        # Select from the db which candidates need forced photometry
        cur.execute("select name from candidate where \
(duration_tot < 21 or duration_tot is null) and \
(hard_reject is NULL or hard_reject = 0)")
        r = cur.fetchall()
        # OK for duration
        ok_dur = list(l[0] for l in r) 

        cur.execute(f"select name from lightcurve \
where jd > {date_end.jd - 14}")
        r = cur.fetchall()
        # OK for alerts light curve
        ok_lc = list(l[0] for l in r) 

        cur.execute(f"select name from lightcurve_forced \
where jd > {date_end.jd - 14}")
        r = cur.fetchall()
        # OK for forced phot light curve
        ok_lc_forced = list(l[0] for l in r)

        # Check which new candidates were already hard rejected
        names_str = "','".join(list(allids))
        cur.execute(f"select name from candidate \
where hard_reject = 1 and name in ('{names_str}')")
        r = cur.fetchall()
        # Bad ones, already rejected
        ko = list(l[0] for l in r)

        names_ok = list(n for n in ok_dur if
                        ((n in ok_lc or n in ok_lc_forced) and not (n in ko)))
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

    if args.doForcePhot and t_for_phot is not None:
        print("Triggering forced photometry...")
        from forcephot import trigger_forced_photometry

        # Trigger forced photometry
        success, _ = trigger_forced_photometry(t_for_phot,
                                               args.targetdir_base,
                                               daydelta_before=7.,
                                               daydelta_after=14.)

        if args.doWriteDb and len(success) > 0:
            # Update the database with forced photometry
            from functions_db import populate_table_lightcurve_forced
            populate_table_lightcurve_forced(con, cur, t_for_phot,
                                             args.targetdir_base)
            print("POPULATED forced photometry table")

            # Update the database with stacked forced photometry
            from functions_db import populate_table_lightcurve_stacked
            populate_table_lightcurve_stacked(con, cur, success)
            print("POPULATED stacked forced photometry table")

            # Close the connection to the db
            cur.close()
            con.close()

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
                                        path_secrets_db=args.path_secrets_db)

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

    print("Done.")
