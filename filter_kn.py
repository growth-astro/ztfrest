'''
Query Kowalski searching for transients
given a set of constraints.

Selection based on alert+prv data, 
                   alert forced photometry, 
                   custom forced photometry (FPZTF pipeline),
                   and binned forced photometry.
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
from astropy.table import vstack, Table

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

    r = k.query(query=q) # response
    if r.get("default").get("data") == []:
        print("No candidates to be checked.")
        return None

    old = []
    objectid_list = []
    
    for info in r.get("default").get("data"):
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
    
    returned_by_query = []
    
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
                                  'candidate.magpsf': {'$gt': 12},
                                  'candidate.isdiffpos': 't'
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
                if r.get("default").get("data") == []:
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
        
        
        returned_by_query += [d['objectId'] for d in r['default']['data'] if d['objectId'] not in returned_by_query]
        
        for info in r.get("default").get("data"):
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
        
        set_objectId_all = set_objectId_all | set_objectId
        #print(f"Cumulative: {len(set_objectId_all)}")

        if verbose is True:
            print("Field", field, len(set_objectId_all))

    print(f"Returned by initial query: {len(set(returned_by_query))}")
    
    return set_objectId_all

###################################################################################################
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
    #parser.add_argument('--days-starthist', dest='days_starthist', type=float, required=False, 
                        #help='Maximum number of days to look back for ndethist', default=30.)
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
    parser.add_argument("--doForcePhot", action="store_true", help="Trigger forced photometry using ForcePhotZTF",
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
    ### new
    parser.add_argument("--dontReadDb", action="store_true", help="Do not read information \
                        from the psql database during selection. Retrieves photometry from kowalski instead", 
                        default = False)
    parser.add_argument("--doAuxFp", action="store_true", 
                        help='Use forced photometry from the \
                        ZTF alert packets', default=False)
    args = parser.parse_args()
# ---------------------------------------------------------------------------------------------
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
    connection_ok = kow.ping()
    if not connection_ok:
       raise KowalskiError('Error: Not connected to Kowalski')
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

    print("Querying Kowalski...")
    for jd in list_jd:    
        # #Query kowalski, return a list of candidate names
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
        #print("Final set:")
        #print(clean_set)
        print(f"Total: {len(clean_set)} candidates found between {date_start.iso} \
and {date_end.iso}")

        #Print results to an output text file
        with open(args.out, 'w') as f:
            f.write(f"#{args} \n")
            f.write("name \n")
            for n in clean_set:
                f.write(f"{n} \n")
        print("----------------------------------------")

#### Get the light curves from Kowalski #################
    print("Getting light curves from the alerts...")
    from get_lc_kowalski import get_lightcurve_alerts, \
                                get_lightcurve_alerts_aux, create_tbl_lc, \
                                get_fphists_alerts_aux, create_tbl_lc_fphists
    
    use_fphists = args.doAuxFp # use forced photometry from alerts?
    
    # If there are candidates at all:
    if clean_set is not None:
        # Get alert data from ZTF_alerts catalog
        light_curves_alerts = get_lightcurve_alerts(username,
                                                    password,
                                                    clean_set) # (SNR >= 5)

        # Add prv_candidates (previous) photometry from ZTF_alerts_aux catalog
        print("Getting light curves from the alerts prv...")
        light_curves_prv = get_lightcurve_alerts_aux(username,
                                                     password,
                                                     clean_set) # (3<=SNR<5)
                                                     
        if use_fphists is True:
            # Add fp_hists (forced photometry from the alerts)
            print("Getting forced photometry from the alerts...")
            light_curves_fp = get_fphists_alerts_aux(username,password,clean_set)
        
    else:
        light_curves_alerts, light_curves_prv, light_curves_fp = None, None, None
        
    # Create a table of alert photometry + alert forced photometry
    if use_fphists is True:
        if light_curves_fp is not None and light_curves_alerts is not None:
            tbl_only_alerts = create_tbl_lc(light_curves_alerts, origin='alert')
            tbl_only_fphists = create_tbl_lc_fphists(light_curves_fp)
            tbl_only_prv = create_tbl_lc(light_curves_prv,origin='prv')
            
            tbllsalerts = []
            for name in clean_set:
                ta = tbl_only_alerts[tbl_only_alerts['name']==name] # per name
                ta.add_index('jd')
                tf = tbl_only_fphists[tbl_only_fphists['name']==name] # per name
                # For each candidate, remove alert datapoint if it has alert FP on the same JD
                # (FP tend to have lower unc)
                for jd in list(ta['jd']):
                    if jd in list(tf['jd']):
                        ta.remove_row(ta.loc_indices[jd])
                tbllsalerts.append(ta)
            tbl_only_alerts = vstack(tbllsalerts)
            
            # For each alert FP point, designate as an upper limit if SNR < 3. Otherwise - detection.
            for row in tbl_only_fphists:
                if row['snr'] < 3.0:
                    row['origin'] = 'magul' # mag upper limit
 
            tbl_alert_fp = vstack([tbl_only_alerts,tbl_only_fphists])
            tbl_alert_fp.sort('jd')
            
        elif light_curves_fp is not None:
            tbl_alert_fp = create_tbl_lc_fphists(light_curves_fp) #already sorted by jd
            for row in tbl_alert_fp:
                if row['snr'] < 3.:
                    row['origin'] = 'magul'
        else:
            tbl_alert_fp = None  
            
    # Create table of alert photometry + prv_candidates photometry, keep track of origin
    # Are there any candidates?
    if light_curves_alerts is not None:
        tbl_alerts = create_tbl_lc(light_curves_alerts, origin='alert')
    else:
        tbl_alerts = None
    if light_curves_prv is not None:
        tbl_prv = create_tbl_lc(light_curves_prv, origin='prv')
    else:
        tbl_prv = None
    
    if tbl_alerts is not None and tbl_prv is not None:
        tbl_lc = vstack([tbl_alerts, tbl_prv])
    elif tbl_alerts is not None:
        tbl_lc = tbl_alerts
    elif tbl_prv is not None:
        tbl_lc = tbl_prv
    else:
        tbl_lc = None
        print("No lightcurves")
    if tbl_lc is not None:
        tbl_lc.sort('jd')


    if args.doWriteDb and tbl_lc is not None:
        print("----------------------------------------")
        print("Updating database...")
        # Connect to the database
        con, cur = connect_database(update_database=args.doWriteDb,
                                    path_secrets_db=args.path_secrets_db, 
                                    dbname='db_kn_2025_admin')
        # make sure you are writing to the correct database
        cur.execute("select current_database()")
        current_dbname = cur.fetchall()[0][0]
        if current_dbname != 'database_kn_2025':
            print("Wrong database!!! terminating connection...")
            con.close()
            cur.close()
        # Add the candidates to the db
        from functions_db import populate_table_candidate
        populate_table_candidate(tbl_lc, con, cur)
        print("POPULATED candidate table")

        # Upload the light curves to the database
        from functions_db import populate_table_lightcurve
        populate_table_lightcurve(tbl_lc, con, cur)
        print("POPULATED alert lightcurves")
        
        # populate DB with alert forced photometry 
        from functions_db import populate_table_lightcurve_alertfp
        populate_table_lightcurve_alertfp(tbl=tbl_alert_fp, con=con, cur=cur)
        print("POPULATED alert forced photometry lightcurves")
        
        # Extinction information
        from functions_db import populate_extinction
        populate_extinction(con, cur)
        print("POPULATED extinctions")

        # Galactic latitude
        from functions_db import populate_gal_lat
        populate_gal_lat(con, cur)
        print("POPULATED galactic latitudes")

        con.close()
        cur.close()
    
    # Get photometry from database or directly from Kowalski?
    if args.dontReadDb:
        read_database = False
    else:
        read_database = True
# -----------------------------------------------------------------------------------------
    # Select based on the variability criteria
    from select_variability_db import select_variability
    # Alert photometry
    if tbl_lc is not None:
        print("--------------------------------------------")
        print("Selecting based on alerts + prv photometry.....")
        selected, rejected, cantsay = select_variability(tbl_lc,
                       hard_reject=[], update_database=args.doWriteDb,
                       read_database=read_database,
                       use_forcephotztf=False, stacked=False,
                       baseline=0.125, var_baseline={'g': 6, 'r': 8, 'i': 10},
                       max_duration_tot=15., max_days_g=7., snr=4,
                       index_rise=-1.0, index_decay=0.3,
                       path_secrets_db=args.path_secrets_db,
                       save_plot=True, path_plot='./plots/',
                       show_plot=False, use_metadata=False,
                       path_secrets_meta='./secrets.csv',
                       save_csv=True, path_csv='./lc_csv',
                       path_forced='./')
    else:
        selected, rejected, cantsay = None, None, None
# -----------------------------------------------------------------------------------------
    if use_fphists:
        # Repeat the selection based on forced photometry from alerts
        print("----------------------------------------------")
        print("Selecting based on alert forced photometry...")
        if tbl_alert_fp is not None:
            selected, rejected, cantsay = select_variability(tbl_alert_fp,
                           hard_reject=[], update_database=False,
                           read_database=read_database,
                           use_forcephotztf=False, use_alertfp=True, stacked=False,
                           baseline=0.125, var_baseline={'g': 6, 'r': 8, 'i': 10},
                           max_duration_tot=15., max_days_g=7., snr=4,
                           index_rise=-1.0, index_decay=0.3,
                           path_secrets_db=args.path_secrets_db,
                           save_plot=True, path_plot='./plots/',
                           show_plot=False, use_metadata=False,
                           path_secrets_meta='./secrets.csv',
                           save_csv=True, path_csv='./lc_csv',
                           path_forced='./') # important that use_forcephotztf is False
        else:
            selected, rejected, cantsay = None, None, None
        
        print("----------------------------------------------")
################################################################################################
    # Determine which candidates need to be passed on to the ForcePhotZTF pipeline.
    # Get candidates that have lightcurves and/or ForcePhotZTF data in the kn database
    
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
        print("----------------------------------------")
        print("Checking alerts...")
        from alert_check import alert_check_complete
        ind_check_alerts = []
        for objid in allids:
            index_check = alert_check_complete(kow, objid)
            ind_check_alerts.append(index_check)
        ind_check_alerts = np.array(ind_check_alerts)
        allids = np.asarray(allids)[ind_check_alerts<2]

    # Check the database for candidates to do forced phot with
    # FIXME add argument to the arg parser --> done

    if read_database:
        print("----------------------------------------------------------------------")
        print("Reading database for lightcurves and ForcePhotZTF forced photometry...")
        # Connect to the database
        con, cur = connect_database(update_database=args.doWriteDb,
			            path_secrets_db=args.path_secrets_db)
        cur.execute("select current_database()")
        current_dbname = cur.fetchall()[0][0]
        if args.doWriteDb is True and current_dbname != 'database_kn_2025':
            print("panic")
            con.close()
            cur.close()
        ####
        # Select from the db which candidates need forced photometry
        # this part relies on the database having been updated by someone 
        cur.execute("select name from candidate where \
(duration_tot < 21 or duration_tot is null) and \
(hard_reject is NULL or hard_reject = 0)")
        r = cur.fetchall()
        # OK for duration
        ok_dur = list(l[0] for l in r)
        print(f"durations Ok: {len(ok_dur)}")

        cur.execute(f"select name from lightcurve where jd between {date_end.jd - 14} and {date_end.jd}")
        r = cur.fetchall()
        # OK for alerts light curve
        ok_lc = list(l[0] for l in r) 
        print("Light curves from alerts from kn database (< 14 days old):", len(ok_lc))

        cur.execute(f"select name from lightcurve_forced where jd between {date_end.jd - 14} and {date_end.jd}")
        r = cur.fetchall()
        # OK for forced phot light curve
        ok_lc_forced = list(l[0] for l in r)
        #print("Forced photometry light curves from kn database (< 14 days old):", len(ok_lc_forced))

        # Check which new candidates were already hard rejected
        names_str = "','".join(list(allids))
        cur.execute(f"select name from candidate \
where hard_reject = 1 and name in ('{names_str}')")
        r = cur.fetchall()
        # Bad ones, already rejected
        ko = list(l[0] for l in r)
        #print("Already rejected in kn database:", len(ko))

        names_ok = list(n for n in ok_dur if
                        ((n in ok_lc or n in ok_lc_forced) and not (n in ko)))
        candidates_for_phot = set(list(n for n in allids if
                                       not n in ko) + names_ok) 
        print("Candidates needing forced photometry:", len(candidates_for_phot))
        # What if there are no candidates?
        if len(candidates_for_phot) == 0:
            print("There are no candidates to do forced photometry with.")
            t_for_phot = None
        else:
            # Get the alerts light curve to improve the location accuracy
            lc_for_phot = get_lightcurve_alerts(username,
                                                password,
                                                candidates_for_phot)
            # Create a table in the right format
            # record candidates needing ForcePhotZTF, so it can be triggered separately -- have inputted outfile argument here
            t_for_phot = create_tbl_lc(lc_for_phot, outfile='../ztfrest/needs_fp.csv', origin = 'fpztf')
    else:
        t_for_phot = tbl_lc
        if tbl_lc is not None:
            candidates_for_phot = set(t_for_phot['name'])
        else:
            candidates_for_phot = None
    
    
    
    # Populate the database with CLU galaxy catalog crossmatch
    if args.doCLU:
        if args.doWriteDb:
            # Connect to the database
            con, cur = connect_database(update_database=args.doWriteDb,
                                        path_secrets_db=args.path_secrets_db, dbname='db_kn_2025_admin')
            cur.execute("select current_database()")
            current_dbname = cur.fetchall()[0][0]
            if current_dbname != 'database_kn_2025':
                print("panic")
                con.close()
                cur.close()

            # Import the relevant functions
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
                                             
    # Run the ForcePhotZTF pipeline on candidates needing forced photometry                                         
    if args.doForcePhot is True and t_for_phot is not None:
        print("----------------------------------------")
        print("Triggering forced photometry...")
        print("WARNING: ForcePhotZTF requires specific environment setup and older package versions. \
               See https://github.com/yaoyuhan/ForcePhotZTF/blob/master")
        from forcephot import trigger_forced_photometry

        # Trigger forced photometry
        try:
            success, _ = trigger_forced_photometry(t_for_phot,
                                                   args.targetdir_base,
                                                   daydelta_before=7.,
                                                   daydelta_after=14.)
        except:
            success, _ = [0], None
        if len(success) == 0:
            print("Failed to run ForcePhotZTF.")
        if args.doWriteDb and len(success) > 0:
            print("----------------------------------------")
            print("Updating database...")
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
##############################################################################
    #if t_for_phot is not None:
    if args.doForcePhot is True and t_for_phot is not None:
        # Repeat the selection based on ForcePhotZTF data 
        print("----------------------------------------------")
        print("Repeating based on ForcePhotZTF photometry...")
        selected, rejected, cantsay = select_variability(t_for_phot,
                   hard_reject=[], update_database=args.doWriteDb,
                   read_database=read_database,
                   use_forcephotztf=True, stacked=False,
                   baseline=0.125, var_baseline={'g': 6, 'r': 8, 'i': 10},
                   max_duration_tot=15., max_days_g=7., snr=4,
                   index_rise=-1.0, index_decay=0.3,
                   path_secrets_db=args.path_secrets_db,
                   save_plot=True, path_plot='./plots/',
                   show_plot=False, use_metadata=False,
                   path_secrets_meta='./secrets.csv',
                   save_csv=True, path_csv='./lc_csv',
                   path_forced='./')
            
        # Repeat the selection based on stacked forced photometry
        print("-----------------------------------------------")
        print("Repeating based on stacked forced photometry...")
        selected, rejected, cantsay = select_variability(t_for_phot,
                   hard_reject=[], update_database=args.doWriteDb,
                   read_database=read_database,
                   use_forcephotztf=True, stacked=True,
                   baseline=0.125, var_baseline={'g': 6, 'r': 8, 'i': 10},
                   max_duration_tot=15., max_days_g=7., snr=4,
                   index_rise=-1.0, index_decay=0.3,
                   path_secrets_db=args.path_secrets_db,
                   save_plot=True, path_plot='./plots/',
                   show_plot=False, use_metadata=False,
                   path_secrets_meta='./secrets.csv',
                   save_csv=True, path_csv='./lc_csv',
                   path_forced='./')  
    
    if read_database:
        con.close()
        cur.close()
    print("Done.")
