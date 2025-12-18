"""
Originally developed at
https://github.com/igorandreoni/kowalski-searches/blob/master/get_lc_kowalski.py

modified by Sarah Teague (steague1@unc.edu)
"""

__author__ = "Igor Andreoni"
__license__ = "MIT, GNU General Public License v3.0"
__email__ = "andreoni@caltech.edu"


from astropy.io import ascii
from astropy.table import Table, unique
import numpy as np
import pdb

from penquins import Kowalski


def get_lightcurve_alerts_aux(username, password, list_names):
    """Query Kowalski for a list of candidates"""
    """ Get 3 < SNR < 5 detections from previous fields in ZTF_alerts_aux"""

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts_aux",
                   "filter": {
                              '_id': {'$in': list(list_names)}
                              },
                   "projection": {}
                       },
         "kwargs": {"hint": "_id_"}
         }

    r = k.query(query=q)

    if r.get("default").get("data") == []:
        print("No candidates to be checked?")
        return None
    out = []
    for l in r.get("default").get("data"):
        if 'prv_candidates' in l.keys():
            with_det = list({'objectId': l['_id'], 'candidate': s} for s in l['prv_candidates'] if 'magpsf' in s.keys())
            out = out + with_det
    return out

def get_fphists_alerts_aux(username, password, list_names):
    """Given a list of candidates, obtain forced photometry from the alerts.
       Updated schema vsns >= 4 now have forced phot as part of the alert packets."""
    k = Kowalski(username=username, password=password, verbose=False)
    proj = {'_id':1, 'fp_hists.jd':1, 'fp_hists.alert_ra':1, 'fp_hists.alert_dec':1,
            'fp_hists.mag':1, 'fp_hists.magerr':1, 'fp_hists.fid':1, 'fp_hists.magzpsci':1,
            'fp_hists.magzpsciunc':1, 'fp_hists.programid':1, 'fp_hists.field':1,
            'fp_hists.rcid':1, 'fp_hists.pid':1, 'fp_hists.limmag5sig':1, 'fp_hists.diffmaglim':1,
            'fp_hists.snr':1, 'fp_hists.forcediffimflux':1, 'fp_hists.forcediffimfluxunc':1}
    q = {"query_type":"find", "query": {"catalog":"ZTF_alerts_aux", 
                                        "filter":{'_id': {'$in': list(list_names)}},
                                        "projection": proj},"kwargs":{"hint": "_id_"}}
    r = k.query(query=q)
    if r.get("default").get("data") == []:
        print("No candidates to be checked?")
        return None
    out = []
    for l in r.get("default").get("data"):
        if 'fp_hists' in l.keys():
            with_det = list({'objectId': l['_id'], 'fp_hists': s} for s in l['fp_hists'])
            out = out + with_det

    return out

def get_lightcurve_alerts(username, password, list_names):
    """Query the light curve for a list of candidates"""
    """ Get the >=5sig alerts from ZTF_alerts"""

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
                                  "candidate.ra": 1,
                                  "candidate.dec": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.isdiffpos": 1,
                                  "candidate.fid": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.programid": 1,
                                  "candidate.magzpsci": 1,
                                  "candidate.magzpsciunc": 1,
                                  "candidate.sgscore1": 1,
                                  "candidate.sgscore2": 1,
                                  "candidate.sgscore3": 1,
                                  "candidate.distpsnr1": 1,
                                  "candidate.distpsnr2": 1,
                                  "candidate.distpsnr3": 1,
                                  "candidate.field": 1,
                                  "candidate.rcid": 1,
                                  "candidate.pid": 1
                                  }
                       },
         "kwargs": {"hint": "objectId_1"}
         }

    r = k.query(query=q)
    try:
        if r.get("default").get("data") == []:
            print("No candidates to be checked?")
            return None
    except KeyError:
        #Try the query one more time
        r = k.query(query=q)
        try:
            if r.get("default").get("data") == []:
                print("No candidates to be checked?")
                return None
        except KeyError:
            return None
    return r.get("default").get("data")


def create_tbl_lc(light_curves, outfile=None, origin='alert'):
    """Create a table with the light curves
    and write a CSV output file"""
    """ To be used for alert and prv photometry"""

    # fid -> filter
    filters = {'1': 'g', '2': 'r', '3': 'i'}

    tbl = Table([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [],
                 [], [], [], [], [], []],
                names=('name', 'ra', 'dec', 'jd', 'isdiffpos', 'magpsf', 'sigmapsf',
                       'filter', 'magzpsci', 'magzpsciunc',
                       'programid', 'field', 'rcid', 'pid',
                       'sgscore1', 'sgscore2', 'sgscore3',
                       'distpsnr1', 'distpsnr2', 'distpsnr3', 'origin'),
                dtype=('S12', 'double', 'double', 'double', 'S',
                       'f', 'f', 'S', 'f', 'f', 'i', 'i', 'i', 'int_',
                       'f', 'f', 'f', 'f', 'f', 'f', 'S'))
    jd_done = []
    for l in light_curves:
        if (l["objectId"], l["candidate"]["jd"]) in jd_done:
            continue
        else:
            jd_done.append((l["objectId"], l["candidate"]["jd"]))
        magzpsci = l["candidate"].get("magzpsci")
        magzpsciunc = l["candidate"].get("magzpsciunc")
        try:
            row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["isdiffpos"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"], l["candidate"]["field"],
               l["candidate"]["rcid"], np.array(l["candidate"]["pid"]).astype('uint64'), 
               l["candidate"]["sgscore1"], l["candidate"]["sgscore2"],
               l["candidate"]["sgscore3"], l["candidate"]["distpsnr1"],
               l["candidate"]["distpsnr2"], l["candidate"]["distpsnr3"], origin]
        except KeyError:
            row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["isdiffpos"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"], l["candidate"]["field"],  
               l["candidate"]["rcid"], np.array(l["candidate"]["pid"]).astype('uint64'), np.nan,
               np.nan, np.nan, np.nan, np.nan, np.nan, origin]
        #print(row[13].dtype)
        tbl.add_row(row)
    # Remove exact duplicates
    tbl = unique(tbl)
    tbl.sort("jd")
    if outfile is not None:
        tbl.write(outfile, format='csv', overwrite=True)
        print("LCs written to csv")

    return tbl

def create_tbl_lc_fphists(light_curves, outfile=None):
    """ Specific to fp_hists - forced photometry from ZTF_alerts_aux """
    
    filters = {'1': 'g', '2': 'r', '3': 'i'}
    # same headers so this table can be stacked w/ the others; but schema are different 
    tbl = Table([[], [], [], [], [], [], [], [],
                 [], [], [], [], [], [], [], [], [], []],
                names=('name', 'ra', 'dec', 'jd', 'magpsf', 'sigmapsf', 'limmag', 'snr',
                       'filter', 'magzpsci', 'magzpsciunc',
                       'programid', 'field', 'rcid', 'pid',
                       'forcediffimflux', 'forcediffimfluxunc', 'origin'),
                dtype=('S12', 'double', 'double', 'double',
                       'f', 'f', 'f', 'f', 'S', 'f', 'f', 'i', 'i', 'i', 'int_',
                       'f', 'f', 'S'))
    jd_done = []
    for l in light_curves:
        if (l["objectId"], l["fp_hists"]["jd"]) in jd_done:
            continue
        else:
            jd_done.append((l["objectId"], l['fp_hists']["jd"]))
        magzpsci = l['fp_hists'].get("magzpsci")
        magzpsciunc = l['fp_hists'].get("magzpsciunc")
        mag = l['fp_hists']['mag'] if 'mag' in l['fp_hists'].keys() else np.nan
        magerr = l['fp_hists']['magerr'] if 'magerr' in l['fp_hists'].keys() else np.nan
        flux = l['fp_hists']['forcediffimflux'] if 'forcediffimflux' in l['fp_hists'].keys() else np.nan
        fluxunc = l['fp_hists']['forcediffimfluxunc'] if 'forcediffimfluxunc' in l['fp_hists'].keys() else np.nan

        if 'snr' in l['fp_hists'].keys():
            snr = l['fp_hists']['snr']
        else: # ??
            continue
        #if snr < 0:
            #continue
        # separate out upper limits; define if snr < 3. (could be adjusted)
        if snr < 3.:
            origin = 'magul'
        else:
            origin = 'alertfp'
        limmag = l['fp_hists']['limmag5sig'] if 'limmag5sig' in l['fp_hists'].keys() else l['fp_hists']['diffmaglim']
        
        row = [l["objectId"], l['fp_hists']["alert_ra"], l['fp_hists']["alert_dec"],
            l['fp_hists']["jd"], mag,
            magerr, limmag, snr, filters[str(l['fp_hists']["fid"])],
            magzpsci, magzpsciunc,
            l['fp_hists']["programid"], l['fp_hists']["field"],
            l['fp_hists']["rcid"], np.array(l['fp_hists']["pid"]).astype('uint64'),
            flux, fluxunc, origin]
        
        tbl.add_row(row)
    # Remove exact duplicates
    tbl = unique(tbl)
    tbl.sort("jd")
    if outfile is not None:
        tbl.write(outfile, format='csv', overwrite=True)

    return tbl

def create_tbl_alert_fp(lcs_alerts, lcs_fphists, outfile=None):
    '''create a table of alert photometry + alert forced photometry.
       include fp with SNR < 3 as upper limits (check threshold) 
       
       current cols: name, ra, dec, jd, magpsf, sigmapsf, limmag5sig,
                     snr, filter, magzpsci, magzpsciunc, programid, 
                     field, rcid, pid, sgscore1, sgscore2, sgscore3,
                     distpsnr1, distpsnr2, distpsnr3, origin     '''
    tbl_alerts = create_tbl_lc(lcs_alerts)
    tbl_fphists = create_tbl_lc_fphists(lcs_fphists)
    
    list_names = list(tbl_alerts['name'])
    tbllsalerts = []
    for name in list_names:
        talert = tbl_alerts[tbl_alerts['name'] == name]
        talert.add_index('jd')
        tfp = tbl_fphists[tbl_fphists['name']== name]
        # check whether there is an alert that has fp on the same jd; if so, remove alert measurement
        for jd in list(talert['jd']):   
            if jd in list(tfp['jd']): 
                #print(f"Removing JD {jd} in alerts for {name} (duplicate jd)")
                talert.remove_row(talert.loc_indices[jd])
        tbllsalerts.append(talert)
    tbl_alerts = vstack(tbllsalerts)
    
    # Designate non-significant FP as upper limit (if SNR < 3). Otherwise, detection: treat it the same as an alert
    for row in tbl_fphists:
        if row['snr'] < 3.0:
            row['origin'] = 'magul' # mag upper limit
    tbl_alert_fp = vstack([tbl_alerts,tbl_fphists])
    tbl_alert_fp.sort('jd')
    
    return tbl_alert_fp

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski to fetch \
transient light curves.')
    parser.add_argument('--n', dest='names', nargs='+', required=False,
                        help='Names of the ZTF candidates; if given, \
the coordinates will be queried from kowalski', default=None)
    parser.add_argument('--f', dest='filename', type=str, required=False,
                        help='Input CSV filename', default=None)
    parser.add_argument('--out', dest='out', type=str, required=False,
                        help='Output filename', default='lightcurves.csv')

    args = parser.parse_args()

    if args.names is None and args.filename is not None:
        a = ascii.read(args.filename, format='csv')
        args.names = list(a['name'])
    elif args.names is None and args.filename is None:
        print("No input candidates. Please use --n and provide ZTF names\
or --file to use a CSV file")
        #exit()

    # Read the secrets
    secrets = ascii.read('secrets.csv', format='csv')
    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]

    # Get the light curves
    light_curves_alerts = get_lightcurve_alerts(username_kowalski,
                                                password_kowalski,
                                                args.names)
    
    # Add prv_candidates photometry to the light curve
    light_curves_aux = get_lightcurve_alerts_aux(username_kowalski,
                                                 password_kowalski,
                                                 args.names)

    light_curves = light_curves_alerts + light_curves_aux

    # Create a Table and output CSV file
    create_tbl_lc(light_curves, outfile=args.out)
