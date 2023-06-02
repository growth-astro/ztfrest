#!/data/ia/anaconda3/bin/python
'''
Given an input skymap and 
a list of coordinates, find into which
probability contour they are enclosed

'''

import numpy as np
import json
import pdb

from astropy.time import Time
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
import requests

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
import healpy as hp

from penquins import Kowalski


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_programidx(program_name, username, password):
    ''' Given a marshal science program name, it returns its programidx'''

    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi', auth=(username, password))
    programs=json.loads(r.text)
    program_dict={p['name']:p['programidx'] for i,p in enumerate(programs)}

    try:
        return program_dict[program_name]
    except KeyError:
        print(f'The user {username} does not have access to the program {program_name}')
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
        try:
            coords=SkyCoord(ra=s['ra']*u.deg, dec=s['dec']*u.deg, frame='icrs')
            sources_out.append({"name":s['name'], "ra":coords.ra, "dec":coords.dec, \
	        "classification":s['classification'], "redshift":s['redshift'], "creation_date":s['creationdate'], \
		"id": s['id'], "candid": s['candid']})
        except ValueError:
            print("Error for source:",s)
    return sources_out

def query_kowalski_coords(username, password, names, jd_trigger=0):
    '''Query kowalski to get the coordinates of given ZTF sources. '''

    k = Kowalski(username=username, password=password, verbose=False)

    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts",
                   "filter": {"objectId": {"$in": names}
                              },
                    "projection": {"_id": 0, 
                        "candid": 1, 
                        "objectId": 1,
                        "candidate.jdstarthist": 1,
                        "candidate.ra": 1,
                        "candidate.dec": 1
                        },
     }
     }
    results_all = k.query(query=q)
    results = results_all['data']
    sources = []
    if results is None:
        return sources
    for n in names:
        source = {}
        source["jdstarthist"] = list(r["candidate"]["jdstarthist"] for r in results if r["objectId"] == n)[0]
        if source['jdstarthist'] <= jd_trigger:
            continue
        source["name"] = n
        source["ra"] = list(r["candidate"]["ra"] for r in results if r["objectId"] == n)[0]
        source["dec"] = list(r["candidate"]["dec"] for r in results if r["objectId"] == n)[0]
        sources.append(source)
        source["jdstarthist"] = list(r["candidate"]["jdstarthist"] for r in results if r["objectId"] == n)[0]
        print(Time(source["jdstarthist"], format='jd').iso)
    return sources


def query_kowalski_clu(username, password, clu):
    '''Query kowalski to get a table of CLU galaxies. '''
    #FIXME: change it to be a cone search

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "general_search", 
        "query": "db['CLU_20180513'].find({},{'distmpc': 1})" 
        }
    r = k.query(query=q)

    return r


def check_marshal_sources(sources_kowalski, marshal_sources, program_name):
    '''Check if the selected sources are present in a given
    science program.  If so, print out the relevant information.'''

    sources_in_program = []
    sources_not_in_program = []
    list_marshal_sources = list(s['name'] for s in marshal_sources)

    for source in sources_kowalski:
        print("-------")
        if source in list_marshal_sources:
            marshal_source = marshal_sources[np.where(np.array(list_marshal_sources) == source)[0][0]]
            try:
                for k in marshal_source.keys():
                    print(f"{k}: {marshal_source[k]}")
                sources_in_program.append(source)
            except:
                pdb.set_trace()
        else:
            print(f"{source} was not saved in CLU")
            sources_not_in_program.append(source)
        print("-------")
    print("Summary:")
    print(f"Sources saved in {program_name}: {sources_in_program}")
    print(f"Sources not saved in {program_name}: {sources_not_in_program}")

    return


def read_skymap(skymap_filename):
    '''Read the healpix skymap'''

    hpx = hp.read_map(skymap_filename, verbose = False)

    return hpx


def do_getfields(healpix, FOV=60/3600.0, ra=None, dec=None, radius=None, level=None, names=None):
    from ligo.skymap import postprocess
    import matplotlib

    ras, decs = list(float(r) for r in ra), list(float(d) for d in dec)

    if (not level is None):
        cls = 100 * postprocess.find_greedy_credible_levels(healpix)
        paths = postprocess.contour(cls, [level], degrees=True, simplify=True)
        paths = paths[0]

        pts = np.vstack((ras, decs)).T
        idx = np.zeros((len(ras)))
        for path in paths:
            polygon = matplotlib.path.Path(path)
            check = polygon.contains_points(pts)
            check = list(map(int, check))
            idx = np.maximum(idx, check)
        idx = np.where(idx == 1)[0]
        ras, decs = np.array(ras), np.array(decs)
        ras, decs = ras[idx], decs[idx]
    print(f"-> Sources included in the {level}% probability contour: {len(ras)}/{len(ra)}")
    if (not names is None) and (len(ras) > 0):
        names = np.array(names)
        names_out = names[idx]   
        print(f"-> Specifically: {list(x for x in names_out)}")
    else:
        names_out = []
    #for rr, dd, nn in zip(ras, decs, names_out):
    #    print(f"{nn}, {rr}, {dd}")
    print("number of candidates: ", len(names_out))
    return ras, decs, names_out


def get_source_autoannotations(sourceid, username, password):
    ''' Fetch a specific source's autoannotations from the GROWTH marshal. '''
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(username, password), data={'sourceid':str(sourceid)
}) 
    summary=json.loads(r.text)
    autoannotations=summary['autoannotations']
    autoannotations_string = ""
    for auto in autoannotations:
        autoannotations_string = f"{autoannotations_string}; \n {auto['username']}, {auto['type']}, {auto['comment']}" 

    return autoannotations_string


def add_autoannotations(name_source, program_name, username, password, event_name, level, skymap):
    ''' Check if the source is present on the GROWTH marshal and,
    if so, add autoannotations '''

    sources = get_candidates_growth_marshal(program_name, username, password)
    id_source=None
    for s, i in zip(sources, np.arange(len(sources))):
        if s['name'] == name_source:
            id_source = s['id']
            break
        else:
            continue

    autoannotations = get_source_autoannotations(s["id"], username, password)

    #Check that the new comment is not already present.  In that case, do not post.
    new_comment = f"{int(level)}% area, {event_name}, {skymap}"
    if new_comment in autoannotations:
        print('Autoannotation already present. Skipping..')

        return
    else:
        r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/add_autoannotation.cgi', auth=(username, password), \
            data={'id': -1, 'datatype': 'STRING','action': 'commit', 'sourceid':str(id_source), 
            'type': 'inside_area', 'comment': new_comment})

        return


def get_jdstarthist_kowalski(source_names, username, password):
    '''Query kowalski, look for one alert and return a list of candidate.jdstarthist '''
    k = Kowalski(username=username, password=password, verbose=False)   
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {"objectId": {"$in": list(source_names)}},
             "projection": {"_id": 0, "objectId": 1, "candidate.jdstarthist": 1}
            }
        }
    r = k.query(query=q)
    jdstarthist_tuples = set(list((s['objectId'], s['candidate']['jdstarthist']) for s in r['result_data']['query_result']))
    #re-output both lists to avoid mixed sorting
    names_out = list(x[0] for x in jdstarthist_tuples)
    jdstarthist_list = list(x[1] for x in jdstarthist_tuples)

    return names_out, jdstarthist_list


def get_hist_kowalski(source_names, username, password):
    '''Query kowalski, look for one alert and return a list of candidate.jdstarthist '''
    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {"objectId": {"$in": list(source_names)}},
             "projection": {"_id": 0, "objectId": 1, "candidate.jdstarthist": 1, "candidate.jdendhist": 1}
            }
        }


    r = k.query(query=q)
    jdstarthist_tuples = set(list((s['objectId'], s['candidate']['jdstarthist']) for s in r['result_data']['query_result']))
    jdendhist_tuples = set(list((s['objectId'], s['candidate']['jdendhist']) for s in r['result_data']['query_result']))
  
    #re-output both lists to avoid mixed sorting
    names_out = list(x[0] for x in jdstarthist_tuples)
    jdstarthist_list = list(x[1] for x in jdstarthist_tuples)
    jdendhist_list = []
    for n in names_out:
        max_jdendhist = 0.
        for tt in jdendhist_tuples:
            if tt[0] == n and tt[1] > max_jdendhist:
                max_jdendhist = tt[1]
        jdendhist_list.append(max_jdendhist)

    return names_out, jdstarthist_list, jdendhist_list


def select_sources_in_level2(sources, skymap_filename, level=90):
    hpx, header = hp.read_map(skymap_filename, h=True, verbose=False)
    i = np.flipud(np.argsort(hpx))
    sorted_credible_levels = np.cumsum(hpx[i])
    credible_levels = np.empty_like(sorted_credible_levels)
    credible_levels[i] = sorted_credible_levels
    npix = len(hpx)
    nside = hp.npix2nside(npix) 
    sources_within = list(s for s in sources
                          if (credible_levels[hp.ang2pix(
                                                        nside,
                                                        0.5 * np.pi - np.deg2rad(s["dec"]),
                                                        np.deg2rad(s["ra"]))
                                             ] <= level/100.
                                              ))
    return sources_within


def select_sources_in_level(sources, skymap_filename, level=90):
    """Select only those sources within a given contour
    level of the skymap"""

    skymap_prob = hp.read_map(skymap_filename, field=0, verbose=False)
    skyimap_prob = skymap_prob / np.sum(skymap_prob)
    sort_idx = np.argsort(skymap_prob)[::-1]
    csm = np.empty(len(skymap_prob))
    csm[sort_idx] = np.cumsum(skymap_prob[sort_idx])
    ipix_keep = np.where(csm <= level/100.)[0]
    nside = hp.npix2nside(len(skymap_prob))
    try:
        sources_within = list(s for s in sources
                              if (hp.ang2pix(
                                             nside,
                                             0.5 * np.pi -
                                             np.deg2rad(s["dec"]),
                                             np.deg2rad(s["ra"])
                                             ) in ipix_keep
                                      )
                              )
    except:
        sources_within = list(s for s in sources
                              if (hp.ang2pix(
                                             nside,
                                             0.5 * np.pi -
                                             np.deg2rad(s["dec"].value),
                                             np.deg2rad(s["ra"].value)
                                             ) in ipix_keep
                                      )
                              )

    return sources_within


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski.')
    parser.add_argument('--skymap', dest='skymap_filename', type=str, required=False, \
    help='Skymap filename', default = None)
    parser.add_argument('--level', dest='level', type=float, required=False, \
    help='Enclosed probability', default = 90)
    parser.add_argument('--program', dest='program', type=str, required=False, \
    help='All the sources present in the given marshal science program \
    will be checked.', default = None)
    parser.add_argument('--event', dest='event_name', type=str, required=False, \
    help='Name of the multi-messenger event (e.g., GW170817)', default = None)
    parser.add_argument('--annotate', dest='annotate', type=str2bool, required=False, \
    help='If True, autoannotations will be posted to the GROWTH marshal.', default = False)
    parser.add_argument('--fov', dest='fov', type=float, required=False, \
    help='Field of view of each cone (radius, in arcmin)', default = 60)
    parser.add_argument('--ra', dest='ra_center', nargs='+', required=False, \
    help='Right ascension of the center (array, in degrees)')
    parser.add_argument('--dec', dest='dec_center', nargs='+', required=False, \
    help='Declination of the center (array, in degrees)')
    parser.add_argument('--n', dest='names', nargs='+', required=False, \
    help='Names of the ZTF candidates; if given, the coordinates\
    will be queried from kowalski', default=None)
    parser.add_argument('--pix', dest='pix', type=str2bool, required=False, \
    help='Check by using the pixels instead of the contours', default=True)
    parser.add_argument('--radius', dest='radius', type=float, required=False, \
    help='Search radius (min)', default = 1.)
    parser.add_argument('--after-trigger', dest='after_trigger', type=str2bool, required=False, \
    help='Query only alerts whose first detection occurred after a certain date. \
    If this boolean value is True, then --jd-trigger must be given.', default = False)  
    parser.add_argument('--jd-trigger', dest='jd_trigger', type=float, required=False, \
    help='Julian Day of the trigger', default = -1.)            
    parser.add_argument('--min-days', dest='min_days', type=float, required=False, \
    help='Minimum time (days) between the first and last alert', default = 0.)
    parser.add_argument('--max-days', dest='max_days', type=float, required=False, \
    help='Maximum time (days) between the first and last alert', default = 30.)
    parser.add_argument('--out', dest='out', type=str, required=False, \
    help='Output filename', default = 'results.txt')
    
    args = parser.parse_args()

    #Check
    if args.annotate == True and args.event_name == None:
        print("-> ERROR: if you want autoannotation, please provide a name for the multi-messenger event using --event")
        exit()

    #Read the secrets
    secrets = ascii.read('secrets.csv', format = 'csv')

    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]

    #Check if the candidates are present in a list of science programs on the Marshal
    username_marshal = secrets['marshal_user'][0]
    password_marshal= secrets['marshal_pwd'][0]

    if (args.program is not None):
        #Get all the candidates from a certain science program on the marshal
        program_name = args.program
        marshal_sources = get_candidates_growth_marshal(program_name, username_marshal, password_marshal)    
        if marshal_sources is None:
            print("-> No saved candidates found.  Is the program name spelled correctly?")
            print("-> Exiting...")
            exit()
        ra_center = list(s['ra'].value for s in marshal_sources)
        dec_center = list(s['dec'].value for s in marshal_sources)


    #Read the skymap and check the level
    if args.skymap_filename is not None:
        healpix = read_skymap(args.skymap_filename)
        with fits.open(args.skymap_filename) as smp:
            skymap_header = smp[1].header

        if args.pix is True:
            if (args.program is not None):
                sources = marshal_sources
            elif (args.names is not None):
                sources = query_kowalski_coords(username_kowalski, password_kowalski, args.names, jd_trigger=args.jd_trigger)
            else:
                sources = []
                for r, d in zip(args.ra_center, args.dec_center):
                    sources.append({"ra":float(r), "dec":float(d)})
            sources_within = select_sources_in_level(sources, args.skymap_filename, args.level)
            string_out=""
            for source in sources:
                if source in sources_within:
                    print(f"{source['name']}, {source['ra']}, {source['dec']}")
                    string_out = f"{string_out} --radec {source['ra']} {source['dec']}"
                else:
                    #print(f"{source['name']}, {source['ra']}, {source['dec']}")
                    print("OUT", source)

            print(string_out)

        if (args.program is None):
            if (args.names is None):
                ra_center, dec_center, names_out = do_getfields(healpix, FOV=args.fov/60., ra=args.ra_center, dec=args.dec_center, level=args.level)
        else:
            names_list = list(s['name'] for s in marshal_sources)
            ra_center, dec_center, names_out = do_getfields(healpix, FOV = args.fov/60., ra=ra_center, dec=dec_center, level = args.level, names = names_list)

            '''Autoannotations '''
            if args.annotate: 
                skymap = args.skymap_filename.split('/')[-1]

                ''' Which sources are ZTF? '''	           
                names_ZTF = []
                names_not_ZTF = []
                for n in names_out:
                    if n[0:3] == 'ZTF':
                        names_ZTF.append(n)
                    else:
                        names_not_ZTF.append(n)
                '''ZTF only'''
                #Query kowalski to find out the jd of the history start.
                if names_ZTF != []:
                    print("ZTF names: ", len(names_ZTF))
                    names_out, jdstarthist_list, jdendhist_list = get_hist_kowalski(names_ZTF, username_kowalski, password_kowalski)
                    for name_source, jdstarthist, jdendhist in zip(names_out,jdstarthist_list, jdendhist_list):
                        #if (jdstarthist > Time(skymap_header['MJD-OBS'], format = 'mjd').jd) and (jdstarthist < Time(skymap_header['MJD-OBS'], format = 'mjd').jd + 30.) and (jdendhist - jdstarthist > args.min_days) and (jdendhist - jdstarthist < args.max_days):
                        add_autoannotations(name_source, program_name, username_marshal, password_marshal, args.event_name, args.level, skymap)
                        print(f"{name_source} passed")
                        #else:
                        #    print(f"{name_source} NOT passed")
                if names_not_ZTF != []:
                         for name_source in names_not_ZTF:
                             print(f"{name_source} passed (not ZTF)")
                             add_autoannotations(name_source, program_name, username_marshal, password_marshal, args.event_name, args.level, skymap)
    else:
        if args.program is not True:
            ra_center, dec_center = args.ra_center, args.dec_center 

    exit()


    program_names = ['Electromagnetic Counterparts to Gravitational Waves', 'Afterglows of Fermi Gamma Ray Bursts', 'Census of the Local Universe']    
    for program_name in program_names:
        marshal_sources = get_candidates_growth_marshal(program_name, username_marshal, password_marshal)    

        #For each transient check if it is present in the CLU science program
        check_marshal_sources(sources_kowalski, marshal_sources, program_name)

    print("Done.")
