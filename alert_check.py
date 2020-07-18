
import os
import requests

import numpy as np
from astropy.time import Time
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
import healpy as hp
from astropy.io import fits

from penquins import Kowalski
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import astropy.stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import io

import optparse
__version__ = 0.1

text = ['No alerts around, safe','Less than 3 alerts around','More than 3 alerts around, needs checking','No data retrievbed, try again, maybe change rb to drb']

collection_ZTF_alerts = 'ZTF_alerts'
collection_ZTF_alerts_aux = 'ZTF_alerts_aux'


def get_radec(kow, n):
    q = {"query_type": "find",
         "query": {
             "catalog": collection_ZTF_alerts,
             "filter": {"objectId": n},
             "projection": {"_id": 0, "cutoutScience": 0, "cutoutTemplate": 0, "cutoutDifference": 0},
         }
         }
    r1 = kow.query(query=q)
    alerts1 = r1['result_data']['query_result']
    radec = (alerts1[0]['candidate']['ra'],alerts1[0]['candidate']['dec'])
    return (radec[0],radec[1]),alerts1[0]['candidate']['jdstarthist'],alerts1[0]['candidate']['ndethist']

def get_alerts_around(kow,radec,r,drb=0.25):
    q = {"query_type": "cone_search",
         "object_coordinates": {
             "radec": '[({}, {})]'.format(radec[0],radec[1]), 
             "cone_search_radius": str(r),
             "cone_search_unit": "arcsec"
         },
         "catalogs": {
             "ZTF_alerts": {
                 "filter": {"candidate.rb": {"$gt": drb}}, #try to modify to drb
#                  "filter": {"candidate.drb": {"$gt": 0.25}},
                     "projection": {
                     "objectId": 1,
                     "candidate.ra": 1,
                     "candidate.dec": 1,
                     "candidate.rcid": 1,
                     "candidate.drb": 1,
                     "_id": 0
                 }
             }

         },
         "kwargs": {
             "limit": 30
         }
         }
    return kow.query(query=q)
    

def check_alerts(kow, radec,r=25,n='ztfcand', names_plot = False):
    if radec == None:
        radec,_,_ = get_radec(kow, n)
#         print(radec,n)
    r = get_alerts_around(kow, radec,r)
    data = r['result_data']['ZTF_alerts']
    coord_name = list(data.keys())[0]
    
    # getting alerts names
    id_names=[]
    if len(data[coord_name])>0:
#         print('coordinates:',radec[0],radec[1])
        for i in range(len(data[coord_name])):
#             print(data[coord_name][i]['objectId'])
            id_names.append(data[coord_name][i]['objectId'])
    else:
        print(n,'no alerts, probably bad name')
        print(radec,len(n))
        return plt.figure(figsize=(15,10)),[],[],[],[],[],[],[],[]
        
    id_names = np.unique(id_names)
    jd_start,n_det=[],[]
    ra,dec=[],[]
    for name in id_names:
        radec_,jdstart_,ndet_ = get_radec(kow, name)
        jd_start.append(jdstart_)
        n_det.append(ndet_)
        ra.append(radec_[0])
        dec.append(radec_[1])
    jd_start=np.asarray(jd_start)
    ra=np.asarray(ra)[np.argsort(jd_start)]
    dec=np.asarray(dec)[np.argsort(jd_start)]
    id_names = id_names[np.argsort(jd_start)]
    jd_start = jd_start[np.argsort(jd_start)]
    n_det = np.asarray(n_det)[np.argsort(jd_start)]
    cords = SkyCoord(ra, dec, unit="deg")
    dists = np.array([cords[0+i].separation(cords[0+i+1]).to(u.arcsec).value for i in range(len(cords)-1)])
    dist_from_alert = SkyCoord(radec[0], radec[1], unit="deg").separation(cords).to(u.arcsec).value
#     print(dist_from_alert)
    dt = np.array([jd_start[0+i] - jd_start[0+i+1] for i in range(len(jd_start)-1)])
    v = dists/dt 
    
    ra_temp,dec_temp=[],[]
    for i in range(len(cords)):
        ra_temp.append((cords[i].ra.arcsec))
        dec_temp.append((cords[i].dec.arcsec))

    f = plt.figure(figsize=(15,10))
    ax1 = plt.subplot(221)
    for i in range(len(cords)):
        ttt = str(i+1)
        if names_plot == True:
            ttt = str(i)+' '+id_names[i][7:]
        ax1.text((cords[i].ra.arcsec-(min(ra_temp)-2)),(cords[i].dec.arcsec-(min(dec_temp)-2)),ttt)

        if id_names[i] == n:
            ax1.scatter((cords[i].ra.arcsec-(min(ra_temp)-2)),(cords[i].dec.arcsec-(min(dec_temp)-2)), c='red')

    ra_max,ra_min = (max(ra_temp)+2),(min(ra_temp)-2)
    dec_max,dec_min = (max(dec_temp)+2),(min(dec_temp)-2)
    d_ra,d_dec = ra_max-ra_min,dec_max-dec_min
    ax1.set_xlim(0,max(d_ra,d_dec))
    ax1.set_ylim(0,max(d_ra,d_dec))
    ax1.set_xlabel('arcsec')
    ax1.set_ylabel('arcsec')
    ax1.set_title(n)

    ax2 = plt.subplot(222)
    bin_10 =  np.array([np.sum(dist_from_alert<=(10*i)) for i in range(1,11)])
    ax2.scatter(np.arange(1,11)*10,bin_10)
    ax2.set_xlabel('arcsec from candidate')
    ax2.set_ylabel('#')
    ax2.set_ylim(0,max(10,len(np.unique(id_names))))
    ax2.set_title('Cumulative # in Distance')

    axins3 = inset_axes(ax2,width=1.3, height=0.9)
    bin_1 =  np.array([np.sum(dist_from_alert<=i) for i in range(1,20)])
    axins3.scatter(np.arange(1,20),bin_1)
    axins3.set_xlabel('arcsec')
    axins3.set_ylim(0,max(5,np.amax(bin_1)))
    
    return f,id_names,jd_start,dists,dt,v,cords,dist_from_alert,n_det

def get_radec_jd(kow, n):
    
    q = {"query_type": "find",
     "query": {
         "catalog": collection_ZTF_alerts,
         "filter": {"objectId": n},
         "projection": {"_id": 0, "cutoutScience": 0, "cutoutTemplate": 0, "cutoutDifference": 0},
     }
     }
    r1 = kow.query(query=q)
    alerts1 = r1['result_data']['query_result']
    ra,dec,jd = [],[],[]
    for al in alerts1:
        ra_t,dec_t,jd_t = al['candidate']['ra'],al['candidate']['dec'],al['candidate']['jd']
        ra.append(ra_t)
        dec.append(dec_t)
        jd.append(jd_t)
    
    return np.array(ra),np.array(dec),np.array(jd)

def plot_science(kow,ax1,name,s_min,s_max):
    q = {"query_type": "general_search", "query": "db['ZTF_alerts'].find({'objectId': {'$eq': '"+name+"'}})" }
    r = kow.query(query=q,timeout=30)


    candidate = r['result_data']['query_result'][0]

    sdir = './temp/'

    if not os.path.isdir(sdir):
        os.makedirs(sdir)

    f=open(sdir+candidate['objectId']+'-sci.fits','wb').write(io.BytesIO(candidate['cutoutScience']['stampData']).getvalue()) 		#;f.close()
    
    cm = plt.cm.cubehelix
    stamp_ext = 'fits'

    ax1.margins(0.05)           # Default margin is 0.05, value 0 means fit
    itype='-sci.'
    img_data = fits.getdata(sdir+candidate['objectId']+itype+stamp_ext)
    img_data[np.isnan(img_data)]=np.nanmedian(img_data)
    # ax1.imshow((img_data), cmap=cm)
    ax1.set_xticks([])
    ax1.set_yticks([])

    imgd_scaled = np.log10(img_data)
    vmax = np.sort(imgd_scaled[28:36,28:36].flatten())[-3]
    npixel = (len(img_data)+1)**2
    imgd_flat = img_data.flatten()
    imgd_scaled[imgd_scaled<0]=np.nanmedian(img_data)
    v_onesig = np.log10(np.nanmedian(img_data) - astropy.stats.mad_std(imgd_flat[np.isnan(imgd_flat)==False])*1.4826)
    vmin= max(v_onesig, np.nanmin(imgd_scaled))
    ax1.imshow(imgd_scaled, cmap=cm, vmax=s_max[0]*vmax, vmin=vmin*s_min[0])

    ax1.axhline(y=31,xmin=0.5-0.15,xmax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
    ax1.axhline(y=31,xmin=0.5+0.05,xmax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)
    ax1.axvline(x=31,ymin=0.5-0.15,ymax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
    ax1.axvline(x=31,ymin=0.5+0.05,ymax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)

def alert_check_complete(kow, nn, plots = True,work_output = './output/'):
    f,id_names,jd_start,dists,times,v,cords,dist_from_alert,n_det = check_alerts(kow, None,r=15,n=nn)
    if len(id_names) == 0:
#         nodata.append(nn)
        index_check = 3
        print(nn,text[index_check])
        return index_check

    if plots:
        if not os.path.isdir(work_output):
            os.makedirs(work_output)

        ax3 = f.add_subplot(2,3,4)
        ax4 = f.add_subplot(2,3,5)
        ax5 = f.add_subplot(2,3,6)
    
        names_to_plot = id_names[np.argsort(dist_from_alert)]
        ndet_to_plot = n_det[np.argsort(dist_from_alert)]
        dist_to_plot = np.sort(dist_from_alert)

        plot_science(kow,ax4,nn,[1.001,1.001],[1.001,1.001])
        if len(ndet_to_plot) > 1:
            plot_science(kow,ax3,names_to_plot[1],[1.001,1.001],[1.001,1.001])
            ax3.set_title(str(int(dist_to_plot[1]))+' arcsec, ' + str(ndet_to_plot[1]) +' det')
        if len(ndet_to_plot) > 2: 
            plot_science(kow,ax5,names_to_plot[2],[1.001,1.001],[1.001,1.001])
            ax5.set_title(str(int(dist_to_plot[2]))+' arcsec, ' + str(ndet_to_plot[2]) +' det')
    
        ax4.set_title(nn,fontsize=14,loc= 'left',fontweight='bold')
        plt.savefig(work_output+nn+'.pdf')
    
    if np.sum(dist_from_alert < 20) > 3:
#         weird.append(nn)
        index_check = 2
        print(nn,text[index_check])
        return index_check
    elif np.sum(dist_from_alert < 20) == 1:
#         super_safe.append(nn)
        index_check = 0
        print(nn,text[index_check])
        return index_check
    else:
#         safe.append(nn)
        index_check = 1
        print(nn,text[index_check])
        return index_check
    
def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-n", "--name", help="Candidate name.", default='ZTF19acbjvwh')
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")
    parser.add_option("--doPlots",  action="store_true", default=False)
    parser.add_option("-o", "--outputDir", help="output file",default="./output/")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running alert_checking..."
        print >> sys.stderr, "version: %s"%__version__
        print >> sys.stderr, ""
        print >> sys.stderr, "***************** PARAMETERS ********************"
        for o in opts.__dict__.items():
            print >> sys.stderr, o[0]+":"
            print >> sys.stderr, o[1]
        print >> sys.stderr, ""

    return opts
    

if __name__ == "__main__":
    
    opts = parse_commandline()

    rank  = alert_check_complete(opts.name,plots = opts.doPlots,work_output = opts.outputDir)

