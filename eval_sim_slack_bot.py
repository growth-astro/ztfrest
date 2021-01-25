# Import the relevant packages
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u

from slack import RTMClient, WebClient
import numpy as np
import logging
import matplotlib.pyplot as plt
import io
import os
import sys
from astropy.time import Time
import traceback
import time

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default="sim")
    parser.add_argument("-c", "--channel", type=str, default="partnership")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-oc", "--onlycaltech", action="store_true", default=False)
    parser.add_argument("-np", "--noplots", action="store_true", default=False)
    parser.add_argument("-s", "--start_day", type=str, default="2020-9-20")
    parser.add_argument("-e", "--end_day", type=str, default="2021-01-20")

    cfg = parser.parse_args()

    outdir = os.path.join(cfg.outdir, cfg.channel)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if cfg.channel == 'partnership':
        channel = 'C01B3ME3GEQ'
        program_ids = [1,2]
    elif cfg.channel == 'caltech':
        channel = 'C01AUCVUKTP'
        program_ids = [1,2,3]
    elif cfg.channel == 'test':
        channel = 'G01A2AUV8Q2'
        program_ids = [1,2,3]
    else:
        print('Sorry, I do not know that channel...')
        exit(0)

    filenames = sorted(glob.glob(os.path.join(outdir,'*.dat')))
    transients = {}
    transients_unique = {}
    transients_all = []
    alllen, uniquelen = [], []
    for filename in filenames:
        key = filename.split("/")[-1].split(".")[0]
        with open(filename) as f:
            content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        content = [x.strip() for x in content] 
        transients[key] = content
        transients_unique[key] = list(set(content) - set(transients_all))

        transients_all = transients_all + content
        transients_all = list(set(transients_all))
 
        alllen.append(len(transients[key]))
        uniquelen.append(len(transients_unique[key]))

    print(len(alllen))
    print(np.max(alllen), np.mean(alllen), np.median(alllen), np.min(alllen), np.std(alllen))
    print(np.max(uniquelen), np.mean(uniquelen), np.median(uniquelen), np.min(uniquelen), np.std(uniquelen))

