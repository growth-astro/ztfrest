# Author: Igor Andreoni

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import ascii
import psycopg2

from functions_db import connect_database

def get_lc(name, con, cur, forced=True, stack=False,
           writecsv=True, program_ids=[1,2,3], out_path='./',
           date_start=None, date_end=None):

    if forced is False:
        plot_alerts = True
        lc = pd.DataFrame(columns=['jd', 'mag', 'mag_unc', 'filter',
                                   'limmag', 'forced', 'programid'])
        # FIXME do retrieve the alert lc!
    else:
        if stack is True:
            table = 'lightcurve_stacked'
        elif forced is True:
            table = 'lightcurve_forced'
        if name is not None:
            wherename = f"name = '{name}'"
        else:
            wherename = f"name IS NOT NULL"
        if date_start is not None:
            # FIXME This takes jds with upper limits into account!!
            wherestart = f"jd > {Time(date_start).jd-30}"
        else:
            wherestart = f"jd > 0"
        if date_end is not None:
            # FIXME This takes jds with upper limits into account!!
            whereend = f"jd < {Time(date_end).jd+30}"
        else:
            whereend = f"jd < 100000"
        lc = pd.read_sql_query(f"SELECT name, jd, mag, mag_unc, filter, limmag, programid FROM {table} WHERE {wherename} AND {wherestart} AND {whereend}", con)
        lc["forced"] = np.ones(len(lc))

        if date_start is not None or date_end is not None:
            lc_selected = pd.DataFrame(columns=("jd", "mag", "mag_unc", "filter", "limmag", "programid", "forced"))
            for n in set(lc['name']):
                firstdet = np.min(lc["jd"][(lc['mag'] < 30) & (lc['mag']>0)])
                if firstdet < Time(date_start.jd):
                    continue
                else:
                    lastdet = np.max(lc["jd"][(lc['mag'] < 30) & (lc['mag']>0)])
                    if lastdet > Time(date_end.jd):
                        continue
                    else:
                        lc_selected.append(lc[lc['name'] == n], ignore_index=True)
            lc = lc_selected

    if writecsv is True:
        if forced is True:
            forcedbool = 1
        else:
            forcedbool = 0
        if stack is True:
            stackbool = 1
        else:
            stackbool = 0
        if writecsv is True:
            lc.to_csv(f"{out_path}/lc_{name}_forced{forcedbool}_stacked{stackbool}.csv")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query db')
    parser.add_argument('-n', dest='name', nargs='+', type=str,
                        help='ZTF name(s)')
    parser.add_argument('-f', action='store_true',
                        required=False,
                        help='Forced photometry')
    parser.add_argument('-s', action='store_true',
                        required=False,
                        help='Stacked forced photometry')
    parser.add_argument('-o', type=str, dest='out_path',
                        required=False,
                        help='Output directory path',
                        default='./')
    parser.add_argument('-date-start', dest='date_start', type=str,
                        help='Start date, e.g. 2020-12-01',
                        default=None)
    parser.add_argument('-date-end', dest='date_end', type=str,
                        help='End date, e.g. 2020-12-05',
                        default=None)

    args = parser.parse_args()

    # Connect to the database
    con, cur = connect_database(update_database=True,
                                path_secrets_db='db_access.csv')

    if len(args.name) > 0:
        for name in args.name:
            get_lc(name, con, cur, forced=args.f, stack=args.s,
                   writecsv=True, program_ids=[1,2,3],
                   out_path=args.out_path,
                   date_start=args.date_start, date_end=args.date_end)
    else:
        get_lc(None, con, cur, forced=args.f, stack=args.s,
               writecsv=True, program_ids=[1,2,3],
               out_path=args.out_path,
               date_start=args.date_start, date_end=args.date_end)
