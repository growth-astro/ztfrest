import psycopg2
from astropy.io import ascii
from functions_db import connect_database
import pandas as pd

con, cur = connect_database(update_database=True, path_secrets_db='db_access.csv', dbname='db_lens_admin')

#lc = pd.read_sql_query("select * from lightcurve where programid = 1 or programid = 2", con)
#lc.to_csv("lc_lens_alerts_20201121_p12.csv")

lc = pd.read_sql_query("select name, \
                        ls_ra, ls_dec, ls_sep_arcsec, ls_z_spec, \
                        ls_z_phot_median, ls_z_phot_std, ls_type, \
                        ls_z_phot_l95, ls_z_phot_u95, ls_fluxz, ls_sep_arcsec from crossmatch", con)
lc.to_csv("photoz_lens_alerts_20201121_p12.csv")
