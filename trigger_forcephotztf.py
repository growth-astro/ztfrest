'''
Given a csv file of candidates needing ForcePhotZTF photometry:
trigger ForcePhotZTF (in old python environment) and write
results to database.
'''

if __name__ == '__main__':
    from astropy.io import ascii
    import argparse
    
    parser = argparse.ArgumentParser(description='trigger fpztf')
    # parser.add_argument("--doForcePhot", action="store_true", help="Trigger forced photometry using ForcePhotZTF",
                        # default=False)
    parser.add_argument("--doWriteDb",  action='store_true',
                        help='Write information to the psql database \
                        (needs admin privileges)',
                        default=False)
    parser.add_argument("--dontReadDb", action="store_true", help="Do not read information \
                        from the psql database", default = False)
    parser.add_argument('--path-secrets-db', dest='path_secrets_db', type=str,
                        required=False,
                        help="Path to the CSV file including the credentials \
                        to access the psql database", default='db_access.csv')
    parser.add_argument('--path-input-table', dest='path_input_table', type=str,
                        required=False,
                        help="Path to the CSV file including the table of cnadidates needing ForcePhotZTF",
                        default='/data/ia/ztfrest_sarah/needs_fp.csv')
    parser.add_argument('--targetdir-base', dest='targetdir_base', type=str,
                        required=False,
                        help='Directory for the forced photometry',
                        default='./forced_photometry/')
    args = parser.parse_args()
    
    t_for_phot = ascii.read(args.path_input_table)

    if len(t_for_phot) > 0:
        print("----------------------------------------")
        print("Triggering forced photometry pipeline...")
        from forcephot import trigger_forced_photometry
        # Trigger forced photometry
        success, _ = trigger_forced_photometry(t_for_phot,
                                               args.targetdir_base,
                                               daydelta_before=7.,
                                               daydelta_after=14.)

        if args.doWriteDb and len(success) > 0:
            from functions_db import connect_database

            con, cur = connect_database(path_secrets_db='/data/ia/ztfrest_sarah/db_access.csv', dbname='db_kn_2025_admin')
            print("----------------------------------------")
            print("Updating database...")
            cur.execute("select current_database()")
            current_dbname = cur.fetchall()[0][0]
            if current_dbname != 'database_kn_2025':
                print("Wrong database!!! terminating connection...")
                con.close()
                cur.close()
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
            
    else:
        print("There are no candidates needing ForcePhotZTF.")
    
    # if len(t_for_phot) > 0:
        # if args.dontReadDb:
            # read_database = False
        # else:
            # read_database = True
        # from select_variability_db import select_variability
        # # Repeat the selection based on forced photometry from ForcePhotZTF
        # print("----------------------------------------------")
        # print("Repeating based on ForcePhotZTF photometry...")
        # selected, rejected, cantsay = select_variability(t_for_phot,
                   # hard_reject=[], update_database=args.doWriteDb,
                   # read_database=read_database,
                   # use_forcephotztf=True, stacked=False,
                   # baseline=0.125, var_baseline={'g': 6, 'r': 8, 'i': 10},
                   # max_duration_tot=15., max_days_g=7., snr=4,
                   # index_rise=-1.0, index_decay=0.3,
                   # path_secrets_db=args.path_secrets_db,
                   # save_plot=True, path_plot='./plots/',
                   # show_plot=False, use_metadata=False,
                   # path_secrets_meta='./secrets.csv',
                   # save_csv=True, path_csv='./lc_csv',
                   # path_forced='./')
            
        # # Repeat the selection based on stacked forced photometry (ForcePhotZTF from psql db)
        # print("-----------------------------------------------")
        # print("Repeating based on stacked forced photometry...")
        # selected, rejected, cantsay = select_variability(t_for_phot,
                   # hard_reject=[], update_database=args.doWriteDb,
                   # read_database=read_database,
                   # use_forcephotztf=True, stacked=True,
                   # baseline=0.125, var_baseline={'g': 6, 'r': 8, 'i': 10},
                   # max_duration_tot=15., max_days_g=7., snr=4,
                   # index_rise=-1.0, index_decay=0.3,
                   # path_secrets_db=args.path_secrets_db,
                   # save_plot=True, path_plot='./plots/',
                   # show_plot=False, use_metadata=False,
                   # path_secrets_meta='./secrets.csv',
                   # save_csv=True, path_csv='./lc_csv',
                   # path_forced='./')
                   
    # print("...Done.")
