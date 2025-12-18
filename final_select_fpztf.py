'''
performs final selection pass on ForcePhotZTF data.
Necessary because ForcePhotZTF only works in old python env
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
    parser.add_argument('--targetdir-base', dest='targetdir_base', type=str,
                        required=False,
                        help='Directory for the forced photometry',
                        default='./forced_photometry/')
    args = parser.parse_args()
    
    
    t_for_phot = ascii.read('needs_fp.csv')

    if len(t_for_phot) > 0:
        if args.dontReadDb:
            read_database = False
        else:
            read_database = True
        from select_variability_db import select_variability
        # Repeat the selection based on forced photometry from ForcePhotZTF
        print("----------------------------------------------")
        print("Selecting based on ForcePhotZTF photometry...")
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
            
        # Repeat the selection based on stacked forced photometry (ForcePhotZTF from psql db)
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
        
    else:
       print("No FPZTF. Skipping selection steps. \n ...Done") 
