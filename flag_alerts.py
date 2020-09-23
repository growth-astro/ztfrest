# Author: Igor Andreoni
# Wrapper for Tomas Ahumada's alert_check function

import argparse

import numpy as np
from astropy.io import ascii

from penquins import Kowalski
from alert_check import alert_check_complete


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check neighboring alerts')
    parser.add_argument('names', metavar='ZTF objectIds', type=str, nargs='+',
                        help='ZTF names of the candidates to check')
    args = parser.parse_args()

    # Read the secrets
    secrets = ascii.read('./secrets.csv', format='csv')
    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    kow = Kowalski(username=username, password=password)

    if len(args.names) > 0:
        print("Checking alerts...")
        ind_check_alerts = []
        for objid in args.names:
            index_check = alert_check_complete(kow, objid)
            ind_check_alerts.append(index_check)
        ind_check_alerts = np.array(ind_check_alerts)
        allids = np.asarray(args.names)[ind_check_alerts < 2]
        print("Safe:", allids)
        print("Flagged:", [a for a in args.names if not (a in allids)])
