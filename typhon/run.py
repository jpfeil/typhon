#!/usr/bin/env python2.7

import sys
sys.path.append('/opt/hydra/')

import argparse
import numpy as np
import pandas as pd

import library.analysis as hydra


ens_to_hugo = {}
with open('/opt/typhon/data/EnsGeneID_Hugo_Observed_Conversions.txt') as f:
    header = next(f)
    for line in f:
        hugo, ens = line.strip().split('\t')
        ens_to_hugo[ens] = hugo


def main():
    """
    Typhon Pipeline
    """
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--diagnosis',
                        help='Patient diagnosis',
                        required=True)

    parser.add_argument('--RSEM',
                        help='Path to N-of-1 RSEM file.',
                        required=True)

    args = parser.parse_args()

    available_models = ['MYCN-NA Neuroblastoma',
                        'Osteosarcoma',
                        'Ewing Sarcoma',
                        'Synovial Sarcoma']

    if args.diagnosis not in available_models:
        msg = "Please select one of the following diagnoses:\n%s" % '\n'.join(available_models)
        raise ValueError(msg)

    data = pd.read_csv(args.RSEM, sep='\t')

    data['hugo'] = data['gene_id'].map(ens_to_hugo)

    tpm = data.reindex(['hugo', 'TPM'], axis=1).groupby('hugo').sum()

    exp = np.log2(tpm + 1)

if __name__ == '__main__':
    main()

