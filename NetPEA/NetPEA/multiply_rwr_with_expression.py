"""
Return rwr where probability of each gene times by its expression value

"""

import argparse
import numpy as np
import pandas as pd
from datetime import datetime


def parse_parameter():

    parser = argparse.ArgumentParser(description = "Return rwr where probability of each gene times by its expression value")

    parser.add_argument("-r", "--rwr_path",
                        required = True,
                        help = "path to RWR-probability matrix with cell by gene")
    parser.add_argument("-e", "--exp_path",
                        required = True,
                        help = "path to expression file with cell by gene (log2-transformed)")
    parser.add_argument("-o", "--output_path",
                        required = True,
                        help = "path to output file")
    parser.add_argument("-debug", "--DEBUG",
                        default = False,
                        type = bool,
                        help = "display print message if True")
    return parser.parse_args()

if __name__ == "__main__":
    # timer
    datetimeFormat = '%Y-%m-%d %H:%M:%S.%f'
    start_time = datetime.now()

    # get args
    args = parse_parameter()

    # load data
    print('loadding data.....')
    rwr_df = pd.read_csv(args.rwr_path, header=0, index_col=0, sep="\t")
    exp_df = pd.read_csv(args.exp_path, header=0, index_col=0, sep="\t")

    # find overlapping cells 
    overlap_cList = sorted(list( set(rwr_df.index)&set(exp_df.index) ))
    rwr_df = rwr_df.loc[overlap_cList]
    exp_df = exp_df.loc[overlap_cList]
    if args.DEBUG == True:
        print('    overlap cells={:}'.format(len(overlap_cList)))

    # find overlapping genes
    overlap_gList = sorted(list( set(rwr_df.columns)&set(exp_df.columns) ))
    if args.DEBUG == True:
        print('    overlap genes={:}'.format(len(overlap_gList)))

    # subsetting to include overlap cell
    print('subsetting to include overlap cells.....')
    rwr_df = rwr_df.loc[overlap_cList]
    exp_df = exp_df.loc[overlap_cList]
    if args.DEBUG == True:
        print('    rwr={:}, exp={:}'.format(rwr_df.shape, exp_df.shape))

    # split rwr by columns
    print('split rwr into columns with/without ppi genes.....')
    rwr_has_exp_col_df = rwr_df[overlap_gList] 
    rwr_has_no_exp_col_df = rwr_df[rwr_df.columns[~rwr_df.columns.isin(overlap_gList)]]
        
    print('    #expressed genes in rwr={:}/{:}'.format(rwr_has_exp_col_df.shape[1], rwr_df.shape[1]))

    # multiply rwr with exp
    print('multiply rwr with exp for columns with ppi genes.....')
    new_rwr_df = rwr_df[overlap_gList]*exp_df[overlap_gList]
    if args.DEBUG == True:
        print(rwr_df[overlap_gList].iloc[:3,:3])
        print(exp_df[overlap_gList].iloc[:3,:3])
        print(new_rwr_df.iloc[:3,:3])

    # merge back to rwr_df
    print('merging files.....')
    df = pd.concat([new_rwr_df, rwr_has_no_exp_col_df], axis=1)
    print('    df={:}'.format(df.shape))
    
    # saving to file
    df.to_csv(args.output_path+'.RWR_timesEXP.txt', header=True, index=True, sep="\t")
   
    # stop timer
    spend = datetime.strptime(str(datetime.now()), datetimeFormat) - datetime.strptime(str(start_time),datetimeFormat)
    print('[Finished in {:}]'.format(spend))
