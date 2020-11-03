"""
Return distribution plot 
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joypy as jpy
import seaborn as sns
sns.set(style="whitegrid")

def parse_parameter():
    parser = argparse.ArgumentParser(description = "Plot distribution")
    parser.add_argument("-i", "--input_path",
                        required = True,
                        help = "path to input file with feature and label")
    parser.add_argument("-o", "--output_path",
                        required = True,
                        help = "path to output files")
    return parser.parse_args()

if __name__ == "__main__":
    # get args
    args = parse_parameter()

    # load data 
    df = pd.read_csv(args.input_path, header=0, index_col=0, sep="\t") # cell (or drug) by pathways
    print(df.head())
    print(df.max(), df.min())
    # plot joyplot
    fig, ax = jpy.joyplot(df, figsize=(10,25), overlap=1)

    # save to file
    fig.savefig(args.output_path+'.joyplot.png', dpi=200, bbox_inches='tight')
    
