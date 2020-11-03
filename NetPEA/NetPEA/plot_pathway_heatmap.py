"""
Return heatmap of enrichment score/P-value/FDR for pathways


Note: for ssGSEA, enrichment score is NES, while NetPEA would be zscore
"""

import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")

def heatmap(df, fignameStr, barlblStr):
    """
    :param df: dataframe with index at the 1st column, header at the 1st row
    :param fignameStr: string representing figure name
    :param barlblStr: string representing label string
    :return fig:
    """
    # colormap
    cmapDict = {'ygb': "YlGnBu"}
    # figure
    fig, ax = plt.subplots(1,1, figsize=(12,12))
    figname = fignameStr + '.heatmap.png'
    # plotting
    ax = sns.heatmap(df, cmap=cmapDict['ygb'], cbar_kws={'label': barlblStr})
    # save to file
    fig.savefig(figname, dpi=200, bbox_inches='tight')
    # return
    print( 'plotting heatmap, figure={:}'.format(figname) )
    return fig

def parse_parameter():
    parser = argparse.ArgumentParser(description = "Plot heatmap")
    parser.add_argument("-s", "--score_path",
                        required = True,
                        help = "path to enrichment score file (cell by pathways)")
    parser.add_argument("-p", "--permutation_path",
                        required = True,
                        help = "path to permutation result folder")
    parser.add_argument("-m", "--type_str", 
                        choices = ['ssGSEA', 'NetPEA'],
                        required = True,
                        help = "string representing type of erichment analysis")
    parser.add_argument("-o", "--output_path",
                        required = True,
                        help = "path to output files")
    return parser.parse_args()

if __name__ == "__main__":
    # get args
    args = parse_parameter()

    # load data
    score_df = pd.read_csv(args.score_path, header=0, index_col=0, sep="\t") # cell (or drug) by pathways

    # collect per-cell/per-drug data
    pval_df_list =[]
    fdr_df_list = []
    zscore_df_list = []
    if args.type_str == 'ssGSEA':
        for idx in score_df.index:
            # get data
            f = args.permutation_path + '/' + idx + '/gseapy.ssgsea.gene_sets.report.txt'
            idx_df = pd.read_csv(f, skiprows=2, header=0, index_col=0, sep="\t") 
            pval_df = idx_df[['pval']]
            fdr_df = idx_df[['fdr']]
            # replace column name
            pval_df.columns = [idx]
            fdr_df.columns = [idx]
            # append to result list
            pval_df_list.append(pval_df)
            fdr_df_list.append(fdr_df)
        # merge
        all_pval_df = pd.concat(pval_df_list, axis=1)
        all_fdr_df = pd.concat(fdr_df_list, axis=1)
    elif args.type_str == 'NetPEA':
        for idx in score_df.index:
            # get data
            f = args.permutation_path + '/' +  args.output_path + '.' + idx + '.NetPEA.background_result.txt'
            idx_df = pd.read_csv(f,  header=0, index_col=0, sep="\t")
            zscore_df = idx_df[['zscore']]
            pval_df = idx_df[['pvalue']]
            # replace column name
            pval_df.columns = [idx]
            zscore_df.columns = [idx]
            # append to result list
            pval_df_list.append(pval_df)
            zscore_df_list.append(zscore_df)
        # merge
        all_pval_df = pd.concat(pval_df_list, axis=1)
        all_zscore_df = pd.concat(zscore_df_list, axis=1)

    else:
        print('ERROR')
    
    # sort pathways
    if args.type_str == 'ssGSEA':
        all_score_df = score_df.T.sort_index() # pathways by cell/drug
        all_pval_df = all_pval_df.sort_index() # pathways by cell/drug
        all_fdr_df = all_fdr_df.sort_index() # pathways by cell/drug    
        # make heatmaps
        heatmap(score_df.T, args.output_path+'.NES', 'Normalized Enrichment Score (NES)')
        heatmap(all_pval_df, args.output_path+'.P-value', 'P-value')
        heatmap(all_fdr_df, args.output_path+'.FDR', 'FDR')
        # save to file
        all_pval_df.to_csv(args.output_path+'.P-value.txt', header=True, index=True, sep="\t")
        all_fdr_df.to_csv(args.output_path+'.FDR.txt', header=True, index=True, sep="\t")
    elif args.type_str == 'NetPEA':
        all_score_df = score_df.T.sort_index() # pathways by cell/drug
        all_pval_df = all_pval_df.sort_index() # pathways by cell/drug
        # make heatmaps
        heatmap(score_df.T, args.output_path+'.Zscore', 'Z-score')
        heatmap(all_pval_df, args.output_path+'.P-value', 'P-value')
        # save to file
        all_pval_df.to_csv(args.output_path+'.P-value.txt', header=True, index=True, sep="\t")

