"""
This is an implemntation of NetPEA proposed by Liu & Ruan (2013).
The network-based pathway enrichment analysis consists of the following steps:
1. perform random walk with restart (i.e., input gene sets) on a protein-protein interaction network
2. perform permutation test to estimate the significance of enrichment


Reference:
Liu, L., & Ruan, J. (2013). Network-based pathway enrichment analysis. In 2013 IEEE International Conference on Bioinformatics and Biomedicine (pp. 218–221).
"""


# import built-in pkgs
import sys
import argparse
import pandas as pd
from datetime import datetime

# import main module
import RWR as rwr
import NetPEA as pea

def times_expression(rwr, exp):
    """
    :param rwrDf: dataframe of cell by gene probability matrix
    :param expDf: dataframe of cell by gene expression matrix
    :return rwr_timesexp_df: dataframe of cell by gene probability matrix,
                             in which genes are multiplied with expression values

    Note: this function assumes cells are all overlapped while gene maybe not
    """
    cell_list = sorted(list(set(rwr.index) & set(exp.index)))
    gene_list = sorted(list(set(rwr.columns)&set(exp.columns)))

    if len(cell_list) == 0:
        print('ERROR! no overlapping cell lines')
        sys.exit(1)
    if len(gene_list) == 0:
        print('ERROR! no overlapping genes')
        sys.exit(1)

    # multiply with gene expression for overlapping cell, gene
    rwr_timesexp =  rwr.loc[cell_list, gene_list]*exp.loc[cell_list, gene_list]

    # concat with other gene
    out_gene_list = list(set(rwr.columns)-set(gene_list))
    out_df = pd.concat([rwr_timesexp, rwr[out_gene_list]], axis=1)
    return out_df

def parse_parameter():
    parser = argparse.ArgumentParser(description = "Return NetPEA zscore")

    #parser.add_argument("-r", "--rwr_path",
    #                    required = True,
    #                    help = "path to random walk ouput with cell by ppi gene")
    parser.add_argument("-l", "--log_transform",
                        type = bool,
                        default = False,
                        help = "log-transforme input by log(data+exp(1))")
    parser.add_argument("-r", "--restart_path",
                        required = True,
                        help = "path to restart file with headers=[cell, gene]")
    parser.add_argument("-ppi", "--ppi_path",
                        required = True,
                        help = "path to protein protein interaction (i.e., edge list), headers=[source, target, weight]")
    parser.add_argument("-p", "--pathway_path",
                        required = True,
                        help = "path to pathway database in gmt format")
    parser.add_argument("-e", "--expression_path", 
                        required = False,
                        default = None,
                        help = "path to expression file, will be aggregated into rwr probability if given")
    parser.add_argument("-s", "--seed_int",
                        type = int,
                        required = False,
                        default = 42,
                        help = "seed integer, default=42")
    parser.add_argument("-n", "--permutation_int",
                        type = int,
                        required = False,
                        default = 1000,
                        help = "permutation integer, default=1000")
    parser.add_argument("-c", "--cpu_int",
                        type = int, 
                        required = False,
                        default = 5,
                        help = "cpu integer, default=5")
    parser.add_argument("-o", "--out_path",
                        required = True,
                        help = "path to output")
    parser.add_argument("-debug", "--DEBUG",
                        default = True,
                        type = bool,
                        help = "display print message if True")
    return parser.parse_args()


if __name__ == "__main__":
    # timer
    datetimeFormat = '%Y-%m-%d %H:%M:%S.%f'
    start_time = datetime.now()
    
    # get args
    args = parse_parameter()

    # perform Random Walk
    print(datetime.now(), 'performing random walk with restart')
    rwr_df = rwr.RWR(args.ppi_path, args.restart_path, restartProbFloat=0.5, convergenceFloat=0.00001, normalize='l1', weighted=True, outPathStr=args.out_path).get_prob()
    print(rwr_df)
    # multiply with gene expression
    if args.expression_path != None:
        print(datetime.now(), 'multiplying gene expression with random walk probability for genes were expressed')
        exp_df = pd.read_csv(args.expression_path, header=0, index_col=0, sep="\t")
        rwr_df = times_expression(rwr_df, exp_df)
    rwr_df.to_csv(args.out_path+'.RWR.txt', header=True, index=True, sep="\t")
    # perform Pathwa Enrichment Analysis
    print(datetime.now(), 'performing network-based pathway enrichment')
    cell_pathway_df = pea.NetPEA(rwr_df, args.pathway_path, log_transform=args.log_transform, permutation=args.permutation_int, seed=args.seed_int, n_cpu=args.cpu_int, out_path=args.out_path)
    spend = datetime.strptime(str(datetime.now()), datetimeFormat) - datetime.strptime(str(start_time),datetimeFormat)
    print( '[Finished in {:}]'.format(spend) )








