# NetPEA
Re-implementation of Network-Based Pathway Enrichment Analysis by Liu et al., 2013

## Reference
Liu, L., & Ruan, J. (2013). Network-based pathway enrichment analysis. In 2013 IEEE International Conference on Bioinformatics and Biomedicine (pp. 218–221).



## USAGE
$ python NetPEA/run_netpea.py -h
usage: run_netpea.py [-h] [-l LOG_TRANSFORM] -r RESTART_PATH -ppi PPI_PATH -p
                     PATHWAY_PATH [-e EXPRESSION_PATH] [-s SEED_INT]
                     [-n PERMUTATION_INT] [-c CPU_INT] -o OUT_PATH
                     [-debug DEBUG]

Return NetPEA zscore

optional arguments:
  -h, --help            show this help message and exit
  -l LOG_TRANSFORM, --log_transform LOG_TRANSFORM
                        log-transforme input by log(data+exp(1))
  -r RESTART_PATH, --restart_path RESTART_PATH
                        path to restart file with headers=[cell, gene]
  -ppi PPI_PATH, --ppi_path PPI_PATH
                        path to protein protein interaction (i.e., edge list),
                        headers=[source, target, weight]
  -p PATHWAY_PATH, --pathway_path PATHWAY_PATH
                        path to pathway database in gmt format
  -e EXPRESSION_PATH, --expression_path EXPRESSION_PATH
                        path to expression file, will be aggregated into rwr
                        probability if given
  -s SEED_INT, --seed_int SEED_INT
                        seed integer, default=42
  -n PERMUTATION_INT, --permutation_int PERMUTATION_INT
                        permutation integer, default=1000
  -c CPU_INT, --cpu_int CPU_INT
                        cpu integer, default=5
  -o OUT_PATH, --out_path OUT_PATH
                        path to output
  -debug DEBUG, --DEBUG DEBUG
                        display print message if True
