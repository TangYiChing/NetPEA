# NetPEA
Re-implementation of Network-Based Pathway Enrichment Analysis by Liu et al., 2013

## Reference
Liu, L., & Ruan, J. (2013). Network-based pathway enrichment analysis. In 2013 IEEE International Conference on Bioinformatics and Biomedicine (pp. 218–221).



## USAGE
1. perform NetPEA with gene expression
$ python /NetPEA/run_netpea.py -r /example_data/GeneMutation.Mat.txt -ppi ./example_data/9606.protein_name.links.v11.0.pkl -p ./example_data/c2.cp.pid.v7.1.symbols.gmt -c 1 -n 10 -o ./example_run


