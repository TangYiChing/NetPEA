# NetPEA
Re-implementation of Network-Based Pathway Enrichment Analysis by Liu et al., 2013

## Reference
Liu, L., & Ruan, J. (2013). Network-based pathway enrichment analysis. In 2013 IEEE International Conference on Bioinformatics and Biomedicine (pp. 218–221).



## USAGE
1. perform NetPEA with gene expression (use #cpu=1, permutation test=10 times, output folder&prefix=./example_run)

```python
$ python /NetPEA/NetPEA/run_netpea.py -r /NetPEA/NetPEA/example_data/GeneMutation.Mat.txt -ppi /NetPEA/NetPEA/example_data/9606.protein_name.links.v11.0.pkl -p /NetPEA/NetPEA/example_data/c2.cp.pid.v7.1.symbols.gmt -e /NetPEA/NetPEA/example_data/GeneExpression.Mat.txt -c 1 -n 10 -o ./example_run
```

2. perform NetPEA without gene expression (use #cpu=1, permutation test=10 times, output folder&prefix=./example_run)

```python
$ python /NetPEA/NetPEA/run_netpea.py -r /NetPEA/NetPEA/example_data/GeneMutation.Mat.txt -ppi /NetPEA/NetPEA/example_data/9606.protein_name.links.v11.0.pkl -p /NetPEA/NetPEA/example_data/c2.cp.pid.v7.1.symbols.gmt -c 1 -n 10 -o ./example_run
```
