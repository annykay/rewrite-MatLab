# rewrite_MatLab
Here is a program in python, that calculates p- and q- values for given expression values with respect to the treshold dstributions. 

Treshold distribution array derived from DetermineTresholdDestribution.m (in this repository - Treshold_CoTFcRF.m), gene expression values matrix is derived from HSC_bulk_exp.tsv.

To calculate these values use BooleanizeFromFile.py function. It takes two arguments: "Inputdata" and "Outputfile". "Inputdata" variable schould contain the way to the .tsv file with expression data. In our example file with bulk expression data contained of three columnes: ensemble id, hgnc symbol and expression level. Also, file with treshhold distribution is needed. File background_genes schould contain names of genes, expressed in our sample. The results will be written to "Outputfile". It would contain two columnes: the name of the gene and value 1 or 0. 

Before importing this package, you schould add the way to this files to PATH. 

```python
python -c 'from rewrite_MatLab.BooleanizeFromFile import BooleanizeFromFile as BooleanizeFromFile; BooleanizeFromFile('HSC_bulk_exp.tsv', 'test.txt')'
```
