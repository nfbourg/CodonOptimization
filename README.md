# CodonOptimization
Codon Optimization

This is a codon optimization algorithm to optimize the payload in a gene editing/gene therapy program. The primary metric used is the BAI, and the algorithm works via Dynamic Programming to generate the sequecnce with the Highest BAI. 

The optimizer is the main function to import the optimizer whereas the metric generator is used to compare sequences after optimization. 

Example is with the co_template.ipynb

## Sources
### Source for CAI/BAI calculation:
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC340524/

### The reference tables used for codon optimization were downloaded and split from TissueCoCoPUTS:
- https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=tisspec
- https://www.sciencedirect.com/science/article/abs/pii/S0022283620300413?via%3Dihub

### Source for SpliceAI toolkit, used for Splice Acceptor/Donor prediciton:
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8360004/
- https://github.com/Illumina/SpliceAI

### Source for using TissueCoCoPUTS:
- https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-00968-8#Tab3

### Source for usage of bicodons over single codons: 
- https://microbialcellfactories.biomedcentral.com/articles/10.1186/s12934-021-01696-y#Sec1
