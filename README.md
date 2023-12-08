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

Copyright (c) 2023 Nathanael Bourgeois-Tchir nfbourg@me.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
