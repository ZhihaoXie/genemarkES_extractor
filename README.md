# genemarkES\_extractor

Tools for taking GeneMarkES output(gff format) and pulling out the relevant sequences from the original fasta file. 

Will generate 3 files about gene:

+ out_prefix.fnn (nucl)
+ out_prefix.faa (prot)
+ out_prefix.gff (only gene info)


# Installation

`git clone git@github.com:ZhihaoXie/genemarkES_extractor.git`


# Usage

```
python3 genemarkES_extractor.py <genemarkES_out(gff)> <genome_fasta> <out_prefix>
```

# Requests

Python3 and Biopython


# Copyright

Zhihao Xie, xiezhihao1122@outlook.com
