# pangenome

Under active development. This software is a sequence and gene aware pangenome construction toolkit. The approach used here is iterative.

0. The assembly with the most sequence is considered the "current working panreference" at the start.
1. Take the longest, unprocessed assembly and GFF3 file, when available.
2. Compare using nucmer to the current working panreference. Add in all non-redundant sequences as separate contigs to the current working panreference. 
3. If a GFF3 file is provided, extended any non-redundant sequences to existing feature (gene) boundaries so that we do not carry partial features to the next iteration.
4. Add any features, and rename when appropriate, to the current work panreference GFF3 file. 
5. If there are more unprocessed assemblies, GOTO 1.

Until recently, most pan-genomes have consisted solely of the genes, while ignoring the rest of the genomic context. There are advantages to identifying all sequence within a population, such as:

* Expanded SNP Calling
* Identify large gene families
* Identify sequence of structural variants
* Expanded sequence 
* Population modalities of sequences

Work based off of a single reference may exclude many genes and sequences from downstream population studies. The construction of a panreference reduces this problem. 

TODO:
* Web interface
* Better statistic tracking
* Load into Neo4j database
* Identify sequences that are singletons, and which sequence chunks are found in multiple, but not all, population members
* Identify conserved sequences
* Identify structural variations such as translocations, inversion, etc...

## Usage

Not yet

## License

Copyright © 2017 Joseph Guhlin

GPLv3
