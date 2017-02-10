# mgblast

Modification of the megablast code from a 2006 version of NCBI C Toolkit
which is specifically designed to efficiently perform an all-vs-all search of DNA sequences
in a multi-FASTA file.

The modified mgblast program also provides customized tab-delimited output which 
can be used for clustering and even assembly of similar DNA sequences (TGICL 
is such a pipeline).

# Installation

Use git to clone this repository - e.g. in a directory called mgblast.
Then:
```
 cd mgblast
 ./build_mgblast.sh
```

