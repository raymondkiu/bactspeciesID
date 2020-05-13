# speciesID
Fast microbial species identification (16S rRNA gene-based approach) using genome assemblies

## Dependencies - can be installed using Conda
* ABRicate v1.0.1 (https://github.com/tseemann/abricate) with all its dependecies such as Blast+ v2.2.30, any2fasta, Emboss etc, also  need to build a 16S rRNA database (recommend SILVA database)
* Barrnap v0.7 (https://github.com/tseemann/barrnap)
* Bedtools
* Samtools
* SILVA 16S database can be downloaded from https://zenodo.org/record/3731176/files/silva_species_assignment_v138.fa.gz?download=1 with format
```
>HG530238.1.1461 Paucibacter toxinivorans
TCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGG
```

## Usage
```
$ speciesID.sh -h
This script identifies bacterial species using genome assemblies

Usage: ./speciesID.sh [options] FASTA
Option:
 -i BLAST identity (default:99)
 -d ABRicate database (default:SILVA-16S)
 -h print usage and exit
 -a print author and exit
 -v print version and exit

Version 1.0 (2020)
Author: Raymond Kiu Raymond.Kiu@quadram.ac.uk
```
For example
```
$ speciesID.sh -i 99 -d SILVA-16S FASTA

will identify with BLAST identity 99% with ABRicate database SILVA-16S
identifying 16S rRNA genes from genome assembly CA-20.fna ...
16S rRNA genes have now been identified
extracting 16S rRNA gene from CA-20.fna ...
index file CA-20.fna.fai not found, generating...
....still extracting...
16S rRNA gene is now extracted
removing intermediary files...
16S rRNA gene sequence is now stored in CA-20.fna-16S.fna
identifying species using ABRicate SILVA-16S database at identity >99%...
Species identified is now stored in CA-20.fna.species
Removing intermediary file...
Programme will now exit
The species identified for CA-20.fna is 
Bifidobacterium bifidum
```
