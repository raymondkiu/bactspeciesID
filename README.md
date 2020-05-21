# speciesID
Fast microbial species identification (16S rRNA gene-based approach) using genome assemblies. This software is run via [ABRicate](https://github.com/tseemann/abricate) gene screening on contigs.

## Dependencies - can be installed using Conda
* ABRicate v1.0.1 (https://github.com/tseemann/abricate) with all its dependecies such as Blast+ v2.2.30, any2fasta, Emboss etc, also  need to build a 16S rRNA database (recommend SILVA database)
* Barrnap v0.7 (https://github.com/tseemann/barrnap)
* Bedtools
* Samtools

## 16S rRNA Database
* SILVA 16S database can be [downloaded from here](https://zenodo.org/record/3731176/files/silva_species_assignment_v138.fa.gz?download=1) with format
```
>HG530238.1.1461 Paucibacter toxinivorans
TCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGG
```

## To set up sequence database in ABRicate
Please refer to this section in [ABRicate repository](https://github.com/tseemann/abricate#making-your-own-database) and rename the database as SILVA-16S if you use the SILVA database (alternatively, any name you like).

## Usage
### Options
```
$ speciesID.sh -h
This script identifies bacterial species using genome assemblies

Usage: ./speciesID.sh [options] FASTA
Option:
 -i BLASTn identity (default:99)
 -c BLASTn coverage (default:50)
 -d ABRicate database (default:SILVA-16S)
 -h print usage and exit
 -a print author and exit
 -v print version and exit

Version 1.1 (2020)
Author: Raymond Kiu Raymond.Kiu@quadram.ac.uk
```
### Input
Multi-fasta genome assemblies, one file at a time. Can be any bacterial species.

### Run the software
You can specify BLASTn identity and ABRicate database if you like, 16S rRNA species boundary is recommended at 98.6%, so 99% is to play safe (default parameter anyway). SILVA database has more than 100K sequences and manually curated so it is the recommended database to use. I have tested on >70 samples from multiple species e.g. *Bifidobacterium breve, Bifidobacterium longum, Staphylococcus spp, E. coli spp, Citrobacter spp* etc, and achieved 100% accuracy based on ANI (>95%)support (compared with type strains). Can be used as a quick preliminary analysis.
This script will only extract one 16S rRNA gene in the genomes even if there is multiple 16S genes in one genome, so could not detect contamination and is not built for that purpose.
-i and -d are optional, if not specified it will run at default parameters.
```
$ speciesID.sh -i 99 -d SILVA-16S FASTA

will identify with BLASTn identity 99% and BLASTn coverage 50% with ABRicate database SILVA-16S
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
### Output
The result will be shown on stdout, also it will be saved into a file automatically called FASTA.fna.species (where FASTA is your genome assembly's name)

## Issues
This script has been tested on Linux OS, it should run smoothly if dependencies are properly installed. Please report any issues to the [issues page](https://github.com/raymondkiu/speciesID/issues).

## License
[GPLv3](https://github.com/raymondkiu/speciesID/blob/master/LICENSE)

## Author
Raymond Kiu Raymond.Kiu@quadram.ac.uk
