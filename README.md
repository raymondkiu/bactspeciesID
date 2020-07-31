# BACTspeciesID
Fast microbial species identification (16S rRNA gene-based approach) using genome assemblies. This software is run using [ABRicate](https://github.com/tseemann/abricate) for gene screening on contigs. BACTspeciesID also checks for potential contaminations on the whole genome assemblies.

## Dependencies - can be installed using Conda
* ABRicate v1.0.1 (https://github.com/tseemann/abricate) with all its dependecies such as Blast+ v2.2.30, any2fasta
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
$ bactspeciesID.sh -h
bactspeciesID identifies bacterial species/potential contaminations using whole genome assemblies

Usage: ./bactspeciesID2.sh [options] FASTA

Options:
 -i BLASTn identity (default:99)
 -d ABRicate database (default:SILVA-16S)
 -c BLASTn coverage (default:50)
 -m contamination check TRUE/FALSE (default:FALSE)
 -r removal of intermediary files TRUE/FALSE (default:TRUE)
 -h print usage and exit
 -a print author and exit
 -v print version and exit

Version 1.2 (2020)
Author: Raymond Kiu Raymond.Kiu@quadram.ac.uk
```
### Input
Multi-fasta genome assemblies, one genome assembly at a time. Can be any bacterial species.

### Run the software
You can specify BLASTn identity and ABRicate database if you like, 16S rRNA species boundary is recommended at 98.6%, so 99% is to play safe (default parameter anyway). SILVA database has more than 100K sequences and manually curated so it is the recommended database to use. I have tested on >70 samples from multiple species e.g. *Bifidobacterium breve, Bifidobacterium longum, Staphylococcus spp, E. coli spp, Citrobacter spp* etc, and achieved 100% accuracy based on ANI (>95%)support (compared with type strains). Can be used as a quick preliminary analysis.
Importantly, bactsepciesID extracts all 16S sequences if you use -m option, it will then tell you whether this genome is contaminated based upon 16S gene sequence comparison. If there are 16S originated from >1 species it is deemed as contaminated genome.
-i and -d are optional, if not specified it will run at default parameters.
```
$ bactspeciesID.sh -i 99 -d SILVA-16S -m TRUE FASTA.fasta

will identify with BLAST identity 99% and coverage 50% with ABRicate database SILVA-16S
identifying 16S rRNA genes from genome assembly PH102.fna ...
16S rRNA genes have now been identified
extracting 16S rRNA gene from PH102.fna ...
[checking for potential contamination]
index file PH102.fna.fai not found, generating...
....still extracting...
16S rRNA gene is now extracted
[identifying species using ABRicate SILVA-16S database at identity >99%...]
Species identified is now stored in PH102.fna.species
[intermediary files have been removed]

-----------------------------------------------------
The species identified for genome assembly PH102.fna is: 
Bifidobacterium bifidum
Clostridium perfringens
------------------------------------------------------------------
[contamination check: this genome is potentially contaminated :( ]
```
### Output
The result will be shown on stdout, also it will be saved into a file automatically called FASTA.fna.species (where FASTA is your genome assembly's name)

## Issues
This script has been tested on Linux OS, it should run smoothly if dependencies are properly installed. Please report any issues to the [issues page](https://github.com/raymondkiu/bactspeciesID/issues).

## Citation
If you use BACTspeciesID for results in your publication, please cite:
* Kiu R, *BACTspeciesID*, **Github** `https://github.com/raymondkiu/bactspeciesID`

## License
[GPLv3](https://github.com/raymondkiu/bactspeciesID/blob/master/LICENSE)

## Author
Raymond Kiu | Raymond.Kiu@quadram.ac.uk | [@raymond_kiu](https://twitter.com/raymond_kiu)
