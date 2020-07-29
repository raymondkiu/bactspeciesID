#!/bin/bash

#print the options
usage () {
  echo ""
  echo "BactspeciesID identifies bacterial species/potential contaminations using whole genome assemblies"
  echo ""
  echo "Usage: $0 [options] FASTA"
  echo ""
  echo "Options:"
  echo " -i BLASTn identity (default:99)"
  echo " -d ABRicate database (default:SILVA-16S)"
  echo " -c BLASTn coverage (default:50)"
  echo " -m contamination check TRUE/FALSE (default:FALSE)"
  echo " -r removal of intermediary files TRUE/FALSE (default:TRUE)"
  echo " -h print usage and exit"
  echo " -a print author and exit"
  echo " -v print version and exit"
  echo ""
  echo "Version 1.2 (2020)"
  echo "Author: Raymond Kiu Raymond.Kiu@quadram.ac.uk"
  echo "";
}
version () { echo "version 1.2 (2020)";}
author () { echo "Author: Raymond Kiu Raymond.Kiu@quadram.ac.uk";}

#dafault value for BLAST identity and database dor ABRicate

IDENTITY=99
DATABASE=SILVA-16S
COVERAGE=50
CONTAMINATION=FALSE
REMOVAL=TRUE
b=1
a=0

# parse the options
while getopts 'i:d:c:m:r:hav' opt;do
  case $opt in
    i) IDENTITY=$OPTARG ;;
    d) DATABASE=$OPTARG ;;
    c) COVERAGE=$OPTARG ;;
    m) CONTAMINATION=$OPTARG ;;
    r) REMOVAL=$OPTARG ;;
    h) usage; exit;;
    a) author; exit;;
    v) version; exit;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    *) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
   esac
done

# skip over the processed options
shift $((OPTIND-1))
# check for mandatory positional parameters, only 1 positional argument will be checked
if [ $# -lt 1 ]; then
   echo "Usage: $0 [options] FASTA"
   echo "e.g. $0 -m TRUE -d SILVA-16S fasta.fna"
   echo ""
   echo "Options:"
   echo " -i BLASTn identity (default:$IDENTITY)"
   echo " -d ABRicate database (default:$DATABASE)"
   echo " -c BLASTn coverage (default:$COVERAGE)"
   echo " -m contamination check TRUE/FALSE (default:$CONTAMINATION)"
   echo " -r removal of intermediary files TRUE/FALSE (default:$REMOVAL)" 
   echo ""
   exit 1
fi

fasta=$1

echo "------------------------------------------------------------------------------------"
echo "BactspeciesID will identify with the following parameters: "
echo "BLAST identity $IDENTITY% and coverage $COVERAGE% with ABRicate database $DATABASE "
echo "Contamination option is set to $CONTAMINATION "
echo "Intermediary file removal option is set to $REMOVAL "
echo "BactspeciesID will start identifying 16S rRNA genes from genome assembly $fasta "
echo "------------------------------------------------------------------------------------"

barrnap --quiet $1 > $fasta-gff

if [[ "$CONTAMINATION" == "TRUE" ]]
then
        grep '16S' $fasta-gff > $fasta-16S.gff;
else    
        grep -m 1 '16S' $fasta-gff > $fasta-16S.gff;
fi
        # check if there is 16S sequences
        NUM1=$(cat $fasta-16S.gff|wc -l)
        if [ "$NUM1" -eq "$a" ]
        then
                echo "No 16S sequences found, exiting..."
                echo "Thank you for using bactspeciesID!"
                exit 0; 
        else
                echo "$NUM1 16S sequence(s) found, continue..."
        fi

        bedtools getfasta -fi $fasta -bed $fasta-16S.gff -fo $fasta-16S-fasta.fna;
        
        grep ">" $fasta-16S-fasta.fna|sed 's/>//g' > $fasta-16S-id.txt;
        xargs samtools faidx $fasta-16S-fasta.fna < $fasta-16S-id.txt > $fasta-16S.fna

        echo "$NUM1 sequence(s) extracted..."
        
        echo "[Identifying species using ABRicate $DATABASE database at identity >$IDENTITY%...]"
        abricate --quiet --db $DATABASE --mincov=$COVERAGE --minid=$IDENTITY $fasta-16S.fna|grep -v "#FILE" | awk '{print $13" "$14}' > $fasta.species

        echo "Species identity is now stored in $fasta.species"
        if [[ "$REMOVAL" == "TRUE" ]]
        then
                rm $fasta-gff
                rm $fasta-16S.gff
                rm $fasta-16S-fasta.fna
                rm $fasta.fai
                rm $fasta-16S-fasta.fna.fai
                rm $fasta-16S-id.txt
                rm $fasta-16S.fna
                echo "[Intermediary files have been removed]"
        else
                echo "[Intermediary files are NOT removed]"
        fi
        echo ""
        echo "-----------------------------------------------------"
        echo "The species identified for genome assembly $fasta : "
        cat $fasta.species
        
        # checking for potential contamination -gt greater than 
        NUM=$(cat $fasta.species|wc -l)
        if [ "$NUM" -eq "$a" ]
        then 
                echo "--------------------------------------------------------------------------------------]"
                echo "[Unfortunately bactspeciesID did not identify any 16S gene sequence in this genome :( ]"
                echo "Thank you for using bactspeciesID!"
                exit 0; 
        else
                # Checking contamination outcome if > 1 species is identified
                if [[ "$CONTAMINATION" == "TRUE" ]]
                then
                        if [ "$NUM" -gt "$b" ]
                        then
                                echo "------------------------------------------------------------------"
                                echo "[contamination check: this genome is potentially contaminated :( ]"
                                echo "Thank you for using bactspeciesID!"
                        else
                                echo "----------------------------------------------------------------------"
                                echo "[contamination check: this genome is NOT known to be contaminated :) ]"
                                echo "Thank you for using bactspeciesID!"
                        fi
                else
                        echo "Thank you for using bactspeciesID!"
                        exit 0;
                fi
        fi
exit 0;
