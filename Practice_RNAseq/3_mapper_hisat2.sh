#!/bin/sh
 
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference
##              First need to use gffread to convert annotation file from .gff3 to .gft formate
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2  Indexing  InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome 
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort         Input: Alignment files  .sam
##                                                  Output: Sorted .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
## For running the script on the Alabama Super Computer.
##  For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
##    queue: class
##    core: 6
##    time limit (HH:MM:SS): 04:00:00 
##    Memory: 12gb
##    run on dmc
###############################################

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load hisat2
module load stringtie/2.2.1
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread/
module load gffcompare/


#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## Replace the [#] with paths to define these variable
MyID=[1]          ## Example: MyID=aubtss

WD=[2]                                                ## Example:/scratch/$MyID/PracticeRNAseq  
CLEAND=[3]                                            ## Example:/scratch/$MyID/PracticeRNAseq/CleanData20   #   *** This is where the cleaned paired files are located
REFD=[4]                                              ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome_6                      # this directory contains the indexed reference genome for the garter snake
MAPD=[5]                                              ## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2_6
COUNTSD=[6]                                           ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie_6

RESULTSD=[7]                                          ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S_6

REF=DaphniaPulex_RefGenome_PA42_v3.0                  ## This is what the "easy name" will be for the genome reference

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
cp /home/shared/schwartz_class/$REF.fasta .
cp /home/shared/schwartz_class/$REF.gff3 .

###  Identify exons and splice sites
gffread $REF.gff3 -T -o $REF.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
extract_splice_sites.py $REF.gtf > $REF.ss
extract_exons.py $REF.gtf > $REF.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss $REF.ss --exon $REF.exon $REF.fasta DpulPA42_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd $CLEAND  #### This is where our clean paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
cd $MAPD

## move the list of unique ids from the original files to map
mv $CLEAND/list  . 

while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "$REFD"/DpulPA42_index       \
    -1 "$CLEAND"/"$i"_1_paired.fastq  -2 "$CLEAND"/"$i"_2_paired.fastq      \
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  ### This works on ASC

    ###  This is sorting the bam
  samtools sort -@ 6  "$i".bam    "$i"_sorted

    ### Index the BAM and get mapping statistics
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  ###
mkdir "$COUNTSD"/"$i"
stringtie -p 6 -e -B -G  "$REFD"/"$REF".gtf -o "$COUNTSD"/"$i"/"$i".gtf -l "$i"   "$MAPD"/"$i"_sorted.bam

done<list

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt $RESULTSD

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix
cd $COUNTSD
python /home/$MyID/class_shared/prepDE.py -i $COUNTSD

### copy the final results files (the count matricies that are .cvs to your home directory)
cp *.csv $RESULTSD
