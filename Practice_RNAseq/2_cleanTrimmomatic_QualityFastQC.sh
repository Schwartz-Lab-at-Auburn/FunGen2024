#! /bin/bash

######## FunGen Course Instructions ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##       and then use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)       
## FASTQC output is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##		Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics.
##			The last line of this script will make a tarball of the output directory to bring back to your computer
## For running the script on the Alabama Super Computer.
                ##For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
        ## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
        ## then run the script by using "run_script [script name]"
        ## suggested paramenters are below to submit this script.
                ## queue: class  
                ## core: 6
                ## time limit (HH:MM:SS): 02:00:00  (may need to increase, if so run on medium queue)
                ## Memory: 12gb
                ## run on dmc
###############################################

## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the read data.
## Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
## Output Data: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## More Information: http://www.usadellab.org/cms/?page=trimmomatic

# Modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load trimmomatic/0.39
module load fastqc/0.10.1

## STOP. You need to replace the [number] with YOUR paths to 
##       make variables for your ASC ID so the directories are automatically made in YOUR directory
MyID=[1]                        ## Example: MyID=aubtss

# Variables: raw data directory (DD), working directory(WD), cleaned status (CS), name of file containing the adpaters.
WD=[2]                          ## Example: WD=/scratch/$MyID/PracticeRNAseq
DD=[3]                          ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
CD=[4]                          ## Example: CD=/scratch/$MyID/PracticeRNAseq/CleanData
CS=PostCleanQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. 
					## You will likely need to edit this for other projects based on how your libraries 
					## were made to search for the correct adapters for your project
				
## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
mkdir $CD
mkdir $WD/$CS

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd $DD

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
cp /home/$MyID/class_shared/AdaptersToTrim_All.fa . 

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77

        java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar  PE -threads 6 -phred33 \
        "$i"_1.fastq "$i"_2.fastq  \
        $CD/"$i"_1_paired.fastq $CD/"$i"_1_unpaired.fastq  $CD/"$i"_2_paired.fastq $CD/"$i"_2_unpaired.fastq \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
                ## requiredQuality: specifies the average quality required.

	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

fastqc $CD/"$i"_1_paired.fastq --outdir=$WD/$CS
fastqc $CD/"$i"_2_paired.fastq --outdir=$WD/$CS

done<list			# This is the end of the loop

#########################  Now compress your results files from the Quality Assessment by FastQC 
## move to the directory with the cleaned data
cd $WD/$CS

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate
tar cvzf $CS.tar.gz $WD/$CS/*
