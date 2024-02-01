#!/bin/sh

############################################################################################
###  	Title: PracticeScript.sh 						
### 	Author: Tonia Schwartz
###		Date: January, 2021
###		BIOL6950: Functional Genomics, Auburn University
###		Purpose: Learn and Practice to make scratch directory, move files, assign variables,
###			and check error files.
##############################################################################################

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load <tool name >


######  make a directory in SCRATCH for holding YOUR raw data
##  Example: mkdir /scratch/YOUR_ID/fastqc
mkdir /scratch/aubtss
mkdir /scratch/aubtss/test2024

##### Assign that directory a variable
DATADIR=/scratch/aubtss/test

### Change directory to the scratch directory for the data
### Example: cd /scratch/YOUR_ID/test
cd ${DATADIR}

### move a practice datafiles here (.)
cp /home/aubtss/class_shared/Mini_HeatStress_Data/HS06_GATCAG_*.fastq.gz  .

###########  Check for errors in transfer
### Calculate the md5sum values of the file in the original folder and read into a text file.
md5sum /home/aubtss/class_shared/Mini_HeatStress_Data/HS06_*.fastq.gz > md5sum.txt

### Calculate the md5sum values of all the files and read into a text file
md5sum ./* >> md5sum.txt

### Make a directory for this project and results in my home folder
mkdir -p /home/aubtss/GarterSnake/Results
### move the text file back to your home folder
mv md5sum.txt /home/aubtss/GarterSnake/Results
