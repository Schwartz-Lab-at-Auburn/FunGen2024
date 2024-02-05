#!/bin/sh

############################################################################################
###  	Title: Practice Script 1						
### 	Author: Tonia Schwartz
###		Date: January, 2024
###		BIOL6850: Functional Genomics, Auburn University
###		Purpose: Learn to make a directory (folder) in scratch, move files, assign variables,
###			and check error files. There are 4 errors in this file. Can you find them and make it run correctly?
##############################################################################################

#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load <tool name >

######### Your first goal is to make a directory in scratch for you to conduct your work #######
##### Assign a variable to the directory name that you plan to make. 
	## This will make it easier in following steps so you don't have to write out the whole directory everytime
	##  IMPORTANT! 'aubtss' is MY (Tonia Schwartz) identifier. You need to replace 'aubtss' with YOUR ID
DATADIR=/scratch/aubtss/test2
SHAREDIR=/home/aubtss/class_shared/

######  Now use that variable to make the directory in SCRATCH for holding your data
###  Example: mkdir /scratch/YOUR_ID/fastqc
###  -p means to make all directories above if needed
mkdir -p ${DATADIR}

### Change (move) to the scratch directory you just made
### Example: cd /scratch/YOUR_ID/test
cd ${DATADIR}
echo "Yea, directory made! Line 30 completed" > Notes.txt

####### Move the practice datafiles (all files with .fastq) from our shared directory to where you are, here (.)
cp ${SHAREDIR}/*.fastq  .


###########  Check for errors in transfer. Calculate the md5sum for the file.
### Calculate the md5sum values of the file in the original folder and read (put) into a text file.
md5sum {SHAREDIR}*.fastq > md5sum_Original.txt

### Calculate the md5sum values of the files you just moved and read (put) into a text file in our a text file
md5sum ./*.fastq >> md5sum_New.txt

######  Make a directory for this project and results in your home folder
mkdir home/aubtss/Practice_Code_2024/md5sum_files

##### Move the md5sums text files back to your home folder
mv md5sum.txt /home/aubtss/Practice_Code_2024/md5sum_files

mv md5sum_Original.txt /home/aubtss/Practice_Code_2024/md5sum_files

