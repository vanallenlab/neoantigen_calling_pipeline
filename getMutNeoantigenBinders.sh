# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# 21 March 2017
# getMutNeoantigenBinders.sh
# 
# Summary: Shell script that calls all python scripts to run neoantigen pipeline on batch of 
#	   samples. Intended to be submitted to UGER from the command line with task array 
#	   of samples. 
#	   Sample usage: qsub getMutNeoantigenBinders.sh
#
# *NOTE*: If you want to run this script, go through and verify that the paths to relevant files
# are in the correct format for your cohort. You will need to change out_dir.txt, among other 
# things, to make the script specific to your cohort.

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Specify shell / UGER preferences

#!/bin/bash

#$ -cwd
#$ -q long
#$ -m e
#$ -l h_vmem=10g 
#$ -t 1-2
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Use statements

source /broad/software/scripts/useuse
reuse Python-2.7
reuse MySQL-5.6
reuse R-3.3

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Set directory paths

patient_dirs=./patient_dirs.txt # File should contain sample names (which double as directory names), one per line
hla_paths=./hla_paths.txt # File should contain paths to HLA allele files for each sample, one per line, in same order as patient_dirs.txt file
snv_maf_paths=./snv_maf_paths.txt # File should contain paths to MuTect files for each sample, one per line, in same order as patient_dirs.txt file
indel_maf_paths=./indel_maf_paths.txt # File should contain paths to Strelka files for each sample, one per line, in same order as patient_dirs.txt file
PAT_DIR=$(cat $patient_dirs | head -n $SGE_TASK_ID | tail -n 1)
HLA_PATH=$(cat $hla_paths | head -n $SGE_TASK_ID | tail -n 1)
SNV_MAF_PATH=$(cat $snv_maf_paths | head -n $SGE_TASK_ID | tail -n 1)
INDEL_MAF_PATH=$(cat $indel_maf_paths | head -n $SGE_TASK_ID | tail -n 1)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run mafToFasta.py for each patient for both SNVs and indels
# ( Converts mutations in maf file to mutant peptides, generates wild-type peptides as well, and 
# writes both to outfile )

echo 'Running mafToFasta.py script for both SNVs and indels.'
python mafToFastaV2.py $SNV_MAF_PATH 0 9,10 $PAT_DIR ./$PAT_DIR  # Only change last argument (./$PAT_DIR) to contain whatever output path desired
python mafToFastaV2.py $INDEL_MAF_PATH 1 9,10 $PAT_DIR ./$PAT_DIR  # Only change last argument to contain whatever output path you want

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run runNetMHCpan.py for each patient
# ( Runs netMHCpan program to get predicted binding affinities for each peptide based on patient 
# HLA type )

echo 'Running runNetMHCpan.py script.'
python runNetMHCpan.py ./$PAT_DIR/len9pep_FASTA_snv.txt,./$PAT_DIR/len9pep_FASTA_indel.txt,./$PAT_DIR/len10pep_FASTA_snv.txt,./$PAT_DIR/len10pep_FASTA_indel.txt $HLA_PATH 9,9,10,10 1 ./$PAT_DIR  # Change "./$PAT_DIR" parts if you wrote to a different output path above

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Run mutationPostProcess.py for each patient
# ( Processes netMHCpan output to make a more user-friendly file incorporating both mutant and 
# wild-type data for each peptide )

echo 'Running mutationPostProcess.py script.'
python mutationPostProcess.py ./$PAT_DIR/NETMHCpan_out_9SNV.xls,./$PAT_DIR/NETMHCpan_out_9InDel.xls,./$PAT_DIR/NETMHCpan_out_10SNV.xls,./$PAT_DIR/NETMHCpan_out_10InDel.xls ./$PAT_DIR/len9pep_headermap_snv.txt,./$PAT_DIR/len9pep_headermap_indel.txt,./$PAT_DIR/len10pep_headermap_snv.txt,./$PAT_DIR/len10pep_headermap_indel.txt 9,9,10,10 $PAT_DIR 1 ./$PAT_DIR/  # Change "./$PAT_DIR" parts if you wrote to a different output path above

# ----------------------------------------------------------------------------------------------- #

