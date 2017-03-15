# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# 2 February 2017
# getMutNeoantigenBinders.sh
# 
# Summary: Shell script that preprocesses patient .maf files and runs netMHCpan on them. 
# 	   ** Updated to fix errors involving translating into/out of intronic regions as opposed to
#	   using coding sequences only. ** 
#
# *NOTE*: If you want to run this script, go through and verify that the paths to relevant files
# are in the correct format for your cohort. You will need to change out_dir.txt, among other 
# things, to make the script specific to your cohort. You can also change preferences for running
# netMHCIpan vs. netMHCIIpan. 

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

patient_dir=test_dirs.txt
PAT_DIR=$(cat $patient_dir | head -n $SGE_TASK_ID | tail -n 1)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run mafToFasta.py for each patient for both SNVs and indels
# ( Converts mutations in maf file to mutant peptides, generates wild-type peptides as well, and 
# writes both to outfile )

echo 'Running mafToFasta.py script for both SNVs and indels.'
python mafToFastaV2.py ../../test_mafs/testmaf.maf 0 9,10 $PAT_DIR ./$PAT_DIR
python mafToFastaV2.py ../../test_mafs/indel.maf 1 9,10 $PAT_DIR ./$PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run runNetMHCpan.py for each patient
# ( Runs netMHCpan program to get predicted binding affinities for each peptide based on patient 
# HLA type )
# *NOTE*: 1 for NetMHCPan, 2 for NetMHCIIPan. Must run script twice if you want both. 

echo 'Running runNetMHCpan.py script.'
python runNetMHCpan.py ./$PAT_DIR/len9pep_FASTA_snv.txt,./$PAT_DIR/len9pep_FASTA_indel.txt,./$PAT_DIR/len10pep_FASTA_snv.txt,./$PAT_DIR/len10pep_FASTA_indel.txt ./$PAT_DIR/hla_alleles.txt 9,9,10,10 1 ./$PAT_DIR

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Run mutationPostProcess.py for each patient
# ( Processes netMHCpan output to make a more user-friendly file incorporating both mutant and 
# wild-type data for each peptide )
# *NOTE*: 1 for NetMHCPan, 2 for NetMHCIIPan

echo 'Running mutationPostProcess.py script.'
python mutationPostProcess.py ./$PAT_DIR/NETMHCpan_out_9SNV.xls,./$PAT_DIR/NETMHCpan_out_9InDel.xls,./$PAT_DIR/NETMHCpan_out_10SNV.xls,./$PAT_DIR/NETMHCpan_out_10InDel.xls ./$PAT_DIR/len9pep_headermap_snv.txt,./$PAT_DIR/len9pep_headermap_indel.txt,./$PAT_DIR/len10pep_headermap_snv.txt,./$PAT_DIR/len10pep_headermap_indel.txt 9,9,10,10 $PAT_DIR 1 ./$PAT_DIR/

# ----------------------------------------------------------------------------------------------- #

