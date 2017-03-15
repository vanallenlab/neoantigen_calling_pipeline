S
# Claire Margolis
# 3 October 2016
# getMutNeoantigenBindersMHCClassII.sh
# 
# Summary: Shell script that preprocesses patient .maf files and runs NetMHCIIPan on them
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
use MySQL-5.6

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
python /xchip/cga_home/margolis/mutationsToNeoantigen/PipelineUpdated/mafToFasta.py ../../test_mafs/testmaf.maf 0 14,15 $PAT_DIR ./$PAT_DIR
python /xchip/cga_home/margolis/mutationsToNeoantigen/PipelineUpdated/mafToFasta.py ../../test_mafs/indel.maf 1 14,15 $PAT_DIR ./$PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run runNetMHCpan.py for each patient
# ( Runs netMHCpan program to get predicted binding affinities for each peptide based on patient 
# HLA type )
# *NOTE*: 1 for NetMHCPan, 2 for NetMHCIIPan. Must run script twice if you want both. 

echo 'Running runNetMHCpan.py script.'
python /xchip/cga_home/margolis/mutationsToNeoantigen/PipelineUpdated/runNetMHCpan.py ./$PAT_DIR/len14pep_FASTA_snv.txt,./$PAT_DIR/len14pep_FASTA_indel.txt,./$PAT_DIR/len15pep_FASTA_snv.txt,./$PAT_DIR/len15pep_FASTA_indel.txt ./$PAT_DIR/hla_alleles_mhcII.txt 14,14,15,15 2 ./$PAT_DIR

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Run mutationPostProcess.py for each patient
# ( Processes netMHCpan output to make a more user-friendly file incorporating both mutant and 
# wild-type data for each peptide )
# *NOTE*: 1 for NetMHCPan, 2 for NetMHCIIPan

echo 'Running mutationPostProcess.py script.'
python /xchip/cga_home/margolis/mutationsToNeoantigen/PipelineUpdated/mutationPostProcess.py ./$PAT_DIR/NETMHCIIpan_out_14SNV.xls,./$PAT_DIR/NETMHCIIpan_out_14InDel.xls,./$PAT_DIR/NETMHCIIpan_out_15SNV.xls,./$PAT_DIR/NETMHCIIpan_out_15InDel.xls ./$PAT_DIR/len14pep_headermap_snv.txt,./$PAT_DIR/len14pep_headermap_indel.txt,./$PAT_DIR/len15pep_headermap_snv.txt,./$PAT_DIR/len15pep_headermap_indel.txt 14,14,15,15 $PAT_DIR 2 ./$PAT_DIR/

# ----------------------------------------------------------------------------------------------- #

