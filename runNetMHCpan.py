# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# runNetMHCpan.py
#
# Summary: Takes in one or more FASTA files containing all of the peptides upon which netMHCpan
# or netMHCIIpan is to be run. Runs whichever version of netMHCpan is requested and returns
# the results in an appropriately-named output file.
#
# Input format: python runNetMHCpan.py len9peptides.txt,len10peptides.txt HLAalleles.txt 1 outpath
# 	Options for specifying which netMHCpan version: 
#	1 = netMHCIpan
#	2 = netMHCIIpan
# *RELEVANT*: HLA allele input file can be in one of two formats: 
#	1. Polysolver winners_hla.txt output file
# 		example line from file: HLA-A   hla_a_02_01_01_01       hla_a_32_01_01
# 	2. Already processed, one allele per line in netMHC compatible format
#		example line from file: HLA-A02:01
#
# Output: netMHCpan output .xls file(s)
#
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: runNetMHCIpan
# Inputs: FASTA file of peptide sequences, patient HLA alleles (these are automatically given 
# by Polysolver and come in a .txt file that needs to be pre-processed into the correct format for 
# netMHCpan), peptide length outpath 
# Returns: None (netMHCpan will automatically write output to a .xls file)
# Summary: Pre-processes patient HLA alleles, runs netMHCIpan. 
def runNetMHCIpan(pepfile, hlafile, length, outpath):
	# Determine whether we're dealing with a snv or indel file (for naming the outfile)
	varianttype = ''
	if pepfile.split('_FASTA_')[1].split('.')[0] == 'snv':
		varianttype = 'SNV'
	if pepfile.split('_FASTA_')[1].split('.')[0] == 'indel':
		varianttype = 'InDel'
	# Read in HLA alleles file and process
	with open(hlafile) as f:
		hlalines  = f.read().splitlines()
	hlaalleles = []
	# Determine which input format the hla allele file is in
	if len(hlalines[0].split('\t')) <= 1:  # In already pre-processed format
		hlaalleles = hlalines
	else:  # Polysolver output file
		for line in hlalines:
			split = line.split('\t')
			# Reformat each allele (2 for each type of HLA A, B, and C)
			for i in range(1, 3):
				currallele = 'HLA-'
				allele = split[i]
				components = allele.split('_')
				currallele += components[1].upper() + components[2] + ':' + components[3]
				hlaalleles.append(currallele)
	hlaalleles = list(set(hlaalleles))  # Remove duplicate alleles if there are any
	hlastring = ','.join(hlaalleles)
	# Run netMHCI pan
	command =  'export NHOME=/xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0; export NETMHCpan=/xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64; /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/bin/netMHCpan -a '+hlastring+' -f '+pepfile+' -inptype 0 -l '+str(length)+' -s -xls -xlsfile '+outpath+'/NETMHCpan_out_'+str(length)+varianttype+'.xls -allname /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/data/allelenames -hlapseudo /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/data/MHC_pseudo.dat -t 500 -version /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/data/version -tdir /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/scratch/XXXXXX -rdir /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/ > '+outpath+'/netMHCpanoutlen_'+str(length)+varianttype+'.txt'
	subprocess.call(command, shell=True)
	
	# Catch case where peptide file was empty (create dummy file) 
	dummyfile = outpath+'/NETMHCpan_out_'+str(length)+varianttype+'.xls'
	open(dummyfile, 'a').close()

	return

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: runNetMHCIIpan
# Inputs: FASTA file of peptide sequences, patient HLA alleles (these are automatically given 
# by Polysolver and come in a .txt file that needs to be pre-processed into the correct format for 
# netMHCIIpan), peptide length outpath 
# Returns: None (netMHCIIpan will automatically write output to a .xls file)
# Summary: Pre-processes patient HLA alleles, runs netMHCIIpan
def runNetMHCIIpan(pepfile, hlafile, length, outpath):
	# Determine whether we're dealing with a snv or indel file (for naming the outfile)
        varianttype = ''
        if pepfile.split('_FASTA_')[1].split('.')[0] == 'snv':
                varianttype = 'SNV'
        if pepfile.split('_FASTA_')[1].split('.')[0] == 'indel':
                varianttype = 'InDel'
        # Read in HLA alleles file and process
        with open(hlafile) as f:
                hlalines = f.read().splitlines()
        hlaalleles = []
        # Determine which input format the hla allele file is in
        if len(hlalines[0].split('\t')) <= 1:  # In already pre-processed format
                hlaalleles = hlalines
        else:  # PHLAT output file
		# DQA1
		DQA1a = hlalines[4].split('\t')[1].split('*')[1][0:5]
		DQA1a = DQA1a.split(':')[0]+DQA1a.split(':')[1]
		DQA1b = hlalines[4].split('\t')[2].split('*')[1][0:5]
		DQA1b = DQA1b.split(':')[0]+DQA1b.split(':')[1]
		# DQB1
		DQB1a = hlalines[5].split('\t')[1].split('*')[1][0:5]
		DQB1a = DQB1a.split(':')[0]+DQB1a.split(':')[1]
		DQB1b = hlalines[5].split('\t')[2].split('*')[1][0:5]
		DQB1b = DQB1b.split(':')[0]+DQB1b.split(':')[1]
		# Concatenate DQA/DQB alleles to be in correct format
		DQA1B1a = 'HLA-DQA1'+DQA1a+'-DQB1'+DQB1a
		DQA1B1b = 'HLA-DQA1'+DQA1b+'-DQB1'+DQB1b
		# DRB1
		DRB1a = hlalines[6].split('\t')[1].split('*')[1][0:5]
		DRB1a = DRB1a.split(':')[0]+DRB1a.split(':')[1]
		DRB1b = hlalines[6].split('\t')[2].split('*')[1][0:5]
		DRB1b = DRB1b.split(':')[0]+DRB1b.split(':')[1]
		# Format DRB1 alleles
		DRB1a = 'DRB1_'+DRB1a
		DRB1b = 'DRB1_'+DRB1b
		# Add alleles to list
		hlaalleles.append(DQA1B1a)
		hlaalleles.append(DQA1B1b)
		hlaalleles.append(DRB1a)
		hlaalleles.append(DRB1b)
        hlaalleles = list(set(hlaalleles))  # Remove duplicate alleles if there are any
        hlastring = ','.join(hlaalleles)


	# Run netMHCIIpan
	command = 'export NHOME=/xchip/cga_home/margolis/Packages/netMHCIIPan/netMHCIIpan-3.1; export NETMHCpan=/xchip/cga_home/margolis/Packages/netMHCIIPan/netMHCIIpan-3.1/Linux_x86_64; /xchip/cga_home/margolis/Packages/netMHCIIPan/netMHCIIpan-3.1/netMHCIIpan -a '+hlastring+' -f '+pepfile+' -inptype 0 -length '+str(length)+' -fast -filter 1 -affF 500 -rankF 2.0 -s -xls -xlsfile '+outpath+'/NETMHCIIpan_out_'+str(length)+varianttype+'.xls rdir /xchip/cga_home/margolis/Packages/netMHCIIPan/netMHCIIpan-3.1/Linux_x86_64/ > '+outpath+'/netMHCIIpanoutlen_'+str(length)+varianttype+'.txt'
	subprocess.call(command, shell=True)

	# Catch case where peptide file was empty (create dummy file) 
        dummyfile = outpath+'/NETMHCIIpan_out_'+str(length)+varianttype+'.xls'
        open(dummyfile, 'a').close()	

	return

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
	 # Check to make sure we have the right number of inputs
	if len(sys.argv) != 6:
		print 'Error: incorrect number of inputs.'
		print 'Please input FASTA file(s), a HLAalleles.txt file, the peptide length(s), a netMHCpan version, and an outpath.'
		sys.exit()
	# Parse inputs
	fastas = sys.argv[1]
	alleles = sys.argv[2]
	peplengths = sys.argv[3]
	versionchoice = sys.argv[4]
	outpath = sys.argv[5]
	# Split FASTA files and peptide lengths
	fastalist = fastas.split(',')
	lengthslist = peplengths.split(',')
	if len(fastalist) != len(lengthslist):
		print 'Error: Please make sure your peptide lengths correspond to the fasta files and are in the same order.'
		sys.exit()
	# Run whichever netMHC version is desired
	if versionchoice == '1':
		for i in range(0, len(fastalist)):
			runNetMHCIpan(fastalist[i], alleles, lengthslist[i], outpath)
	else:
		for i in range(0, len(fastalist)):
			runNetMHCIIpan(fastalist[i], alleles, lengthslist[i], outpath)
	
	return

if __name__ == '__main__':
    main()

# ----------------------------------------------------------------------------------------------- #


