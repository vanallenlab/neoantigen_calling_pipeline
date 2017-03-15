# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# 17 October 2016
# mutationPostProcess.py
#
# Summary: Takes in NETMHC_out.xls file(s) and does postprocessing to create a more user-friendly 
# output format.
# Input format: python mutationPostProcess.py NETMHCpan_out9snv.xls,NETMHCpan_out9indel.xls,NETMHCpan_out10snv.xls,NETMHCpan_out10indel.txt 
#		len9pep_headermap_snv,len9pep_headermap_indel,len10pep_headermap_snv,len10pep_headermap_indel patientID outpath
# Output format: processedcombinedNETMHCpan_out.txt
#
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
import os
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: processSingleFileOutput
# Inputs: netMHCpan output .xls file (tab-delimited), header map file, patient ID, version
# Returns: An ndarray with rows corresponding to distinct binder peptides and columns corresponding 
# to binder features and metadata.
# Summary: Postprocesses the netMHCpan output to eliminate useless rows, change data format from 
# wide to long, add allele name columns. Also incorporates information about the sequence from the
# .maf file (which is stored in the header map file input parameter).
def processSingleFileOutput(netMHCfile, mapfile, length, patID, version):
	# Catch case where netMHCfile is empty
	if os.path.getsize(netMHCfile) == 0:
		return np.empty(shape=[0, 17])			
	# Parse netMHC filename to get whether file is SNVs or InDels
	snvorindel = 0
	if "InDel" in netMHCfile:
		snvorindel = 1
	length = int(length)
	# Read in first line of file to get number and names of alleles
        with open(netMHCfile, 'r') as f:
        	alleles = f.readline().strip().split('\t')
        alleles = filter(None, alleles)  # Remove empty strings just in case
	# Read in rest of file (skip HLA alleles at the top and file header
	data = np.loadtxt(netMHCfile, dtype=str, delimiter='\t', skiprows=2, ndmin=2)
        nrow = data.shape[0]
        ncol = data.shape[1]
	# Move columns so that data is in long form
        listofarrays = []  # Will store all allele-specific arrays
        initcols = data[:,0:3]  # Initial three columns that are common to all HLA alleles
	if version == 1:
        	for i in range(0, len(alleles)):
                	currstartcol = (3*(i+1))+i
                	currendcol = currstartcol+4
                	currarray = data[:,currstartcol:currendcol]
                	listofarrays.append(currarray)
        		datav2 = np.vstack(tuple(listofarrays))
	else:
		for i in range(0, len(alleles)):
			currstartcol = (3*(i+1))
			currendcol = currstartcol+4
			currarray = data[:,currstartcol:currendcol]
			listofarrays.append(currarray)
			datav2 = np.vstack(tuple(listofarrays))
        # Add initial columns and allele column into data frame
        # Allele column
        allelevec = []
        for i in range(0, len(alleles)):
                currnewcol = [alleles[i]]*nrow
                allelevec.extend(currnewcol)
        datav2 = np.insert(datav2, 1, allelevec, axis=1)  # Add allele column to datalong
        # Initial columns
        initcollist = []
        for i in range(0, len(listofarrays)):
                initcollist.append(initcols)
        initcolstoappend = np.vstack(tuple(initcollist))
        datav3 = np.concatenate((initcolstoappend, datav2), axis=1)
	datav3 = np.delete(datav3, 3, axis=1)
	# Create mutant / WT dictionary for SNVs
	if snvorindel == 0:
		newnrow = datav3.shape[0]
		newncol = datav3.shape[1]
		mutWTdict = {}
		WTindices = []
		for i in range(0, newnrow):
			if datav3[i,2].strip().split('_')[2] == 'mut':  # For each mutant row, do:
				currkey = datav3[i,2]+'|'+datav3[i,0]+'|'+datav3[i,3]
				# Find the corresponding wild-type row
				for j in range(i+1, newnrow):
					if datav3[j,0] == datav3[i,0] and datav3[j,2][0:5] == datav3[i,2][0:5]:
						currval = datav3[j,2]+'|'+datav3[j,1]+'|'+datav3[j,4]+'|'+datav3[j,5]+'|'+datav3[j,6]
						mutWTdict[currkey] = currval
						WTindices.append(j)
						break
		# Delete WT rows
		datav4 = np.delete(datav3, WTindices, axis=0)
	else:
		datav4 = datav3
	# Eliminate any rows that have a rank above 2%
	toremove = []
	newnrow2 = datav4.shape[0]
	for i in range(0, newnrow2):
		if float(datav4[i,6]) > 2:
			toremove.append(i)
	datav5 = np.delete(datav4, toremove, 0)	
	# Read in map file and create map dictionary
	headerdict = {}
	with open(mapfile, 'r') as f:
	 	lines = f.read().splitlines()
	for line in lines:
		key = line.split('\t')[0][1:]
		val = line.split('\t')[1]
		headerdict[key] = val
	# Initialize new ndarray columns (will eventually use np.hstack to stack them all theogether into a numpy array)
	patID,sample,transcript,chrom_loc,gene,gene_num,cdna_change,prot_change,pep_pos,pep_length,hla,pep_mut,aff_mut,rank_mut,pep_wt,aff_wt,rank_wt = ([] for i in range(17))
	# For every row in current data array, use WT dict and header dict to find metainformation and save in a new ndarray
	newnrow3 = datav5.shape[0]
	for i in range(0, newnrow3):
		currrow = datav5[i,:]
		seqnum = currrow[2].split('_')
		headerdictkey = seqnum[0]+'_'+seqnum[1]
		headervals = headerdict[headerdictkey].split('|')
		patID.append(headervals[0])
		sample.append(headervals[1])
		transcript.append(headervals[3])
		chrom_loc.append(headervals[2])
		gene.append(headervals[4])
		gene_num.append(headervals[5])
		cdna_change.append(headervals[6])
		prot_change.append(headervals[7])
		pep_pos.append(currrow[0])
		pep_length.append(length)
		hla.append(currrow[3])
		pep_mut.append(currrow[1])
		aff_mut.append(currrow[5])
		rank_mut.append(currrow[6])
		if snvorindel == 0:
			WTdictkey = currrow[2]+'|'+currrow[0]+'|'+currrow[3]
			WTvals = mutWTdict[WTdictkey].split('|')
			pep_wt.append(WTvals[1])
			aff_wt.append(WTvals[3])
			rank_wt.append(WTvals[4])
		else:
			pep_wt.append('n/a')
			aff_wt.append('n/a')
			rank_wt.append('n/a')
	# Join all lists into new numpy array
	datafull = np.column_stack((patID,sample,transcript,chrom_loc,gene,gene_num,cdna_change,prot_change,pep_pos,pep_length,hla,pep_mut,aff_mut,rank_mut,pep_wt,aff_wt,rank_wt))

	return datafull
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: writeToOutfile
# Inputs: final numpy array, patient ID, outpath
# Returns: None (writes to file)
# Summary: Takes the full numpy array and writes it, plus appropriate header, to a tab-delimited
# file in the specified outpath.
def writeToOutfile(array, patID, version, outpath):
	suffix = ''
	if version == 1:
		suffix = '_processedcombinedNETMHCpan_out.txt'
	else:
		suffix = '_processedcombinedNETMHCIIpan_out.txt'
	outfile = outpath+'/'+patID+suffix
	headerstring = 'patient\tsample\ttranscript\tchrom_loc\tgene\tgene_num\tcdna_change\tprot_change\tmut_pos\tpep_length\tHLA\tpep_mut\taff_mut\trank_mut\tpep_wt\taff_wt\trank_wt'
	np.savetxt(outfile, array, fmt='%s', delimiter='\t', header = headerstring, comments = '')

	return

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
         # Check to make sure we have the right number of inputs
        if len(sys.argv) != 7:
                print 'Error: incorrect number of inputs.'
                print 'Please input netMHC outfile(s), corresponding header map(s), corresponding length(s), the patient ID, version of netMHCpan that was run, and an outpath.'
                sys.exit()
	# Parse inputs
	netmhcoutfiles = sys.argv[1]
	headermapfiles = sys.argv[2]
	lengths = sys.argv[3]
	patientID = sys.argv[4]
	version = int(sys.argv[5])
	outfilepath = sys.argv[6]
	# Split FASTA files and peptide lengths, making sure there are no leading/trailing whitespaces
	netmhclist = list(map(str.strip, netmhcoutfiles.split(',')))
	headerlist = list(map(str.strip, headermapfiles.split(',')))
	lengthslist = list(map(str.strip, lengths.split(',')))
        if len(netmhclist) != len(headerlist):
                print 'Error: Please make sure your header map files correspond to the netMHC outfiles and are in the same order.'
                sys.exit()
  	# Loop through each netMHC file and add in relevant information
	procarrays = []
	for i in range(0, len(netmhclist)):
		curroutputprocessed = processSingleFileOutput(netmhclist[i], headerlist[i], lengthslist[i], patientID, version)
		procarrays.append(curroutputprocessed)
	# If there is more than one processed array, concatenate them together then write to outfile
	if len(procarrays) > 1:
		fullarray = np.concatenate(tuple(procarrays), axis=0)
		writeToOutfile(fullarray, patientID, version, outfilepath)
	else:  # Otherwise, write single array to file
		writeToOutfile(procarrays[0], patientID, version, outfilepath)

	return

if __name__ == '__main__':
        main()
# ----------------------------------------------------------------------------------------------- #
