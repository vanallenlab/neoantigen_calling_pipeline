# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# mafToFasta.py
# 
# Summary: Takes in a .maf file of either SNVs or InDels and translates mutations to peptide 
# sequences of (desired length*2 - 1). These peptides will be fed into netMHCpan and will result
# in getting the binding affinities of peptide of desired length with at least one AA overlapping
# with the mutated nucleotide(s). Also outputs a map file of headers which will be used to process 
# and annotate the netMHC output.
#
# Input format: python mafToFasta.py maffile maffiletype peptidelengths patientID outpath
# 	maffiletype = 0 for SNVs, 1 for indels
# 	peptidelengths = comma-separated list of lengths (e.g., 9,10 for netMHCpan)
#
# Outputs:
# 	For each length of peptide desired:
#		len#pep_headermap_(snv/indel).txt
#		len#pep_FASTA_(snv/indel).txt
#
# ------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
from Bio.Seq import Seq
from ConfigParser import ConfigParser
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: DNASeqToProtein
# Inputs: List of DNA sequences, corresponding list of headers, peptide length
# Returns: List of protein sequences, corresponding list of headers
# Summary: Translates sequences in-frame to get proteins
def DNASeqToProtein(nucs, headers, length):
        # Initialize codon table and list of peptides
        codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
        peptides = []
        pepheaders = []
        # Translate nucleotide sequences
        for n in range(0, len(nucs)):
                seq = nucs[n]
                fullprotein = ''
                for s in xrange(0, len(seq), 3):
                        codon = seq[s:s+3]
                        # Break if the sequence ends with trailing AAs (should never happen)
                        if len(codon) != 3:
                                break
                        # Find corresponding AA to codon
                        AA = codontable[codon]
                        fullprotein += AA
                # Stop at stop codon, if there is one
                if '*' in fullprotein:
                        substrings = ' '.join(fullprotein.split('*')).split()
                        if len(substrings[0]) >= int(length):
                                peptides.append(substrings[0])
                                pepheaders.append(headers[n])
                else:  # Case when there is no stop codon in full-length protein
                        peptides.append(fullprotein)
                        pepheaders.append(headers[n])
        return peptides, pepheaders
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: MutationsToDNASeq
# Inputs: maffile, peptide length, patient ID, outpath, snv/indel indicator
# Returns: List of mutant/wild type DNA sequences of desired length, list of corresponding headers
# Summary: Takes in maf file, peptide length, patient ID, and outpath. Finds ORF orientation at 
# mutation location from .maf file Codon_change field, based on ORF orientation calculates nucleotide
# window to yield correct number of peptides flanking the mutation, calls twoBitToFa function to 
# get nucleotide sequence, writes header map to output file.
def MutationsToDNASeq(maf, length, patID, outpath, indicator):
	# Read in maf file (desired columns only)
        mafarray = np.loadtxt(maf, dtype=str, delimiter='\t', skiprows=1, usecols=(0,1,4,5,6,8,12,15,35,36,39,40,41,9,38,10), ndmin=2)
        # Create dictionary containing lengths to go backward and forward based on ORF orientation 
        distancedict = {0:[3,0], 1:[2,1], 2:[1,2]}
        # Open header map file for writing
	headermapfile = ''
	# Map file will be labeled differently for SNVs and InDels
	if indicator == 0:
        	headermapfile = open(outpath+'/len'+str(length)+'pep_headermap_snv.txt', 'w')
	else:
		headermapfile = open(outpath+'/len'+str(length)+'pep_headermap_indel.txt', 'a')
	# Translate length from AAs to nucleotides
        length = int(length)*3
	# Initialize sequence list, header list
	seqlist = []
	headerlist = []
	counter = 1
	# Loop through maf array and generate DNA sequences for each sequence
	isnonstop = 0 # Set indicator for nonstop mutations
	nonstopcounter = -1
	nonstopalphabet = 'abcdefghijklmnopqrstuvwxyz'
	for row in mafarray:
		
		# Check to make sure mutation is one that we care about and skip to next sequence if not
		classification = row[5]
                if not (classification == 'Missense_Mutation' or classification == 'Frame_Shift_Ins' or classification == 'Frame_Shift_Del' or classification == 'Nonstop_Mutation' or classification == 'In_Frame_Ins' or classification == 'In_Frame_Del'):
                        continue
		# Go through special case for nonstop mutation
		if (classification == 'Nonstop_Mutation'):
			isnonstop = 1
			nonstopcounter += 1
			nonstopheadermapfile = open(outpath+'/len'+str(length/3)+'pep_headermap_indel.txt','a')
		# Calculate coding strand start and end positions
		orig_start = 0
		orig_end = 0
		if indicator == 0:  # SNVs
			orig_start = int(((row[10].split('.')[1]).split('>')[0])[0:-1])-1 # Subtract because MAFs are 1-indexed but python is 0-indexed
			orig_end = orig_start # SNVs only affect one position
		else:  # InDels
			if row[13] == 'DEL': # deletion
				orig_start = int(((row[10].split('.')[1]).split('del')[0]).split('_')[0])-1
				if '_' in row[10]:
					orig_end = int((row[10].split('del')[0]).split('_')[1])-1
				else:
					orig_end = orig_start
			else: # insertion
				orig_start = int(((row[10].split('.')[1]).split('ins')[0]).split('_')[0])
				if '_' in row[10]:
					orig_end = int((row[10].split('ins')[0]).split('_')[1])-1
				else:
					orig_end = orig_start

		# Calculate mutation length
		mut_length = len(row[15].strip())

                # Calculate ORF orientation at mutation start site
		if indicator == 0:  # If SNV, use codon_change .maf field
                	codonchange = (row[11].split(')')[1]).split('>')[0]
                	orfpos = 0  # Initalize variable
                	for i in range(0, len(codonchange)):
                        	if codonchange[i].isupper():
                                	orfpos = i  # Set ORF variable
		else:  # If InDel, use codon_change .maf field but do further processing
			codonstartnum = (((row[11].split('('))[1]).split('-'))[0]
			cdnastartnum = 0
			if row[13] == 'DEL':
				cdnastartnum = ((((row[10].split('c.'))[1]).split('del'))[0]).split('_')[0]
			else:
				cdnastartnum = ((((row[10].split('c.'))[1]).split('ins'))[0]).split('_')[0]
				cdnastartnum = int(cdnastartnum)+1
			orfpos = (int(cdnastartnum) - int(codonstartnum)) % 3

		# Set new start and end positions for chromosome region with appropriate nucleotide region around mutation site based on ORF orientation
		start = 0
		end = 0
		if indicator == 0:  # For SNVs, do this:
			snvlength = orig_end-orig_start+1
			start = orig_start - (length - distancedict[orfpos][0])
			end = orig_end + (length - distancedict[orfpos][1])
		else:  # For InDels, do this:
			start = orig_start - (length - distancedict[orfpos][0])
			end = orig_start + mut_length + (length - distancedict[orfpos][1])
		
		# Get output from R script that will contain the coding sequence for transcript of interest
		annot_transcript = row[8]
		CONFIG_FILENAME = 'fasta_paths.config'
		config = ConfigParser()
		config.read(CONFIG_FILENAME)
		if isnonstop == 0:
			ref_37_path = config.get('Reference Paths','GRCh37cds')
		else:
			ref_37_path = config.get('Reference Paths', 'GRCh37cdna')
		command = "sed -n -e '/"+annot_transcript+"/,/>/ p' "+ref_37_path+" | sed -e '1d;$d'"
		codingseq = subprocess.check_output(command, shell=True)
		
		# Check to see whether coding sequence was found by biomaRt (if not, continue)
		if len(codingseq) == 0:
			print 'Error: BiomaRt database cannot return coding sequence for transcript '+annot_transcript+'. Skipping this mutation.'
			continue
		
		# Get length of coding sequence plus position of mutation, and get desired sequence start and end indices
		codingseq = codingseq.replace('\n','')
		seqlength = len(codingseq)
		seqstart = 0
		seqend = 0
		if start >= 0:
			seqstart = start
		else:
			seqstart = 0
		if classification == 'Frame_Shift_Ins' or classification == 'Frame_Shift_Del':  # Special case of frameshift mutations
			seqend = seqlength-1
		else:
			if end <= seqlength-1:
				seqend = end
			else:
				seqend = seqlength-1
		
		# Retrieve sequence desired
		sequence = codingseq[seqstart:seqend+1]

		# Substitute in mutation at appropriate position
		disttomut = orig_start - seqstart
		if indicator == 0: # SNV
			mutregion = row[10].split('>')[1]
			mutatedseq = sequence[0:disttomut]+mutregion+sequence[disttomut+snvlength:]
		else:  # InDel
			if row[13] == 'DEL':
				mutatedseq = sequence[0:disttomut]+sequence[disttomut+mut_length:]
			else:
				mutregion = row[10].split('ins')[1]
				mutatedseq = sequence[0:disttomut]+mutregion+sequence[disttomut:]

		# Deal with nonstop mutations in their own way (separately)
		if isnonstop == 1:
			nonstopseqlist = [mutatedseq]
			nonstopheaderlist = ['>seq_'+nonstopalphabet[nonstopcounter]+'_mut']
			nonstoppeptide, nonstoppepheader = DNASeqToProtein(nonstopseqlist, nonstopheaderlist, length/3)		
			nonstopfilehandle = outpath+'/len'+str(length/3)+'pep_FASTA_indel.txt'
			f = open(nonstopfilehandle, 'a')
			f.write(nonstoppepheader[0]+'\n'+nonstoppeptide[0]+'\n')
			nonstopheadermapfile.write('>seq_'+nonstopalphabet[nonstopcounter]+'\t'+patID+'|'+row[7]+'|chr'+row[2]+':'+row[3]+'-'+row[4]+'|'+row[8]+'|'+row[0]+'|'+row[1]+'|'+row[10]+'|'+row[12]+'\n')
			continue

		# Add sequences to lists (mutant and WT for SNV, just mutant for indel)
		seqlist.append(mutatedseq)
		headerlist.append('>seq_'+str(counter)+'_mut')
		if indicator == 0:
			seqlist.append(sequence)
			headerlist.append('>seq_'+str(counter)+'_wt')
		# Write maf annotation information to map file (will be used in netMHC postprocessing)
		headermapfile.write('>seq_'+str(counter)+'\t'+patID+'|'+row[7]+'|chr'+row[2]+':'+row[3]+'-'+row[4]+'|'+row[8]+'|'+row[0]+'|'+row[1]+'|'+row[10]+'|'+row[12]+'\n')
		counter += 1
	headermapfile.close()		

	return seqlist, headerlist
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: writeToOutfile
# Inputs: list of peptides, corresponding headers, current length, outpath, SNV vs. InDel indicator
# Returns: None (writes to file)
# Summary: Writes header+peptide combos to a file, one item per line, that will be an input to netMHC.
# Example output file name: path/len9pep_snv_FASTA.txt (this file would contain SNV peptides of length 9). 
def writeToOutfile(peps, headers, length, outpath, indicator):
	filehandle = ''
	# If SNVs, do this:
	if indicator == 0:
		filehandle = outpath+'/len'+str(length)+'pep_FASTA_snv.txt'
	# If InDels, do this:
	else:
		filehandle = outpath+'/len'+str(length)+'pep_FASTA_indel.txt'
	# Loop through the peptide, header lists and write to filehandle
	f = open(filehandle, 'a')
	for i in range(0, len(peps)):  # The peptide and header lists will always be the same length
		if len(peps[i]) > 0:
			f.write(headers[i]+'\n'+peps[i]+'\n')
	f.close()

	return
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
        # Check to make sure we have the right number of inputs
        if len(sys.argv) != 6:
                print 'Error: incorrect number of inputs.'
                print 'Please input a .maf file, .maf file type, the peptide lengths you want, the patient ID, and an outfile path'
                sys.exit()
        # Store inputs
        maffile = sys.argv[1]
        snvorindel = int(sys.argv[2])
	lengthlist = sys.argv[3].split(',')
        patientID = sys.argv[4]
        outpath = sys.argv[5]
	# For each peptide length in list, do this:
	for length in lengthlist:
		# Convert mutation into nucleotide sequences
		nucleotideseqs, nucleotideheaders = MutationsToDNASeq(maffile, length, patientID, outpath, snvorindel)
		# Convert nucleotide sequences into peptide sequences
		peptideseqs, peptideheaders = DNASeqToProtein(nucleotideseqs, nucleotideheaders, length)
		# Print to outfile
		writeToOutfile(peptideseqs, peptideheaders, length, outpath, snvorindel)

        return

if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------- #
                                                                                                       
