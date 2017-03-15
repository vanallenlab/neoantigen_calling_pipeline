#!/usr/bin/env Rscript

library('biomaRt')

# Get transcript ID passed from python script
functinput <- commandArgs(trailingOnly = TRUE)

# Function to grab coding sequence
getCodingSeq <- function(transcriptID) {
	# Set up mart database
	ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org")
	ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
	# Do getSequence query
	transcriptseq <- paste(getSequence(id=transcriptID, type="ensembl_transcript_id", seqType="coding", mart=ensembl)$coding,
	      getSequence(id=transcriptID, type="ensembl_transcript_id", seqType="3utr", mart=ensembl)$`3utr`, sep="")
	# Return to command line
	cat(transcriptseq)
}

# Call function
getCodingSeq(functinput)
