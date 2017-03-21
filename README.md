# neoantigen_calling_pipeline

This pipeline calls somatic cancer neoantigens generated from genetic mutations in patient tumor DNA.

To run: 
  - Download NetMHCPan-3.0 (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan) and configure paths in runNetMHCpan.py file (line 70).
  - Download GRCh37 Ensembl FASTA files, Homo_sapiens.GRCh37.cds.all.fa and Homo_sapiens.GRCh37.cdna.all.fa (http://grch37.ensembl.org/info/data/ftp/index.html) and update paths in fasta_paths.config file.
  - Change paths in shell script getMutNeoantigenBinders.sh (notes in file comments). 
      - For each sample, pipeline is intended to run on a MuTect SNV maf file, a Strelka InDel maf file, and a list of patient Class I HLA alleles (e.g., HLA-A02:01).
  - Run getMutNeoantigenBinders.sh from command line as an SGE Array Job. This script is a wrapper and will call all other relevant scripts.

Additional notes: 
  - Currently, pipeline works for MHC Class I only. MHC Class II functionality in development.
  - Detailed execution instructions and functionality descriptions can be found in each script header, as well as for each individual function.
  
