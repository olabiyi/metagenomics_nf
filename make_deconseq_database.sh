#!/usr/bin/env bash

# Sample script to Download maize (HOST) and 
# Phix genome for contaminant removal using deconseq

# Pull the required docker images
docker pull biocontainers/entrez-direct:v7.50.20171103_cv5
docker pull dceoy/prinseq:latest


# Go to NCBI to find the URL for the host genome
# Here we fetch data for maize genome
URL=' https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_902167145.1/download?include_annotation_type=GENOME_FASTA'
wget -O genome.zip ${URL}
unzip genome.zip

# Download Phix genome from here https://www.ncbi.nlm.nih.gov/search/all/?term=NC_001422.1
HOST='ncbi_dataset/data/GCF_902167145.1/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna'
PHIX='Phix_NC_001422.1.fna'
# Concating host genome with Phix to remove host contaminating DNA
# and splitting sequences by long repeats of ambiguous base N to aviod fasle positives
seqkit concat --full ${PHIX} ${HOST}  | \
   perl -p -e 's/N\n/N/' | \
   perl -p -e 's/^N+//;s/N+$//;s/N{200,}/\n>split\n/' > host.fna

mkdir database/

mv ${PHIX} ${HOST} host.fna database/


# Tutorials on how to build a Deconseq database can be found by following the links below
#https://deconseq.sourceforge.net/manual.html#DB
#https://sourceforge.net/projects/deconseq/files/


# Filter sequences
docker run  --rm -v ${PWD}:${PWD} -w ${PWD} -u $(id -u):$(id -g) dceoy/prinseq:latest perl prinseq-lite.pl -log -verbose -fasta database/host.fna -min_len 200 -ns_max_p 10 -derep 12345 -out_good database/host_prinseq -seq_id maize_genome -rm_header -out_bad null
# Build database
docker run  --rm -v ${PWD}:${PWD} -w ${PWD} -u $(id -u):$(id -g) quay.io/grbot/deconseq:latest bwa index -p database/host_prinseqDB -a bwtsw database/host_prinseq.fasta &