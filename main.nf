#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


log.info """\
 M E T A G E N O M I C S - N F   P I P E L I N E
 ===================================
 Host Genome : ${params.genomeDB}
 Reads File       : ${params.csv_file}
 PFAM Reference       : ${params.reference}
 """

process Trim_reads {

    tag "Trimming the forward and reverse reads of ${sample} using trim_galore...."
    //publishDir "02.Trim_reads/" , mode: "copy"

    input:
        tuple val(sample), path(forward), path(rev)
    output: 
        tuple val(sample), path("${sample}_R1_val_1.fq"), path("${sample}_R2_val_2.fq")
    script:
        """
        trim_galore --paired --length 20 --dont_gzip  -o ./ ${forward} ${rev}
	"""
}

// Removing sequences that look contaminated using deconseq from the forward reads
process Remove_contaminants_Forward {

    tag "Removing contaminants from ${sample}'s forward reads"
    label 'Deconseq'
    //publishDir "03.Remove_contaminants/" , mode: "copy"
    
    input:
        each path(db)
        each path(bwa64)
        tuple val(sample), path(forward), path(rev)
    output: 
        tuple val(sample), path("${sample}_for_read_clean.fq") 
    script:
        """  
        [ -d database/ ] || mkdir database/ && mv ${db}  database/
        USER=\$(whoami)
        chown \${USER}:\${USER} ./bwa64 && chmod +x ./bwa64
        deconseq.pl -f ${forward} -dbs plant -out_dir ./ -id ${sample}_for_read
	"""
}

// Removing sequences that look contaminated using deconseq from the reverse reads
process Remove_contaminants_Rev {

    tag "Removing contaminants from ${sample}'s reverse reads"
    label 'Deconseq'
    //publishDir "03.Remove_contaminants/" , mode: "copy"

    input:
        each path(db)
        each path(bwa64)
        tuple val(sample), path(forward), path(rev)
    output: 
        tuple val(sample), path("${sample}_rev_read_clean.fq") 
    script:
        """  
        [ -d database/ ] || mkdir database/ && mv ${db}  database/
        USER=\$(whoami)
        chown \${USER}:\${USER} ./bwa64 && chmod +x ./bwa64
        deconseq.pl -f ${rev} -dbs plant -out_dir ./ -id ${sample}_rev_read 
	"""
}

// Rewriting paired end fastq files to make sure that all reads 
// have a mate and to separate out singletons

 process Rewrite_pairs {

    tag "Rewriting pairs for ${sample}"
    //publishDir "04.Rewrite_pairs/" , mode: "copy"

    input:
        tuple val(sample), path(forward), path(rev)
    output: 
        tuple val(sample), path("${sample}_for_read_clean.fq.paired.fq"), path("${sample}_rev_read_clean.fq.paired.fq") 
    script:
        """ 
        fastq_pair ${forward} ${rev}
        """
 }


// Assembly

process Assemble_metagenome {

    tag "Assembling the metagenome for ${sample} using megahit"
    //publishDir "05.Assemble_metagenome/" , mode: "copy"
       
    input:
        tuple val(sample), path(forward), path(rev)
    output: 
        tuple val(sample), path("megahit/${sample}.contigs.fa")
    script:
        """
        megahit \
           -1 ${forward} \
           -2 ${rev} \
	   -t ${task.cpus} \
	   -o megahit/ \
	   --out-prefix ${sample}
	"""
}


// Annotation
// Running Prodigal to find open reading frames (ORF)
process Predict_protiens {

    tag "Finding open reading frames for ${sample}"
    //publishDir "06.Predict_proteins/" , mode: "copy", pattern: "*_proteins*" 
    
    input:
        tuple val(sample), path(assembly)
    output:
        tuple val(sample), path("${sample}_proteins.faa")
    script:
        """
        prodigal \
           -p meta \
	   -f gff \
	   -i ${assembly} \
	   -o "${sample}_proteins.gff" \
	   -a "${sample}_proteins.faa" \
	   -d "${sample}_proteins.fna"
        """
}


process Protein_families {

    tag "Running HMMScan (hmm models) on predicted proteins for ${sample}"
    publishDir "07.Protein_families/" , mode: "copy"

    input:
        each path(reference)
        tuple val(sample), path(proteins)

    output:
        tuple val(sample), path("${sample}_hmmsearch.PFAM.tsv")

    script:
        """
         hmmsearch \
           --domE 1e-5 \
	   --domtblout "${sample}_hmmsearch.PFAM.tsv" \
           --noali --cpu ${task.cpus} ${reference} ${proteins}
	"""
}

workflow {
    // the .collect changes a queue channel to a value channel 
    // thus allowing the channel to be consumed multiple times
    ref_ch = Channel.fromPath(params.reference).collect()
    bwa64_ch = Channel.fromPath(params.bwa64).collect()
    genome_ch =  Channel.fromPath(params.genomeDB, type: 'dir').collect()
    
    
    Channel.fromPath( params.csv_file )
        .splitCsv()
        .map{row -> tuple( "${row[0]}", file("${row[1]}"), file("${row[2]}") )}
        .set{reads_ch}        

    trim_ch = Trim_reads(reads_ch)     

    forward_ch = Remove_contaminants_Forward(genome_ch, bwa64_ch, trim_ch) 
    rev_ch = Remove_contaminants_Rev(genome_ch, bwa64_ch, trim_ch)
    
    clean_ch = forward_ch.join(rev_ch)

    protein_ch =  clean_ch | 
                  Rewrite_pairs |
	          Assemble_metagenome |
	          Predict_protiens 
    
    result_ch = Protein_families(ref_ch, protein_ch)   
}
