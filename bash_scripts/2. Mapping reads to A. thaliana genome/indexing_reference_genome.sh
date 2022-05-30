#!/usr/bin/bash

#CREATING THE REFERENCE GENOME (date: 15/03/2022)
    mkdir /home/gema/CNB/reference_genome
    cd /home/gema/CNB/reference_genome

    #retrieving Arabidopsiss thaliana genome files and the gff file from ensembl
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa.gz 
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa.gz 
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa.gz 
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa.gz 
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa.gz
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.Mt.fa.gz
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.Pt.fa.gz
    wget http://ftp.ensemblgenomes.org/pub/plants/release-52/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gff3.gz
    gunzip *.gz #uncompressing 'fasta' and 'gff' files

    #merging TAIR10 chromosome files into a single file with the reference genome
    cat Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa Arabidopsis_thaliana.TAIR10.dna.chromosome.Mt.fa Arabidopsis_thaliana.TAIR10.dna.chromosome.Pt.fa > genome.fasta
    grep '>' genome.fasta #checking that the new file contains all the chromosomes



#INDEX BUILDING OF THE REFERENCE GENOME

    #HISAT2 INDEXING
    #hisat2 scripts require a gtf file
    conda activate agat
    agat_convert_sp_gff2gtf.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3 -o Arabidopsis_thaliana.TAIR10.52.gtf #converting gff3 to gtf
    
    #creating a lists of known splice sites, known exons and intron size for hisat2
    agat_sp_manage_introns.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3 > Arabidopsis_thaliana.TAIR10.52_exon_intron_info.txt
        #longest intron = 11602 bp

    conda activate seq        
    hisat2_extract_splice_sites.py Arabidopsis_thaliana.TAIR10.52.gtf > Arabidopsis_thaliana.TAIR10.52_splicesites_hisat2.txt
    hisat2_extract_exons.py Arabidopsis_thaliana.TAIR10.52.gtf > Arabidopsis_thaliana.TAIR10.52_exons_hisat2.txt

    #hisat2 index building
    fastafileslist=$(find *.fa | paste -sd,)
    hisat2-build --threads 8 --ss Arabidopsis_thaliana.TAIR10.52_splicesites_hisat2.txt --exon Arabidopsis_thaliana.TAIR10.52_exons_hisat2.txt -f genome.fasta hisat2_athaliana_genome > hisat2-index-build-stdout+stderr.txt 2>&1


    #BOWTIE2 INDEXING
    fastafileslist=$(find *.fa | paste -sd,)
    bowtie2-build --threads 8 -f genome.fasta bowtie2_athaliana_genome > bowtie2-index-build-stdout+stderr.txt 2>&1 #building a small index