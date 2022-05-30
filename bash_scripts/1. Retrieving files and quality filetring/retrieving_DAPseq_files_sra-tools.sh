#!/usr/bin/bash
conda activate seq

#OBTAINING DAP-seq FILES (date: 29/03/2022)
#source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60143

    mkdir /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data
    cd /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data

    #retrieving DAP-seq files of HB21/40/53 from SRA
    SRAnames=("SRR2926068" "SRR2926398" "SRR2926400" "SRR2926404") #creating a list of the control and the HBs' Run accessions
        # SRR2926068 (GSM1924998): control_none.beads_col_v3a; Arabidopsis thaliana
        # SRR2926398 (GSM1925328): HB_tnt.ATHB21_col_a; Arabidopsis thaliana
        # SRR2926400 (GSM1925330): HB_tnt.ATHB40_col_a; Arabidopsis thaliana
        # SRR2926404 (GSM1925334): HB_tnt.ATHB53_col_a; Arabidopsis thaliana) 
    for sra in ${SRAnames[@]}; do fasterq-dump --split-files $sra; done

    #renaming DAP-seq files
    mv SRR2926068_1.fastq DAP_control_1.fq; mv SRR2926068_2.fastq DAP_control_2.fq; mv SRR2926398.fastq DAP_ATHB21.fq; mv SRR2926400.fastq DAP_ATHB40.fq; mv SRR2926404.fastq DAP_ATHB53.fq