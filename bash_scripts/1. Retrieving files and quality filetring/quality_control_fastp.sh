#!/usr/bin/bash
conda activate seq

#fastp (https://github.com/OpenGene/fastp)
##    -q, --qualified_quality_phred -> Default 15 means phred quality >=Q15 is qualified

    ##SINGLE-END DATA
    #GFP:BRC1ind experiment files
    mkdir /home/gema/CNB/RNAseq/2022-03-14_vanEsEtal_sequencing_data_fastp
    cd /home/gema/CNB/RNAseq/2022-03-14_vanEsEtal_sequencing_data_fastp
    RNAfileslist=($(ls /home/gema/CNB/RNAseq/2020_vanEsEtal_single-end_sequencing_data | grep -E 'fq.gz|fastq.gz' | cut -d. -f1)) #storing the name of each file of readings (e.g., GFP-BRC1-2) in an array
    for file in ${RNAfileslist[@]}; do echo "Processing $file with fastp";
    fastp -i /home/gema/CNB/RNAseq/2020_vanEsEtal_single-end_sequencing_data/$file.fastq.gz -o ${file%%_*}_fastp.fq.gz --thread 8 --dont_overwrite -q 20 --disable_adapter_trimming --overrepresentation_analysis --failed_out ${file%%_*}_failed_reads.fq --html $file.fastp-report.html --json $file.fastp-report.json --report_title "fastp report of $file" 2> $file.fastp-report.txt;
    done

    #ChIP-seq experiment files
    mkdir /home/gema/CNB/ChIPseq/2022-03-14_vanEsEtal_sequencing_data_fastp
    cd /home/gema/CNB/ChIPseq/2022-03-14_vanEsEtal_sequencing_data_fastp
    ChIPfileslist=($(ls /home/gema/CNB/ChIPseq/2020_vanEsEtal_single-end_sequencing_data | grep -E 'fq.gz|fastq.gz' | cut -d. -f1)) #storing the name of each file of readings (e.g., ChIP2_S5_R1_001) in an array
    for file in ${ChIPfileslist[@]}; do echo "Processing $file with fastp";
    fastp -i /home/gema/CNB/ChIPseq/2020_vanEsEtal_single-end_sequencing_data/$file.fastq.gz -o ${file%%_*}_fastp.fq.gz --thread 8 --dont_overwrite -q 20 --disable_adapter_trimming --overrepresentation_analysis --failed_out ${file%%_*}_failed_reads.fq --html $file.fastp-report.html --json $file.fastp-report.json --report_title "fastp report of $file" 2> $file.fastp-report.txt;
    done

    #DAP-seq experiment files (O'Malley et al., 2016)
    mkdir /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_fastp
    cd /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_fastp
    SE_DAPfileslist=($(ls /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data | grep -E 'DAP_ATHB' | cut -d. -f1)) #storing the name of each file of readings (e.g., DAP_ATHB21) in an array
    for file in ${SE_DAPfileslist[@]}; do echo "Processing $file with fastp";
    fastp -i /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data/$file.fq -o ${file}_fastp.fq --thread 8 --dont_overwrite -q 20 --disable_adapter_trimming --overrepresentation_analysis --failed_out ${file}_failed_reads.fq --html $file.fastp-report.html --json $file.fastp-report.json --report_title "fastp report of $file" 2> $file.fastp-report.txt;
    done


    ##PAIR-END DATA
    #ABA and HBs experiment files
    mkdir /home/gema/CNB/RNAseq/2022-03-14_Novogene_sequencing_data_fastp
    cd /home/gema/CNB/RNAseq/2022-03-14_Novogene_sequencing_data_fastp
    filenames=($(ls /home/gema/CNB/RNAseq/2022-02-25_Novogene_sequencing_data/ | grep '_1' | cut -c-3)) #storing the ID of each pair of readings (e.g., A01) in an array
    for file in ${filenames[@]}; do echo "Processing ${file%%_*}_1.fq.gz and ${file%%_*}_2.fq.gz with fastp";
    fastp -i /home/gema/CNB/RNAseq/2022-02-25_Novogene_sequencing_data/${file%%_*}_1.fq.gz -I /home/gema/CNB/RNAseq/2022-02-25_Novogene_sequencing_data/${file%%_*}_2.fq.gz -o ${file%%_*}_1_fastp.fq.gz -O ${file%%_*}_2_fastp.fq.gz --thread 8 --dont_overwrite -q 20 --disable_adapter_trimming --overrepresentation_analysis --failed_out ${file%%_*}_failed_reads.fq --html $file.fastp-report.html --json $file.fastp-report.json --report_title "fastp report of $file" 2> $file.fastp-report.txt;
    done

    #DAP-seq control files (O'Malley et al., 2016)
    fastp -i /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data/DAP_control_1.fq -I /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data/DAP_control_2.fq -o DAP_control_1_fastp.fq -O DAP_control_2_fastp.fq --thread 8 --dont_overwrite -q 20 --disable_adapter_trimming --overrepresentation_analysis --failed_out DAP_control_failed_reads.fq --html DAP_control.fastp-report.html --json DAP_control.fastp-report.json --report_title "fastp report of DAP_control" 2> DAP_control.fastp-report.txt