#!/usr/bin/bash
conda activate peak_calling
mkdir /home/gema/CNB/ChIPseq/2022-03-30_vanEsEtal_sequencing_data_peak-calling
mkdir /home/gema/CNB/DAPseq/2022-03-30_OMalleyEtal_sequencing_data_peak-calling


#PEAK CALLING

    ## ChIP-seq FILES
    cd /home/gema/CNB/ChIPseq/2022-03-30_vanEsEtal_sequencing_data_peak-calling

    macs2 callpeak -t /home/gema/CNB/ChIPseq/2022-03-15_vanEsEtal_sequencing_data_mapping/ChIP2.sorted.bam -c /home/gema/CNB/ChIPseq/2022-03-15_vanEsEtal_sequencing_data_mapping/Input2.sorted.bam -g 119481543 --outdir macs2 -q 1 -n "ChIP2_q1" --nomodel --extsize 206 --call-summits > ChIP2_q1-macs2-alg-stdout+stderr.txt 2>&1
    macs2 callpeak -t /home/gema/CNB/ChIPseq/2022-03-15_vanEsEtal_sequencing_data_mapping/ChIP3.sorted.bam -c /home/gema/CNB/ChIPseq/2022-03-15_vanEsEtal_sequencing_data_mapping/Input3.sorted.bam -g 119481543 --outdir macs2 -q 1 -n "ChIP3_q1" --nomodel --extsize 206 --call-summits > ChIP3_q1-macs2-alg-stdout+stderr.txt 2>&1
    macs2 callpeak -t /home/gema/CNB/ChIPseq/2022-03-15_vanEsEtal_sequencing_data_mapping/ChIP7.sorted.bam -c /home/gema/CNB/ChIPseq/2022-03-15_vanEsEtal_sequencing_data_mapping/Input7.sorted.bam -g 119481543 --outdir macs2 -q 1 -n "ChIP7_q1" --nomodel --extsize 206 --call-summits > ChIP7_q1-macs2-alg-stdout+stderr.txt 2>&1



    ##DAP-seq FILES
    cd /home/gema/CNB/DAPseq/2022-03-30_OMalleyEtal_sequencing_data_peak-calling

    macs2 callpeak -t /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping/DAP_ATHB21.sorted.bam -c /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping/DAP_control.sorted.bam -g 119481543 --outdir macs2 -q 1 -n "DAP_ATHB21_q1_es200" --nomodel --extsize 200 --call-summits > DAP_ATHB21_q1_es200-macs2-alg-stdout+stderr.txt 2>&1
    macs2 callpeak -t /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping/DAP_ATHB40.sorted.bam -c /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping/DAP_control.sorted.bam -g 119481543 --outdir macs2 -q 1 -n "DAP_ATHB40_q1_es200" --nomodel --extsize 200 --call-summits > DAP_ATHB40_q1_es200-macs2-alg-stdout+stderr.txt 2>&1
    macs2 callpeak -t /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping/DAP_ATHB53.sorted.bam -c /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping/DAP_control.sorted.bam -g 119481543 --outdir macs2 -q 1 -n "DAP_ATHB53_q1_es200" --nomodel --extsize 200 --call-summits > DAP_ATHB53_q1_es200-macs2-alg-stdout+stderr.txt 2>&1