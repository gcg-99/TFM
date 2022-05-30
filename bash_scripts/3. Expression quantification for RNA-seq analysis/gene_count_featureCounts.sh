#!/usr/bin/bash
conda activate gene_count
SE_bam_files=(/home/gema/CNB/RNAseq/2022-03-15_vanEsEtal_sequencing_data_mapping/*.sorted.bam)
PE_bam_files=(/home/gema/CNB/RNAseq/2022-03-15_Novogene_sequencing_data_mapping/*.sorted.bam)
ABA_experiment=${PE_bam_files[@]:0:18}
HBs_experiment=${PE_bam_files[@]:18:30}


#EXPRESSION QUANTIFICATION FOR DIFFERENTIAL EXPRESSION ANALYSIS

    #featureCounts is part of the Subread package
        ##   -T <int> -> Number of the threads (1 by default)
        ##   -s <int> -> Perform strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
        ##   -p -> If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads (htseq-count does this by default)
        ##   -O -> Assign reads to all their overlapping meta-features***
        ##   -C -> Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands
        ##   -t <string> -> Specify feature type in GTF annotation ('exon' by default)
        ##   -g <string> -> Specify attribute type in GTF annotation ('gene_id' by default)
        ##   -a <string> -> Name of an annotation file (GTF/GFF format by default)
        ##   -o <string> -> Name of the output file including read counts. A separate file including summary statistics of counting results is also included in the output ('<string>.summary')


    ##SINGLE-END MAPPED READS COUNT
    ## ignoring overlapping genes
    cd /home/gema/CNB/RNAseq/2022-03-24_vanEsEtal_sequencing_data_counts
    echo "Counting how many reads of single-end bam files map to each genomic feature with featureCounts"
    featureCounts -T 8 -t gene -g gene_id -a /home/gema/CNB/reference_genome/Arabidopsis_thaliana.TAIR10.52.gtf -o BRC1ind_nooverlapgen_featureCounts.count ${SE_bam_files[@]}

    ## including overlapping meta-features***
    echo "Counting how many reads of single-end bam files map to each genomic feature with featureCounts"
    featureCounts -T 8 -O -t gene -g gene_id -a /home/gema/CNB/reference_genome/Arabidopsis_thaliana.TAIR10.52.gtf -o BRC1ind_overlapgen_featureCounts.count ${SE_bam_files[@]}



    ##PAIR-END MAPPED READS COUNT
    ## ignoring overlapping genes
    cd /home/gema/CNB/RNAseq/2022-03-24_Novogene_sequencing_data_counts
    echo "Counting how many reads of ABA experiment bam files map to each genomic feature with featureCounts"
    featureCounts -T 8 -s 2 -p -C -t gene -g gene_id -a /home/gema/CNB/reference_genome/Arabidopsis_thaliana.TAIR10.52.gtf -o ABA_experiment_nooverlapgen_featureCounts.count ${ABA_experiment}
    echo "Counting how many reads of HBs inducible lines bam files map to each genomic feature with featureCounts"
    featureCounts -T 8 -s 2 -p -C -t gene -g gene_id -a /home/gema/CNB/reference_genome/Arabidopsis_thaliana.TAIR10.52.gtf -o HBs_experiment_nooverlapgen_featureCounts.count ${HBs_experiment}

    ## including overlapping meta-features***
    echo "Counting how many reads of ABA experiment bam files map to each genomic feature with featureCounts"
    featureCounts -T 8 -s 2 -p -O -C -t gene -g gene_id -a /home/gema/CNB/reference_genome/Arabidopsis_thaliana.TAIR10.52.gtf -o ABA_experiment_overlapgen_featureCounts.count ${ABA_experiment}
    echo "Counting how many reads of HBs inducible lines bam files map to each genomic feature with featureCounts"
    featureCounts -T 8 -s 2 -p -O -C -t gene -g gene_id -a /home/gema/CNB/reference_genome/Arabidopsis_thaliana.TAIR10.52.gtf -o HBs_experiment_overlapgen_featureCounts.count ${HBs_experiment}