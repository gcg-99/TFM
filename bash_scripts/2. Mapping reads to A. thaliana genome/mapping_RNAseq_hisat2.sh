#!/usr/bin/bash
conda activate seq

#READ MAPPING
    
    ##SINGLE-END DATA ALIGNMENT
    mkdir /home/gema/CNB/RNAseq/2022-03-15_vanEsEtal_sequencing_data_mapping
    cd /home/gema/CNB/RNAseq/2022-03-15_vanEsEtal_sequencing_data_mapping
    SEfilenames=($(ls /home/gema/CNB/RNAseq/2022-03-14_vanEsEtal_sequencing_data_fastp | grep -E 'fq.gz|fastq.gz' | cut -d. -f1 | awk '{ print substr( $0, 1, length($0)-6 ) }')) #retrieving the name of each file of readings without the fastp extension (i.e., GFP-BRC1-2 instead of GFP-BRC1-2_fastp) and storing them in an array
    for file in ${SEfilenames[@]}; do echo "Mapping $file to the reference genome with hisat2";
    hisat2 -p 8 --very-sensitive --max-intronlen 11700 -t -x /home/gema/CNB/reference_genome/hisat2_athaliana_genome -U /home/gema/CNB/RNAseq/2022-03-14_vanEsEtal_sequencing_data_fastp/${file%%_*}_fastp.fq.gz -S ${file%%_*}_hisat2.sam > $file-hisat2-alg-stdout+stderr.txt 2>&1;
    echo "Converting ${file%%_*}_hisat2.sam to $file.bam";
    samtools view -b -h -@ 8 ${file%%_*}_hisat2.sam > $file.bam; #samtools view (-b is to generate the output in BAM format, -h to include the header in the output, and -@ to allocate additional threads)
    echo "Showing the header of $file.bam";
    samtools view -H -@ 8 $file.bam; #samtools view (-H is to output only the header)
    echo "Providing simple statistics of $file.bam";
    samtools flagstat -@ 8 $file.bam > samstats-$file.err; #providing simple statistics on BAM files
    echo "Sorting $file.bam to $file.sorted.bam";
    samtools sort -@ 8 $file.bam > $file.sorted.bam; #sorting the mappings
    echo "Indexing $file.sorted.bam";
    samtools index -@ 8 $file.sorted.bam; #indexing sorted bam files
    echo "Converting $file.sorted.bam to ${file%%_*}_RPKM.bw";
    bamCoverage -p 8 --binSize 5 --effectiveGenomeSize 119481543 --normalizeUsing RPKM -b $file.sorted.bam -o ${file%%_*}_RPKM.bw;
        #bamCoverage takes an alignment of reads or fragments as input (BAM file) and generates a coverage track (bigWig or bedGraph) as output
            ##   --numberOfProcessors, -p -> Number of processors to use
            ##   --binSize, -bs -> Size of the bins, in bases, for the output of the bigwig/bedgraph file
            ##   --effectiveGenomeSize -> The effective genome size is the portion of the genome that is mappable
    done


    ##PAIR-END DATA ALIGNMENT
    mkdir /home/gema/CNB/RNAseq/2022-03-15_Novogene_sequencing_data_mapping
    cd /home/gema/CNB/RNAseq/2022-03-15_Novogene_sequencing_data_mapping
    PEfilenames=($(ls /home/gema/CNB/RNAseq/2022-03-14_Novogene_sequencing_data_fastp/ | grep '_1' | cut -c-3)) #storing the ID of each pair of readings (e.g., A01) in an array
    for file in ${PEfilenames[@]}; do echo "Mapping $file to the reference genome with hisat2";
    hisat2 -p 8 --very-sensitive --rna-strandness RF --max-intronlen 11700 -t -x /home/gema/CNB/reference_genome/hisat2_athaliana_genome -1 /home/gema/CNB/RNAseq/2022-03-14_Novogene_sequencing_data_fastp/${file%%_*}_1_fastp.fq.gz -2 /home/gema/CNB/RNAseq/2022-03-14_Novogene_sequencing_data_fastp/${file%%_*}_2_fastp.fq.gz -S ${file%%_*}_hisat2.sam > $file-hisat2-alg-stdout+stderr.txt 2>&1;
    echo "Converting ${file%%_*}_hisat2.sam to $file.bam";
    samtools view -b -h -@ 8 ${file%%_*}_hisat2.sam > $file.bam; #samtools view (-b is to generate the output in BAM format, -h to include the header in the output, and -@ to allocate additional threads)
    echo "Showing the header of $file.bam";
    samtools view -H -@ 8 $file.bam; #samtools view (-H is to output only the header)
    echo "Providing simple statistics of $file.bam";
    samtools flagstat -@ 8 $file.bam > samstats-$file.err; #providing simple statistics on BAM files
    echo "Sorting $file.bam to $file.sorted.bam";
    samtools sort -@ 8 $file.bam > $file.sorted.bam; #sorting the mappings
    echo "Indexing $file.sorted.bam";
    samtools index -@ 8 $file.sorted.bam; #indexing sorted bam files
    echo "Converting $file.sorted.bam to ${file%%_*}_RPKM.bw";
    bamCoverage -p 8 --binSize 5 --effectiveGenomeSize 119481543 --normalizeUsing RPKM -b $file.sorted.bam -o ${file%%_*}_RPKM.bw;
    done