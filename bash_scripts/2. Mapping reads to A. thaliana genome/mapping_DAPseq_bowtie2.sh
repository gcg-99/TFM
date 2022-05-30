#!/usr/bin/bash
conda activate seq

#READ MAPPING
    
    ##SINGLE-END DATA ALIGNMENT
    mkdir /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping
    cd /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_mapping
    DAPfileslist=($(ls /home/gema/CNB/DAPseq/2016_OMalleyEtal_sequencing_data | grep -E 'DAP_ATHB' | cut -d. -f1)) #storing the name of each file of readings (e.g., DAP_ATHB21) in an array
    for file in ${DAPfileslist[@]}; do echo "Mapping $file to the reference genome with bowtie2";
    bowtie2 -p 8 --very-sensitive-local -t -x /home/gema/CNB/reference_genome/bowtie2_athaliana_genome -U /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_fastp/${file}_fastp.fq -S ${file}_bowtie2.sam > $file-bowtie2-alg-stdout+stderr.txt 2>&1;
        ##--very-sensitive-local -> Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    echo "Converting ${file}_bowtie2.sam to $file.bam";
    samtools view -b -h -@ 8 ${file}_bowtie2.sam > $file.bam; #samtools view (-b is to generate the output in BAM format, -h to include the header in the output, and -@ to allocate additional threads)
    echo "Showing the header of $file.bam";
    samtools view -H -@ 8 $file.bam; #samtools view (-H is to output only the header)
    echo "Providing simple statistics of $file.bam";
    samtools flagstat -@ 8 $file.bam > samstats-$file.err; #providing simple statistics on BAM files
    echo "Sorting $file.bam to $file.sorted.bam";
    samtools sort -@ 8 $file.bam > $file.sorted.bam; #sorting the mappings
    echo "Indexing $file.sorted.bam";
    samtools index -@ 8 $file.sorted.bam; #indexing sorted bam files
    echo "Converting $file.sorted.bam to ${file%%_*}_BPM.bw";
    bamCoverage -p 8 --binSize 5 --effectiveGenomeSize 119481543 --normalizeUsing BPM -b $file.sorted.bam -o ${file%%_*}_BPM.bw;
        #bamCoverage takes an alignment of reads or fragments as input (BAM file) and generates a coverage track (bigWig or bedGraph) as output
            ##   --numberOfProcessors, -p -> Number of processors to use
            ##   --binSize, -bs -> Size of the bins, in bases, for the output of the bigwig/bedgraph file
            ##   --effectiveGenomeSize -> The effective genome size is the portion of the genome that is mappable
            ##   --normalizeUsing -> Normalize the number of reads per bin using BPM (Bins Per Million mapped reads)
    done


    ##PAIR-END DATA ALIGNMENT
    echo "Mapping DAP_control to the reference genome with bowtie2"
    bowtie2 -p 8 --very-sensitive-local -t -x /home/gema/CNB/reference_genome/bowtie2_athaliana_genome -1 /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_fastp/DAP_control_1_fastp.fq -2 /home/gema/CNB/DAPseq/2022-03-29_OMalleyEtal_sequencing_data_fastp/DAP_control_2_fastp.fq -S DAP_control_bowtie2.sam > DAP_control-bowtie2-alg-stdout+stderr.txt 2>&1
        ##--very-sensitive-local -> Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    echo "Converting DAP_control_bowtie2.sam to DAP_control.bam"
    samtools view -b -h -@ 8 DAP_control_bowtie2.sam > DAP_control.bam #samtools view (-b is to generate the output in BAM format, -h to include the header in the output, and -@ to allocate additional threads)
    echo "Showing the header of DAP_control.bam"
    samtools view -H -@ 8 DAP_control.bam; #samtools view (-H is to output only the header)
    echo "Providing simple statistics of DAP_control.bam"
    samtools flagstat -@ 8 DAP_control.bam > samstats-DAP_control.err #providing simple statistics on BAM files
    echo "Sorting DAP_control.bam to DAP_control.sorted.bam"
    samtools sort -@ 8 DAP_control.bam > DAP_control.sorted.bam #sorting the mappings
    echo "Indexing DAP_control.sorted.bam"
    samtools index -@ 8 DAP_control.sorted.bam #indexing sorted bam files
    echo "Converting DAP_control.sorted.bam to DAP_control_BPM.bw"
    bamCoverage -p 8 --binSize 5 --effectiveGenomeSize 119481543 --normalizeUsing BPM -b DAP_control.sorted.bam -o DAP_control_BPM.bw
        #bamCoverage takes an alignment of reads or fragments as input (BAM file) and generates a coverage track (bigWig or bedGraph) as output
            ##   --numberOfProcessors, -p -> Number of processors to use
            ##   --binSize, -bs -> Size of the bins, in bases, for the output of the bigwig/bedgraph file
            ##   --effectiveGenomeSize -> The effective genome size is the portion of the genome that is mappable
            ##   --normalizeUsing -> Normalize the number of reads per bin using BPM (Bins Per Million mapped reads)