#!/usr/bin/bash
cd /home/gema/CNB/motif_discovery


#MOTIF DISCOVERY
peaksfileslist=($(ls /home/gema/CNB/motif_discovery | grep -E '10-500bp.fa' | cut -d_ -f1)) #storing the name of each file of readings (e.g., DAP_ATHB53) in an array
for file in ${peaksfileslist[@]}; do echo "Finding motifs of $file";
xstreme --p ${file}_top500_best_peaks_10-500bp.fa --oc 2022-05-13_motif-discovery_xstreme_$file --m JASPAR2022_CORE_non-redundant_pfms_meme.txt --dna --minw 6 --maxw 20 --meme-p 16 > $file-xstreme-alg-stdout+stderr.txt 2>&1;
    ##--dna -> Set the alphabet to DNA
    ##--minw 6 --maxw 20 -> Find motifs of 6-20 bp
    ##--meme-p 8 -> Number of processors to use
done

xstreme --p ChIP3_best_summits_regions.fa --oc 2022-05-13_motif-discovery_xstreme_ChIP3 --m JASPAR2022_CORE_non-redundant_pfms_meme.txt --dna --minw 6 --maxw 20 --meme-p 8 > ChIP3-xstreme-alg-stdout+stderr.txt 2>&1


#SEARCH FOR PEAKS THAT HAVE THE MOTIF OF INTEREST
    #TCP motif
    fimo --parse-genomic-coord --verbosity 1 --oc 2022-05-13_motif-discovery_xstreme_BRC1_fimo_out_TCP23 --bgfile 2022-05-13_motif-discovery_xstreme_BRC1/background 2022-05-13_motif-discovery_xstreme_BRC1/TCP23.meme BRC1_peaks.fa

    #G-box
    fimo --parse-genomic-coord --verbosity 1 --oc 2022-05-13_motif-discovery_xstreme_BRC1_fimo_out_1 --bgfile 2022-05-13_motif-discovery_xstreme_BRC1/background --motif 1-GACACGTGTC 2022-05-13_motif-discovery_xstreme_BRC1/streme_out/streme.xml BRC1_peaks.fa

    # HB motif
    fimo --parse-genomic-coord --verbosity 1 --oc 2022-05-13_motif-discovery_xstreme_HB21_fimo_out_1 --bgfile 2022-05-13_motif-discovery_xstreme_HB21/background --motif ACCAATWATTGDDNH 2022-05-13_motif-discovery_xstreme_HB21/meme_out/meme.xml HB21_peaks.fa
    fimo --parse-genomic-coord --verbosity 1 --oc 2022-05-13_motif-discovery_xstreme_HB40_fimo_out_1 --bgfile 2022-05-13_motif-discovery_xstreme_HB40/background --motif CAATTATTGGT 2022-05-13_motif-discovery_xstreme_HB40/meme_out/meme.xml HB40_peaks.fa
    fimo --parse-genomic-coord --verbosity 1 --oc 2022-05-13_motif-discovery_xstreme_HB53_fimo_out_1 --bgfile 2022-05-13_motif-discovery_xstreme_HB53/background --motif HCAATAATTGD 2022-05-13_motif-discovery_xstreme_HB53/meme_out/meme.xml HB53_peaks.fa