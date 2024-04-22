# Load modules
module load meme/5.4.1
module load bedtools
module load R/4.1.2
genome=/home/FCAM/ssun/Genome/hg38.fa
sizes=/home/FCAM/ssun/Genome/hg38.chrom.sizes

# First make this Background file:
fasta-get-markov -m 3 $genome > hg38_bkgrnd.txt # 3-ordered Markov Background Model on hg38 genome file

# GATA3 peak summit file ranked by peak Intensity (normalized via DESeq2)
wc -l GATA_ChIP_summits_final_DESeq2_ranked_intensity.bed
head GATA_ChIP_summits_final_DESeq2_ranked_intensity.bed

# Extend the peak region to 161bp window centered on the summit
slopBed -b 80 -i GATA_ChIP_summits_final_DESeq2_ranked_intensity.bed -g $sizes  | sort -k1,1 -k2,2n > GATA_ChIP_summits_final_DESeq2_ranked_intensity_161bpwindow.bed #96868

# convert to FASTA
fastaFromBed -fi $genome -bed GATA_ChIP_summits_final_DESeq2_ranked_intensity_161bpwindow.bed -fo GATA_ChIP_summits_final_DESeq2_ranked_intensity_161bpwindow.fasta #193736

# Sort peaks by intensity, and select top 5000 peaks. Convert to FASTA.
sort -nrk4,4 GATA_ChIP_summits_final_DESeq2_ranked_intensity_161bpwindow.bed | head -n 5000 > GATA_ChIP_top5000_161window.bed #5000
fastaFromBed -fi $genome -bed GATA_ChIP_top5000_161window.bed -fo GATA_ChIP_top5000_161window.fasta



# MEME-Round1
meme -p 72 GATA_ChIP_top5000_161window.fasta -oc GATA_ChIP_top5000_161window_meme_output -nmotifs 10 -csites 20000 -objfun classic -searchsize 0 -minw 4 -maxw 35 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000

wget https://raw.githubusercontent.com/sysunn/siyu_daily_update/main/December_2023/MEME_individual_from_db_python3.py
python MEME_individual_from_db_python3.py -i GATA_ChIP_top5000_161window_meme_output/meme.txt
#WGATAAVATCW_meme.txt -- top enriched GATA-like motif in round1 (E-value: 3.8e-740)
#GATAADVATCW_meme.txt                                           (E-value: 5.8e-397)
#HTTATCTNYHHATCT_meme.txt                                       (E-value: 1.0e-134)
#ATCWGATAAVN_meme.txt                                           (E-value: 1.0e-083)
mkdir individual_meme
cp WGATAAVATCW_meme.txt individual_meme/GATAmotif1_meme.txt
rm *meme.txt
# MAST
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif1_meme.txt GATA_ChIP_summits_final_DESeq2_ranked_intensity_161bpwindow.fasta > mast_GATA3_PSWM_in_peaks_round1.txt #27756
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round1.txt #This updated_parse_mast_to_coordinates.R is included in the GitHub folder
wc -l mast_GATA3_PSWM_in_peaks_round1.bed #27753
# Use `intersectBed` to find GATA3 peaks without motif 1. And convert to FASTA. 
intersectBed -v -a GATA_ChIP_summits_final_DESeq2_ranked_intensity_161bpwindow.bed -b mast_GATA3_PSWM_in_peaks_round1.bed > without_motifs_1.bed #69112
fastaFromBed -fi $genome -bed without_motifs_1.bed -fo without_motifs_1.fasta



# MEME-round2:
sort -nrk4,4 without_motifs_1.bed | head -n 5000 > without_motifs_1_top5000.bed # select top 5000
fastaFromBed -fi $genome -bed without_motifs_1_top5000.bed -fo without_motifs_1_top5000.fasta
meme -p 72 without_motifs_1_top5000.fasta -oc GATA_ChIP_without_motifs_1_top5000_161window_meme_output -nmotifs 10 -csites 20000 -objfun classic -searchsize 0 -minw 4 -maxw 35 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000
python MEME_individual_from_db_python3.py -i GATA_ChIP_without_motifs_1_top5000_161window_meme_output/meme.txt
#AGATBYTTATC_meme.txt -- top enriched GATA-like motif in round2 (E-value: 3.9e-629)
cp AGATBYTTATC_meme.txt individual_meme/GATAmotif2_meme.txt
rm *meme.txt
# MAST
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif2_meme.txt without_motifs_1.fasta > mast_GATA3_PSWM_in_peaks_round2.txt #16778
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round2.txt
wc -l mast_GATA3_PSWM_in_peaks_round2.bed #16775
# From peaks without motif1 exclude peaks with motif2
intersectBed -v -a without_motifs_1.bed -b mast_GATA3_PSWM_in_peaks_round2.bed > without_motifs_12.bed #52334
fastaFromBed -fi $genome -bed without_motifs_12.bed -fo without_motifs_12.fasta



# MEME-round3:
sort -nrk4,4 without_motifs_12.bed | head -n 5000 > without_motifs_12_top5000.bed # select top 5000
fastaFromBed -fi $genome -bed without_motifs_12_top5000.bed -fo without_motifs_12_top5000.fasta
meme -p 72 without_motifs_12_top5000.fasta -oc GATA_ChIP_without_motifs_12_top5000_161window_meme_output -nmotifs 10 -csites 20000 -objfun classic -searchsize 0 -minw 4 -maxw 35 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000
# MAST
python MEME_individual_from_db_python3.py -i GATA_ChIP_without_motifs_12_top5000_161window_meme_output/meme.txt
#AGATDDDNAGATAAD_meme.txt -- top enriched GATA-like motif in round3 (E-value: 3.6e-259)
cp AGATDDDNAGATAAD_meme.txt individual_meme/GATAmotif3_meme.txt
rm *meme.txt
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif3_meme.txt without_motifs_12.fasta > mast_GATA3_PSWM_in_peaks_round3.txt #9430
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round3.txt
wc -l mast_GATA3_PSWM_in_peaks_round3.bed #9427
intersectBed -v -a without_motifs_12.bed -b mast_GATA3_PSWM_in_peaks_round3.bed > without_motifs_123.bed #42904
fastaFromBed -fi $genome -bed without_motifs_123.bed -fo without_motifs_123.fasta


# MEME-round4:
sort -nrk4,4 without_motifs_123.bed | head -n 5000 > without_motifs_123_top5000.bed
fastaFromBed -fi $genome -bed without_motifs_123_top5000.bed -fo without_motifs_123_top5000.fasta
meme -p 72 without_motifs_123_top5000.fasta -oc GATA_ChIP_without_motifs_123_top5000_161window_meme_output -nmotifs 10 -csites 20000 -objfun classic -searchsize 0 -minw 4 -maxw 35 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000
# MAST
python MEME_individual_from_db_python3.py -i GATA_ChIP_without_motifs_123_top5000_161window_meme_output/meme.txt
#TTATCTBYNNVATCT_meme.txt -- top enriched GATA-like motif in round4 (E-value:5.3e-232)
cp TTATCTBYNNVATCT_meme.txt individual_meme/GATAmotif4_meme.txt
rm *meme.txt
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif4_meme.txt without_motifs_123.fasta > mast_GATA3_PSWM_in_peaks_round4.txt #7434
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round4.txt
wc -l mast_GATA3_PSWM_in_peaks_round4.bed #7431
intersectBed -v -a without_motifs_123.bed -b mast_GATA3_PSWM_in_peaks_round4.bed > without_motifs_1234.bed #35473
fastaFromBed -fi $genome -bed without_motifs_1234.bed -fo without_motifs_1234.fasta



# MEME-round5
sort -nrk4,4 without_motifs_1234.bed | head -n 5000 > without_motifs_1234_top5000.bed
fastaFromBed -fi $genome -bed without_motifs_1234_top5000.bed -fo without_motifs_1234_top5000.fasta
meme -p 72 without_motifs_1234_top5000.fasta -oc GATA_ChIP_without_motifs_1234_top5000_161window_20_meme_output -nmotifs 20 -csites 20000 -objfun classic -searchsize 0 -minw 4 -maxw 35 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000
# MAST
python MEME_individual_from_db_python3.py -i GATA_ChIP_without_motifs_1234_top5000_161window_20_meme_output/meme.txt
#NCTTATCWGAT_meme.txt -- top enriched GATA-like motif in round5 (E-value: 1.2e-065)
cp NCTTATCWGAT_meme.txt individual_meme/GATAmotif5_meme.txt
rm *meme.txt
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif5_meme.txt without_motifs_1234.fasta > mast_GATA3_PSWM_in_peaks_round5.txt #4654
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round5.txt
wc -l mast_GATA3_PSWM_in_peaks_round5.bed #4651
intersectBed -v -a without_motifs_1234.bed -b mast_GATA3_PSWM_in_peaks_round5.bed > without_motifs_12345.bed #30822
fastaFromBed -fi $genome -bed without_motifs_12345.bed -fo without_motifs_12345.fasta



# MEME-round6:
sort -nrk4,4 without_motifs_12345.bed | head -n 5000 > without_motifs_12345_top5000.bed
fastaFromBed -fi $genome -bed without_motifs_12345_top5000.bed -fo without_motifs_12345_top5000.fasta
meme -p 72 without_motifs_12345_top5000.fasta -oc GATA_ChIP_without_motifs_12345_top5000_161window_25_meme_output -nmotifs 25 -csites 20000 -objfun classic -searchsize 0 -minw 10 -maxw 35 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000
# MAST
python MEME_individual_from_db_python3.py -i GATA_ChIP_without_motifs_12345_top5000_161window_25_meme_output/meme.txt
#AGATBNNNNNNVNWGATAA_meme.txt -- top enriched GATA-like motif in round6 (E-value:2.7e-218)
cp AGATBNNNNNNVNWGATAA_meme.txt individual_meme/GATAmotif6_meme.txt
rm *meme.txt
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif6_meme.txt without_motifs_12345.fasta > mast_GATA3_PSWM_in_peaks_round6.txt #3810
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round6.txt
wc -l mast_GATA3_PSWM_in_peaks_round6.bed #3807
intersectBed -v -a without_motifs_12345.bed -b mast_GATA3_PSWM_in_peaks_round6.bed > without_motifs_123456.bed #27014
fastaFromBed -fi $genome -bed without_motifs_123456.bed -fo without_motifs_123456.fasta



# MEME-round7:
sort -nrk4,4 without_motifs_123456.bed | head -n 5000 > without_motifs_123456_top5000.bed
fastaFromBed -fi $genome -bed without_motifs_123456_top5000.bed -fo without_motifs_123456_top5000.fasta
meme -p 72 without_motifs_123456_top5000.fasta -oc try16_72_meme_output -nmotifs 25 -csites 20000 -objfun classic -searchsize 0 -minw 15 -maxw 30 -revcomp -dna -bfile hg38_bkgrnd.txt -maxsize 10000000 -wg 1 -ws 0.1
# MAST
python MEME_individual_from_db_python3.py -i try16_72_meme_output/meme.txt
#NAGATBNNARNNVWGATA_meme.txt -- top enriched GATA-like motif in round7 (E-value:7.2e-075)
cp NAGATBNNARNNVWGATA_meme.txt individual_meme/GATAmotif7_meme.txt
rm *meme.txt
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif7_meme.txt without_motifs_123456.fasta > mast_GATA3_PSWM_in_peaks_round7.txt #3022
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round7.txt
wc -l mast_GATA3_PSWM_in_peaks_round7.bed #3019
intersectBed -v -a without_motifs_123456.bed -b mast_GATA3_PSWM_in_peaks_round7.bed > without_motifs_1234567.bed #23994
fastaFromBed -fi $genome -bed without_motifs_1234567.bed -fo without_motifs_1234567.fasta



# MEME-round8:
sort -nrk4,4 without_motifs_1234567.bed | head -n 2000 > without_motifs_1234567_top2000.bed #select top 2000 #161win

# make a 81bp window top2000 peak file (in R)
module load R/4.1.2
R
library(bigWig)
options(scipen = 999)
GATA3_peak_summits=center.bed(read.table("without_motifs_1234567_top2000.bed", header=FALSE), upstreamWindow = 0, downstreamWindow = 0)
GATA3_peak_81win=center.bed(GATA3_peak_summits, upstreamWindow = 40, downstreamWindow = 40)
write.table(GATA3_peak_81win,file= "without_motifs_1234567_81win_top2000.bed", quote=F,sep="\t",col.names=F,row.names=F)
# Convert to FASTA
fastaFromBed -fi $genome -bed without_motifs_1234567_81win_top2000.bed -fo without_motifs_1234567_81win_top2000.fasta
meme -oc GATA3_without_mot_81bp_top2000_ws_0.5_wg_11_classic_zoops_input3_meme_output -nmotifs 20 -evt 0.05 -objfun classic -csites 20000 -searchsize 0 -minw 4 -maxw 20 -ws 0.5 -wg 11 -revcomp -dna -markov_order 3 -maxsize 100000000 without_motifs_1234567_81win_top2000.fasta
python MEME_individual_from_db_python3.py -i GATA3_without_mot_81bp_top2000_ws_0.5_wg_11_classic_zoops_input3_meme_output/meme.txt
#YTTATCTYYNNHVATCT_meme.txt -- top enriched GATA-like motif in round8 (E-value:7.0e-017)
cp YTTATCTYYNNHVATCT_meme.txt individual_meme/GATAmotif8_meme.txt
rm *meme.txt
# MAST-against the 161bp peaks
mast -mt 0.0005 -hit_list -best individual_meme/GATAmotif8_meme.txt without_motifs_1234567.fasta > mast_GATA3_PSWM_in_peaks_round8.txt #1639
Rscript /home/FCAM/ssun/scripts/updated_parse_mast_to_coordinates.R mast_GATA3_PSWM_in_peaks_round8.txt
wc -l mast_GATA3_PSWM_in_peaks_round8.bed #1636
intersectBed -v -a without_motifs_1234567.bed -b mast_GATA3_PSWM_in_peaks_round8.bed > without_motifs_12345678.bed #22358
fastaFromBed -fi $genome -bed without_motifs_12345678.bed -fo without_motifs_12345678.fasta



