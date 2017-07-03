TAD_file=$1
GENE_file=$2
CNV_file=$3
TAD_CNV_GENE_intersect=$4
# TAD_file format: chrom, start, end, tad_id
# GENE_file format: chrom, start, end, gene_id, gene_strand
# CNV_file format: chrom, start, end, segment_mean, donor_id, specimen_id, sample_id

bedtools window -a $TAD_file -b $GENE_file | sort -k4,4 > temp.TAD_GENE.intersect.txt
bedtools window -a $TAD_file -b $CNV_file | sort -k4,4 > temp.TAD_CNV.intersect.txt

join -t $'\t' temp.TAD_CNV.intersect.txt temp.TAD_GENE.intersect.txt -14 -24 > temp.TAD_CNV_GENE.intersect.txt

cat temp.TAD_CNV_GENE.intersect.txt | awk -F $'\t' '{if($6<$3 && $7<$4){print $0"\t"1} \
else if($3<$6 && $7<$4){print $0"\t"2} \
else if($3<$6 && $4<$7){print $0"\t"3} \
else{print $0"\t"4}}' > temp.TAD_CNV_GENE.TAD_genotype.txt

cat temp.TAD_CNV_GENE.TAD_genotype.txt | awk -F $'\t' '{if($20==1){if($7<$16){print $0"\t"11"\t"1}else{print $0"\t"12"\t"2}} \
else if($20==2){if($7<$16){print $0"\t"21"\t"1}else if($17<$6){print $0"\t"23"\t"1}else{print $0"\t"22"\t"3}} \
else if($20==3){if($17<$6){print $0"\t"31"\t"1}else{print $0"\t"32"\t"2}} \
else{print $0"\t"4"\t"4}}' |sort -k22,22n | awk -F $'\t' 'BEGIN{OFS="\t"}{print $1,$8,$9,$10,$11,$18,$20,$21}'>  $TAD_CNV_GENE_intersect
rm temp.*
