CNV_file=$1;
processed_CNV_file=$2;
threshold=$3;
# $12:chrom ,$13:start ,$14:end ,$10:segment_mean ,$1:icgc_donor_id ,$3:icgc_specimen_id ,$4:icgc_sample_id

cat $CNV_file | awk -F $'\t' '{if(NR>1){print $12"\t"$13"\t"$14"\t"$10"\t"$1"\t"$3"\t"$4}}'| \
sort -k7,7 -k1,1 -k2,2n -k3,3n | awk -F $'\t' '{if(strtonum($4)>0.2 || strtonum($4)<-0.2){print $0}}' > filtered.project_code.tsv
# $1:chrom ,$2:start ,$3:end ,$4:segment_mean ,$5:icgc_donor_id ,$6:icgc_specimen_id ,$7:icgc_sample_id
for SA_id in `cat filtered.project_code.tsv | awk -F $'\t' '{print $7}'| uniq`
do
cat filtered.project_code.tsv | grep -P "$SA_id$" > $SA_id.tmp.tsv
bedtools merge -d 10000 -i $SA_id.tmp.tsv -c 4,5,6,7 -o mean,distinct,distinct,distinct > $SA_id.merged.tmp.tsv
rm $SA_id.tmp.tsv
done

rm filtered.project_code.tsv

cat *.merged.tmp.tsv > $processed_CNV_file
rm *.merged.tmp.tsv