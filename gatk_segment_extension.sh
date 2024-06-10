#!/bin/bash

# This script has been developed at the Institute of Biological Chemistry of the School of Exact and Natural Science (IQUIBICEN)
# of University of Buenos Aires.
#____________________________________________________________________________________________________________________________
# After Common CNVs filtering, using cnv-filter.sh, each CNV interval is extended by 50% length at each side

#Paths to Directories required for this pipeline:

output_dir=${HOME}/output_dir
newfiles= ${output_dir}/newfiles
unique_cnvs=${newfiles}/unique_cnvs
extended_segments=${unique_cnvs}/extended_segments
mergedfiles=${extended_segments}/mergedfiles
qsmergedfiles=${mergedfiles}/qsmergedfiles
cnvsqsfiltered=${qsmergedfiles}/cnvsqsfiltered

# Each CNV interval is extended by 50% length at each side:

if [ ! $(ls -A ${extended_segments}) ];then
for file in ${unique_cnvs}/*.bed;do
	cat $file|awk 'var=0.5*($3-$2){printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"8"\t"}{printf "%.0f\n", var}'| \
	awk 'min=$2-$9 {if (min <0) min=$2;else min=min;print $1"\t"min"\t"$3+$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\\"}'| \
	sed 's/chr//g'|sort  -nk1,1 -nk2,2|sed 's/^/chr/g' \
	> ${extended_segments}/$(basename $file .bed)_ext.bed;
done
else
	echo "Segment extension is already done"
fi

# Step to Merge overlapping intervals resulting from the previous step, using bedtools merge
#conda activate pipeline
for file in ${extended_segments}/*_ext.bed;do 
	bedtools merge -i $file \
	> ${mergedfiles}/$(basename $file _ext.bed)_merged.bed;
done

for file in ${mergedfiles}/*_merged.bed;do
        bedtools intersect -a $file -b  ${newfiles}/$(basename $file _merged.bed)_type.bed -wa -wb| \
        sort -nrk1,8|awk '!a[$2] ++'|sed 's/:/\t/g'|sed 's/chr//g'| \
        sort -nk1,1 -nk2,2|sed 's/^/chr/g'|cut -f1,2,3,7,8,10,11,12 > ${qsmergedfiles}/$(basename $file _merged.bed)_qsmerged.bed;
done

# In this step, a dynamic Quality Score (QS) is defined based on Copy Number and Number of Exons for each CNV.Files are saved in cnvfiltered directory. 
#conda activate pipeline

for file in ${qsmergedfiles}/*_qsmerged.bed;do
         awk '{if ($6==0 && $7<40) {print $0,"\t"400;} else if ($6==0 && $7>=40 && $7<100) {print $0,"\t"$7*10;} \
         else if ($6==0 && $7>=100) {print $0,"\t"1000;} else if ($6==1 && $7<10) {print $0,"\t"100;} \
         else if ($6==1 && $7>=10 && $7<100) {print $0,"\t"$7*10;} else if ($6==1 && $7>=100) {print $0,"\t"1000;} \
         else if ($6 >=3 && $7<13) {print $0,"\t"50;} else if ($6>=3 && $7>=13 && $7<100) {print $0,"\t"$7*4;} \
         else {print $0,"\t"100;}}' $file \
         > ${cnvsqsfiltered}/$(basename $file _qsmerged.bed)_filtered.bed;
done

# In this step, each CNV dynamic QS < QS is filtered out from the sample, just keeping high-quality CNVs for further analysis.

for file in ${cnvsqsfiltered}/*_filtered.bed;do 
        awk '{if ($5> $8) print $0}' $file \
        > ${cnvsqsfiltered}/$(basename $file _filtered.bed).bed;
done 

