#!/bin/bash

# This script has been developed at the Institute of Biological Chemistry of the School of Exact and Natural Science (IQUIBICEN)
# of University of Buenos Aires.
#____________________________________________________________________________________________________________________________
# After gGATK-CNV analysis with gatk_case.sh, this pipeline is used for sample and variant filtering

#1) Sample filtering: remove samples with more than 100 CNVs.
#2) Variant filtering: remove (common) CNVs with frequency > 1%. 

#Paths to Directories required for this pipeline:

output_dir=${HOME}/output_dir
case=${output_dir}/case
case-results=${case}/case-results
newfiles= ${output_dir}/newfiles
unique_cnvs=${newfiles}/unique_cnvs

# gGATK-CNV output are VCF segment files (zcat genotyped-segments-case-sample_0.vcf), one per sample, 
# saved in case-results directory. The file is edited as follows,generating _type.bed files in the newfiles folder.

if [ ! $(ls -A ${newfiles} ];then
	for file in ${case-results}/*;do
		cd $file;
		zcat genotyped-segments-case-sample_0.vcf| \
		grep -v "#"| \
		cut -f1,2,5,6,8,10| \
		sed 's/END=//g'|sed 's/:/\t/g'| \
		awk '{if ($3!=".") print $1"\t"$2"\t"$5"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"}'> \
		${newfiles}/(basename $file)_type.bed;
		cd ../;
done
else
	echo "Case VCF segment files were processed to access CNV calls"
fi
# 1) Sample filtering: Samples with more than 100 CNVs are removed from _type.bed files

for file in ${newfiles}/*_type.bed;do
	count=$(cat $file|wc -l);
	if [ $count -gt 100 ];
	then \
		rm $file;
	fi;
done

# 2) Variant filtering:
# this script concatenates the fields: chr/start/end/cnv_type in intervals with '@', then identical concatenated intervals 
# are sorted and counted. Counts greater than 3 are removed,keeping unfrequent CNVs,saved as unique_cnvs.bed

if [ ! -f  ${newfiles}/unique_cnvs.bed ];then
	for file in ${newfiles}/*_type.bed;do
       		cat $file| \
        	awk '{print $1,$2,$3,$4}'| \
        	sed 's/ /@/g';
	done|sort|uniq -c|sort -nrk1,1|awk -F " " '{if ($1 <= 2) print $2}'|sed 's/@/\t/g'|sed 's/^chr//g'| \
	sort -nk1,1 -nk2,2|sed 's/^/chr/g' > ${newfiles}/unique_cnvs.bed 
fi

# Search and filter samples for unfrequent CNVs, listed in unique_cnvs.bed:

if [ ! $(ls -A ${unique_cnvs}/*) ];then
	for file in ${newfiles}/*_type.bed;do
		 awk '(FNR==NR){a[$1,$2,$3,$4];next} (($1,$2,$3,$4) in a) {print $0}' ${newfiles}/unique_cnvs.bed $file \
		 > ${unique_cnvs}/$(basename $file _type.bed).bed;
	done 
else
	echo "Samples Common CNVs filtering already done "
fi

