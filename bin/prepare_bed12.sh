#!/usr/bin/env bash

bamfile=$1
outpath=$2
output_name=$3
multi=$4

if [ "True" = "${multi}" ]; then

	# Keep only properly paired reads
	samtools view -f 0x2 ${bamfile} > ${outpath}/${output_name}_properpair.sam &&

	# Output sam header in file
	samtools view -H ${bamfile} > ${outpath}/${output_name}_properpair.sorted.sam &&

	# Sort file by name and HI tag to keep the right mates together
	echo 'sorting'
	sort -k1,1 -k13,13 ${outpath}/${output_name}_properpair.sam >> ${outpath}/${output_name}_properpair.sorted.sam &&

	rm ${outpath}/${output_name}_properpair.sam &&

	# convert sam to bam
	samtools view -b ${outpath}/${output_name}_properpair.sorted.sam > ${outpath}/${output_name}_properpair.sorted.bam &&

	rm ${outpath}/${output_name}_properpair.sorted.sam

elif [ "False" = "${multi}" ]; then
        echo 'Not sorting hit index';

	# Keep only properly paired reads
	samtools view -b -f 0x2 ${bamfile} > ${outpath}/${output_name}_properpair.bam &&
	
	# Sort file by name
	echo 'sorting'
	samtools sort -o ${outpath}/${output_name}_properpair.sorted.bam -n ${outpath}/${output_name}_properpair.bam &&

	rm ${outpath}/${output_name}_properpair.bam

else

	echo 'Please set option contains_multi=True or False'
	exit 1

fi &&


# convert paired bam to split bed12
echo 'converting bam to bed12'
pairedBamToBed12 -i ${outpath}/${output_name}_properpair.sorted.bam  > ${outpath}/${output_name}.bed12 &&

rm ${outpath}/${output_name}_properpair.sorted.bam &&

mkdir -p $outpath/${output_name}_chromo/ &&

cd $outpath/${output_name}_chromo/ &&

awk -F '\t' '{print > $1".bed12"}' ${outpath}/${output_name}.bed12 &&

exit

