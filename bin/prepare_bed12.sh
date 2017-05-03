#!/usr/bin/env bash

bamfile=$1
outpath=$2
multi=$3

if [ "True" = "${multi}" ]; then

	# Keep only properly paired reads
	samtools view -f 0x2 $outpath/${bamfile}.bam > $outpath/${bamfile}_properpair.sam &&

	# Output sam header in file
	samtools view -H $outpath/${bamfile}.bam > $outpath/${bamfile}_properpair.sorted.sam &&

	# Sort file by name and HI tag to keep the right mates together
	sort -k1,1 -k13,13 $outpath/${bamfile}_properpair.sam >> $outpath/${bamfile}_properpair.sorted.sam &&

	rm $outpath/${bamfile}_properpair.sam &&

	#convert sam to bam
	samtools view -b $outpath/${bamfile}_properpair.sorted.sam > $outpath/${bamfile}_properpair.sorted.bam &&

	rm $outpath/${bamfile}_properpair.sorted.sam

elif [ "False" = "${multi}" ]; then
    echo 'Not sorting hit index';

	# Keep only properly paired reads
	samtools view -b -f 0x2 $outpath/${bamfile}.bam > $outpath/${bamfile}_properpair.bam &&
	
	# Sort file by name
	samtools sort -o $outpath/${bamfile}_properpair.sorted.bam -n $outpath/${bamfile}_properpair.bam &&

	rm $outpath/${bamfile}_properpair.bam 

else

	echo 'Please set option contains_multi=True or False'
	exit 1

fi &&


# convert paired bam to split bed12

pairedBamToBed12 -i $outpath/${bamfile}_properpair.sorted.bam  > $outpath/${bamfile}.bed12 &&

rm $outpath/${bamfile}_properpair.sorted.bam &&

cd $outpath &&

mkdir -p $outpath/chromo/ &&

awk -F '\t' '{print > "chromo/"$1".bed12"}' $outpath/${bamfile}.bed12 &&

exit

