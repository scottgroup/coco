#!/usr/bin/env bash

bamfile=$1
outpath=$2
output_name=$3
threads=$4

read repair < $(dirname $0)/repair_path.txt;

# Keep only properly paired reads
samtools view -bf 0x2 -@ ${threads} ${bamfile} > ${outpath}/${output_name}_properpair.bam &&

# Sort file by name and HI tag to keep the right mates together
echo 'sorting' &&

$repair -T $threads -i ${outpath}/${output_name}_properpair.bam -o ${outpath}/${output_name}_properpair.sorted.bam &&

rm ${outpath}/${output_name}_properpair.bam &&

# convert paired bam to split bed12
echo 'converting bam to bed12' &&

pairedBamToBed12 -i ${outpath}/${output_name}_properpair.sorted.bam  > ${outpath}/${output_name}.bed12 &&

rm ${outpath}/${output_name}_properpair.sorted.bam &&

mkdir -p $outpath/${output_name}_chromo/ &&

cd $outpath/${output_name}_chromo/ &&

awk -F '\t' '{print > $1".bed12"}' ${outpath}/${output_name}.bed12 &&

exit

