BAM_FILE=$1
GTF=$2
OUTPUT_PREFIX=$3

featureCounts \
--minOverlap 10 \
--largestOverlap \
-s 1 \
-C \
-p \
-T 4 \
-a $GTF \
-o $OUTPUT_PREFIX \
$BAM_FILE
