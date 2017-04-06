BAM_FILE=$1
GTF=$2
OUTPUT_PREFIX=$3
STRAND=$4
THREAD=$5
PAIRED=$6

featureCounts \
--minOverlap 10 \
--largestOverlap \
-s $STRAND \
-C \
-T $THREAD \
$PAIRED \
-a $GTF \
-o $OUTPUT_PREFIX \
$BAM_FILE
