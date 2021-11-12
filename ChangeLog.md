## CoCo 1.0.0 - 2021/11/11
- ca: Added option `-f/--fraction` to merge overlapping embedded genes having the same biotype (probable duplicates). 
By default, coco ca will merge genes having >=85% overlap. The minimal overlap can be set to any values between to 
]0,1] by using the `-f/--fraction` option. This feature can be disabled with `-f/--fraction -1` 

## CoCo 0.2.8
- ca: Better handling of gtf files: faster and lighter reading, compatibility with RefSeq "gene" attribute,
handling of missing attributes

## CoCo 0.2.7
- cc: Added `--countReadPairs` argument for featureCounts v2.0.2 and higher when `-p/--paired` is specified

## CoCo 0.2.6
- Do not check dependencies before printing help
- ca: Do not correct very complex cases
- cb: Added `--bed12` option to keep intermediate bed12 files

## CoCo 0.2.5p1
- Adjusted samtools threads to be "Number of *additional* threads to use"

## CoCo 0.2.5
- Compatibility with pandas 1.0

## CoCo 0.2.4
- ca: Runs even if there are no embedded genes in the gtf file

## CoCo 0.2.3
- cc: Search repair path before running
- ca: Better exon numbering, handles gtf without transcript name

## CoCo 0.2.2  
- cc: More efficient reading of gtf file and better memory handling
- ca: Only show progress bar if specified by `-V/--verbose` option

## CoCo 0.2.1p4
- ca:Do not perform intron correction on the embedded genes overlapping exons from different genes

## CoCo 0.2.1p3
- Use samtools multi-threading
- Avoid negative counts
- Do not consider unmapped reads in multimapped part

## CoCo 0.2.1p2
- Skip multimapped correction when input doesn't have multimapped reads and using option `-c both`

## CoCo 0.2.1p1
- Corrected length of monoexonic genes contained in small noncoding RNA for TPM calculation

## CoCo 0.2.1
- cb: Better sorting
- Backward compatibility with python 3.4 (recommended 3.5.1 or higher)

## CoCo 0.2.0
Added intron correction
- ca : Creation of a second gtf file (.intron.gtf) only having the features overlapping embedded genes
- cc : Added distribute_embedded_counts