#arg[1] = .bam name sorted file
#arg[2] = .gtf annotation file
#arg[3] = .fasta genome file
#arg[4] = output prefix

args <- commandArgs(trailingOnly = TRUE)
BAM=args[1]
GTF=args[2]
#GENOME=args[3]
OUTPUT=args[3]
OUTPUT.count=paste(OUTPUT,'.Rsubread.count',sep='')
OUTPUT.stat=paste(OUTPUT,'.Rsubread.stat',sep='')
#OUTPUT.junctions=paste(OUTPUT,'.Rsubread.junctions',sep='')

library(Rsubread)
featureCounts=featureCounts(BAM, annot.ext=GTF, isGTFAnnotationFile=TRUE,
              GTF.featureType="exon",GTF.attrType="gene_id", countChimericFragments=FALSE,
              largestOverlap=TRUE, isPairedEnd=TRUE, useMetaFeatures=TRUE,
              requireBothEndsMapped=TRUE,strandSpecific=1, autosort=FALSE, allowMultiOverlap=FALSE,
              nthreads=7, reportReads=TRUE, juncCounts=FALSE,fraction=FALSE,countMultiMappingReads=FALSE)

head(featureCounts$counts)

#Junctions=featureCounts$counts_junction[!(is.na((featureCounts$counts_junction)$PrimaryGene)),]


write.table(as.data.frame(featureCounts$counts),
          file=OUTPUT.count, quote=FALSE, col.names=FALSE)

write.table(as.data.frame(featureCounts$stat),
            file=OUTPUT.stat, quote=FALSE)

#write.table(as.data.frame(Junctions),
#            file=OUTPUT.junctions, quote=FALSE, col.names=FALSE)
