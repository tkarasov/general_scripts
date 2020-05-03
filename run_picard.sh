#!/bin/bash
#
#
#
#####################################################################
# Run picard to remove duplicates
#
# Usage: bash run_picard.sh <sample>
# or:
#       qsub -N sample -o plantRemoval.out -e plantRemoval.err run_plantRemoval.sh samplename read1.fq read.fq
# where:
#       sample = "Name of the sample to process, this string
#                     will be used as the read group in the bamfile"
#
####################################################################

#$ -cwd
#
#  Import environment
#$ -V
$
# Request 12 cores
#$ -pe parallel 12
#
# Request 8Gb per process
#$ -l h_vmem=8G
#


# Run picardtools to remove PCR and optical duplicates from alignment
# file
java -Xmx20g -jar \
    $picard \
    MarkDuplicates I="$sample".bam O="$sample"DEDUP.bam \
    M=$sample.logs/metrics.txt REMOVE_DUPLICATES=true 2> $sample.logs/picard.err
