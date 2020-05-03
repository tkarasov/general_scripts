#!/bin/bash
#
#
#
#####################################################################
# Run map short read sequences to a reference
#
# Usage: bash bwa_mapping.sh <sample> <reference> <read1.fq> <read2.fq>
# or:
#       qsub -N sample -o plantRemoval.out -e plantRemoval.err run_plantRemoval.sh samplename read1.fq read.fq
# where:
#       sample = "Name of the sample to process, this string
#                     will be used as the read group in the bamfile"
#       reference = "Full path to reference genome"
#       read1.fq   = "full path to file with first pairs of sequencing"
#       read2.fq   = "full path to file with seconf pairs of sequencing"
#
#
####################################################################
#$ -N run_bwa
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



sample=$1

# map using bwa and output mapped reads

bwa mem -R "@RG\tID:$sample\tSM:$sample" -t 12 ${reference} ${read1} ${read2} 2> $sample.logs/bwa.err |\
    samtools view -b -F 4 | \
    samtools vew -q 20 | \
    samtools/samtools sort -m 970m - > "$sample".bam 2> $sample.logs/samsort.err
