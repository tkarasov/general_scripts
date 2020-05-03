# Run picardtools to remove PCR and optical duplicates from alignment
# file
java -Xmx20g -jar \
    $picard \
    MarkDuplicates I="$sample".bam O="$sample"DEDUP.bam \
    M=$sample.logs/metrics.txt REMOVE_DUPLICATES=true 2> $sample.logs/picard.err
