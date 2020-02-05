#!/usr/bin/env bash
numb_of_files=$1
numbproc=$2
refgen=$3
aux=$4
out=$5
sample=$6
quality=$7

echo 'Aligning reads to the reference genome ...'
if [[ $numb_of_files == 1 ]]; then # IF SE READS
    bwa mem -v 1 -t $numbproc $refgen $aux/r1.2b.aln.fq |
	samtools sort --threads $numbproc -T ${aux}/${sample} > ${out}/${sample}.all.bam
    samtools view --threads $numbproc -b -q $quality ${out}/${sample}.all.bam > ${out}/${sample}.q${quality}.bam
    parallel "samtools index {}" ::: ${out}/${sample}.q${quality}.bam ${out}/${sample}.all.bam
fi
if [[ $numb_of_files == 2 ]]; then # IF PE READS
    bwa mem -v 1 -t $numbproc $refgen $aux/r1.2b.aln.fq $aux/r2.2b.aln.fq |
	samtools sort --threads $numbproc -T ${aux}/${sample} > ${out}/${sample}.all.bam
    samtools view --threads $numbproc -b -q $quality ${out}/${sample}.all.bam > ${out}/${sample}.q${quality}.bam
    parallel "samtools index {}" ::: ${out}/${sample}.q${quality}.bam ${out}/${sample}.all.bam
fi
echo 'Done'
