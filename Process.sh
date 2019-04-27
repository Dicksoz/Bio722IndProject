#!/bin/bash

nThreads=$1
Outdir=$2
Indir=$3

mkdir -p $Outdir;
cd $Outdir;

echo "Running Fastp";
mkdir fastpOut;
cd fastpOut;
for f in $Indir/*_R1_*.fastq.gz; do
   ../../runfastp.sh $f $nThreads ../../Adapters.fasta 
done
cd ..;

mkdir mergeOut;
parallel -j $nThreads 'sample=$(basename {} .log);
    cat {.}* > mergeOut/${sample}.fq.gz' ::: fastpOut/*.log;

echo "Deduplicating";
mkdir -p dedupOut;
parallel -j $nThreads 'sample=$(basename {} .fq.gz);
    zcat {} | perl /usr/local/prinseq/current/prinseq-lite.pl -fastq stdin -derep 14 -out_good stdout -log dedupOut/${sample}.log -out_bad null 2>/dev/null | gzip > dedupOut/${sample}.fq.gz;' ::: mergeOut/*.fq.gz;

echo "Mapping To Human Genome"
mkdir -p hmapOut;
for f in dedupOut/*.fq.gz; do
    sample=$(basename $f .fq.gz);
    bwa mem -t $nThreads ../BWAIndices/GRCh38 $f \
        2>hmapOut/${sample}.log | samtools sort -o hmapOut/${sample}.bam -;
done;

echo "Filtering Out Human Reads"
mkdir -p filterOut;
parallel -j $nThreads 'sample=$(basename {} .bam);
    samtools fastq -f 0x4 {} | gzip >filterOut/${sample}.fq.gz' ::: hmapOut/*.bam;

echo "Mapping Enriched Samples to Baits"
mkdir -p mapOut;
for f in filterOut/*ENR*.fq.gz; do
    sample=$(basename $f .fq.gz);
    bwa mem -t $nThreads ../BWAIndices/Regions $f \
        2>mapOut/${sample}.log | samtools sort -o mapOut/${sample}.bam -;
done;

echo "Mapping Shotgun Samples to Pathogen Genomes"
for f in filterOut/*SG*.fq.gz; do
    sample=$(basename $f .fq.gz);
    bwa mem -t $nThreads ../BWAIndices/WholeGenomes $f \
        2>mapOut/${sample}.log | samtools sort -o mapOut/${sample}.bam -;
done;

echo "Counting Enriched Reads"
mkdir countOut;
parallel -j $nThreads 'samtools view -h {} |
    perl ExactFilterSAM.pl -hSso 0.01 - |
    perl OrgCounting.pl ../HeaderMap.tab - > countOut/{/.}.count' ::: mapOut/*ENR*.bam;

echo "Counting Shotgun Reads"
parallel -j $nThreads 'samtools view -h {} |
    perl ExactFilterSAM.pl -hSs - |
    perl OrgCounting.pl ../HeaderMap.tab - > countOut/{/.}.count' ::: mapOut/*SG*.bam;

echo "Merging Count Files"
ls countOut | head -1 | xargs cut -f1 >Org.quant;
for f in countOut/*; do
    paste Org.quant <(cut -f2 $f) > tmp;
    mv -f tmp Org.quant;
done
parallel -j $nThreads 'echo {/.}' ::: OrgCounts/* | awk 'BEGIN {printf(Organism)}
{printf("\t"$0)} END {printf("\n")}' > tmp;
cat Org.quant >> tmp;
mv -f tmp Org.quant;

cd ..;

