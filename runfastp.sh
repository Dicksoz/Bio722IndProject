#!/bin/bash

if [ $@ -eq 0 ]; then
  echo "Usage: ./runfastp.sh ForwardReads.fq.gz Num_Threads [Adapters.fasta]"
  echo "Note: Forward/reverse is assumed to be denoted with _R1_/_R2_";
  echo "Note: Fastp doesn't use more than 16 threads";
fi

fwd="$1";
rev="${fwd/_R1_/_R2_}";
sample=$(basename "$fwd" .fastq.gz);
sample="${sample%_R1*}"
nThreads=$2
adapters=$3

if [ -z $nThreads ]; then nThreads=2; fi;
if [ $nThreads -gt 16 ]; then nThreads=16; fi; 
if [ -z $adapters ]; then 
    adapters="--detect_adapter_for_pe";
else
    adapters="--adapter_fasta $adapters";
fi


fastp -i "$fwd" -I "$rev" -o "${sample}_1U.fq.gz" -O "${sample}_2U.fq.gz" \
    $adapters --n_base_limit 0 --length_required 30 \
    --merge --merged_out "${sample}_M.fq.gz" --correction --trim_tail1 1 \
    --html "${sample}.html" --json /dev/null -R $sample \
    --thread $nThreads 2> "${sample}.log"
