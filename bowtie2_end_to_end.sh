#!/bin/bash

# Quality control check of raw sequence data with FastQC.
# Mapping with Bowtie2 (end to end mode). Samples must have the following format: "sample_47165_r1.fastq".
# Proccessing the output for use with samtools. Statistics of all samples are saved to flagstat_results.txt.
# Sample range number have to be specified as numeric arguements.
# count_script.sh initiates as soon as bowtie2_end_to_end.sh finish.

# FastQC run

#for i in *.fastq

#do

 #   fastqc $i

#done

# Run bowtie2 (end to end mode)

for ((x = $1; x <= $2; x++))

do 

    r1="sample_${x}_r1.fastq"
    r2="sample_${x}_r2.fastq"
    bowtie2 -p 6 --end-to-end --score-min L,15,-0.6 --mp 7,7 --ignore-quals --rdg 0,7 -I 100 -X 1000 -x /project/bselect/botskar/illumina_POOL_28_CADK3ANXX_8_Truseq_2016_M4138_betama/RefBeet_bowtie2_index/RefBeet_bowtie2 -S sample_${x}.sam -1 $r1 -2 $r2 >> log_file 2>&1

done

    wait
    
# Convert the SAM file into a BAM files with samtools

for ((i = $1; i <= $2; i++))

do

    samfile="sample_${i}.sam"
    bamfile="sample_${i}.bam"
    nohup samtools view -bS $samfile > $bamfile & 

done

    wait    

# Sort the BAM files

for ((y = $1; y <= $2; y++))

do

    bamfile="sample_${y}.bam"
    bamsorted="sample_${y}_sorted.bam"
    nohup samtools sort -o $bamsorted  $bamfile &

done

    wait
    
# Index the BAM sorted files

for ((z = $1; z <= $2; z++))

do

    sortedbam="sample_${z}_sorted.bam"
    nohup samtools index $sortedbam &

done  

    wait

# Use samtools flagstat for simple stats.

for ((e = $1; e <= $2; e++))

do

    sortedbam="sample_${e}_sorted.bam"
    nohup samtools flagstat $sortedbam > flagstat_output_sample_${e}.tmp &    
     
done

wait

# Merging flagstat files. Output is "flagstat_results.txt"

for ((i = $1; i <= $2; i++))

do
    flagstatfile="flagstat_output_sample_${i}.tmp"
    addheader="header_sample_${i}.tmp1"
    echo "SAMPLE_${i}" | cat - $flagstatfile > $addheader
    cat *.tmp1 > flagstat_results.txt
  
done

    wait

# Remove intermediate files

rm *.tmp *.tmp1

# Initiate count_script.sh scirpt

./count_script.sh $1 $2

