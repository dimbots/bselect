#!/bin/bash

# This script performs downstream analysis of a bowtie2 run and demonstrates the following steps. 
# 1) Extracts the properly paired reads (pp) from SAM files and saves it in pp_sample.sam file.
# 2) Calculates the mapping quality of pp reads and save the output in mapping_quality_results. (mapping quality 0-5, 5-10 and >10)
# 3) Calculates the alignment score of pp reads (sum and reads with alignment score = 0) and save the output in align_count.txt.
# 4) Calculates the number of ambiguous bases in the reference and the number of gap opens and gap extensions of properly paired reads.
# 5) Calculates the number of mismatches of all pp reads and save output in mismatches_count.txt.

# 1) Extract properly paired reads from SAM files. And calculte mapping quality.

for ((x = $1; x <= $2; x++))

do
  
    sample="sample_${x}.sam"
    pp="pp_sample_${x}.sam"
    samtools view -f 0x2 $sample > $pp

done

    wait

# Save results for mapping quality 0-5, 5-10 and >10 in mapq_results.tmp.  

for ((x = $1; x <= $2; x++))

do 
    
    sample="pp_sample_${x}.sam"
    reads=$(fgrep -v "@" $sample | awk '$3 != "*"' | wc -l)
    mapq_results=mapq_results_sample${x}.tmp

    quality_5=$(fgrep -v "@" $sample | awk '$3 != "*" && $5 <= 5' | wc -l)
        bc <<< "scale=4; $quality_5/$reads * 100" >> $mapq_results
    quality5_10=$(fgrep -v "@" $sample | awk '$3 != "*" && $5 > 5 && $5 <= 10' | wc -l)
        bc <<< "scale=4; $quality5_10/$reads * 100" >> $mapq_results
    quality_10=$(fgrep -v "@" $sample | awk '$3 != "*" && $5 > 10' | wc -l)
        bc <<< "scale=4; $quality_10/$reads * 100" >> $mapq_results

    wait
    
done

# Merge mapq_results files. Add header for each sample. 

for ((i = $1; i <= $2; i++))

do
    mapqfile="mapq_results_sample${i}.tmp"
    addheader="header_sample_${i}.tmp1"
    echo "SAMPLE_${i}" | cat - $mapqfile > $addheader
    cat *.tmp1 > mapping_quality_results.txt

    wait

done
    
# Remove intermadiate files
rm *.tmp *.tmp1

# 2) Calculates the alignment score of pp reads and save the output in align_count.txt

for ((i = $1; i <= $2; i++))

do

    samfile="pp_sample_${i}.sam"
    fgrep -v "@" $samfile | awk '{print $12}' | grep AS | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> align_count.tmp
    fgrep -v "@" $samfile | awk '{print $12}' | grep AS:i:0 | wc -l >> align_0.tmp
    paste align_count.tmp align_0.tmp >> alignment_count.txt
    rm align_count.tmp align_0.tmp

done

    wait

# Add headers to alignment_score_count.txt

    nl alignment_count.txt | awk '$1=$1 + 47164' > align_count.tmp
    echo -e "sample \t align_score \t align_0" | cat - align_count.tmp > align_count.txt
    rm alignment_count.txt align_count.tmp

# Calculate the number of ambiguous bases in the reference 
# and the number of gap opens and gap extensions of properly paired reads.

# Count the number of ambiguous bases of the reference

for ((i = $1; i <= $2; i++))

do

    samfile="pp_sample_${i}.sam"
    fgrep -v "@" $samfile | awk '{print $13}' | grep XN | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> N_ref1.tmp
    fgrep -v "@" $samfile | awk '{print $14}' | grep XN | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> N_ref2.tmp
    paste N_ref1.tmp N_ref2.tmp | awk '{print ($1 + $2)}' >> N_ref.txt
    rm N_ref1.tmp N_ref2.tmp
    
done

    wait

# Count the number of gap opens and gap extensions

for ((i = $1; i <= $2; i++))

do

    samfile="pp_sample_${i}.sam"
    fgrep -v "@" $samfile | awk '{print $15}' | grep XO | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> N_ref1.tmp
    fgrep -v "@" $samfile | awk '{print $16}' | grep XO | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> N_ref2.tmp    
    paste N_ref1.tmp N_ref2.tmp | awk '{print ($1 + $2)}' >> N_ref1.txt
    rm N_ref1.tmp N_ref2.tmp

done

    wait


for ((i = $1; i <= $2; i++))

do

    samfile="pp_sample_${i}.sam"
    fgrep -v "@" $samfile | awk '{print $16}' | grep XG | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> N_ref1.tmp
    fgrep -v "@" $samfile | awk '{print $17}' | grep XG | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> N_ref2.tmp
    paste N_ref1.tmp N_ref2.tmp | awk '{print ($1 + $2)}' >> N_ref2.txt
    rm N_ref1.tmp N_ref2.tmp

done

    wait


# Merge files into N_gaps.txt. Column 1 : ambiguous bases of the reference. 
# Column 2 : Gap opens. Column 3 : Gap extensions.

    paste N_ref.txt N_ref1.txt N_ref2.txt | nl | awk '$1=$1 + 47164' > N_gaps.tmp
    echo -e "sample \t N_bases \t g_open \t g_extension" | cat - N_gaps.tmp > N_gaps.txt

rm N_ref.txt N_ref1.txt N_ref2.txt N_gaps.tmp

# Calculate the number of mismatches of all properly paired reads
# Results are saved to mismatch_count.txt

for ((i = $1; i <= $2; i++))

do

    samfile="pp_sample_${i}.sam"
    fgrep -v "@" $samfile | awk '{print $14}' | grep XM | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> mis_count.tmp
    fgrep -v "@" $samfile | awk '{print $15}' | grep XM | cut -d ":" -f 3 | awk '{sum+=$1}END{print sum}' >> mis_count1.tmp
    paste mis_count.tmp mis_count1.tmp | awk '{print ($1 + $2)}' >> mis_count.txt
    rm mis_count.tmp mis_count1.tmp

done

    wait

# Add headers to mis_count.txt and create final file mismatch_count.txt
 
    nl mis_count.txt | awk '$1=$1 + 47164' > mis_count.tmp
    echo -e "sample \t mismatches" | cat - mis_count.tmp > mismatch_count.txt
    rm mis_count.tmp mis_count.txt

# end of script
