#!/bin/bash

splan=$1
minReads=$2
protocol=$3

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan | uniq)

## Table headers 
## Barcodes only for indrop 
if [[ $protocol =~  "indrop" ]]
then
    echo "Sample_id,Sample_name,Barcoded,Index 1 and 2 found not 3,Index 1 found not 2 and 3,Index 2 found not 1 and 3,Index 3 found not 1 and 2,No Index Found ~ genomic DNA" > scChIPseq_barcode.csv
fi

## Mapping
if [[ $protocol == "scchip" ]]
then
    echo "Sample_id,Sample_name,Deduplicated reads, Window duplicates,RT duplicates,PCR duplicates,Uniquely mapped not barcoded,Mapped to multiple loci,Unmapped" > scChIPseq_alignments.csv
else
    echo "Sample_id,Sample_name,Deduplicated reads, PCR duplicates,Uniquely mapped not barcoded,Mapped to multiple loci,Unmapped" > scChIPseq_alignments.csv
fi

## Summary table
# The column names have to be the same as the ID column in the multiqcConfig.yaml !!!!! 
echo -e "Sample_id,Sample_name,Tot_frag , Aligned, Aligned_Barcoded, Deduplicated_reads, Cells>minReads, Reads(median)/cell, FRiP" > scChIPseq_table.csv

for sample in $all_samples
do
    ## sample name
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan | uniq)

    # BOWTIE2
    if [[ $protocol =~  "indrop" ]]; then
        total_frag=`grep "reads; of these:" index/${sample}_indexBBowtie2.log | cut -f1 -d' ' `
    elif [[ $protocol == "sccut_10X" ]]; then
        total_frag=`grep "reads; of these:" index/${sample}Bowtie2.log | cut -f1 -d' ' `
    else # cellenone
        total_frag=`grep "reads; of these:" index/${sample}Bowtie2.log | cut -f1 -d' ' `
    fi

    match_index_1=$(grep -e "## Number of matched indexes 1:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_2=$(grep -e "## Number of matched indexes 2:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_1_2=$(grep -e "## Number of matched indexes 1 and 2:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_3=$(grep -e "## Number of matched indexes 3:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_barcode=$(grep -e "## Number of matched barcodes:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')

    # Remove RT & PCR duplicats
    uniquely_mapped_and_barcoded=$(grep -e "## Number of reads mapped and barcoded:" allDup/${sample}_allDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    pcr_duplicates=$(grep -e "## Number of pcr duplicates:" allDup/${sample}_allDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    # scChip ::: 
    rt_duplicates=$(grep -e "## Number of rt duplicates:" allDup/${sample}_allDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    R1_mapped_R2_unmapped=$(grep -e "## Number of R1 mapped but R2 unmapped:" allDup/${sample}_allDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    # sont-ils similaires ???????????:
    reads_after_pcr_rt_rm=$(grep -e "## Number of reads after PCR and RT removal (not R1 unmapped R2):" allDup/${sample}_allDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    window_dup=$(grep -e "## Number of duplicates:" removeWindowDup/${sample}_removeWindowDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')

    if [[ $protocol == "scchip_indrop" ]]
    then
        unique_reads=$(grep -e "## Number of reads after duplicates removal:" removeWindowDup/${sample}_removeWindowDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
        unique_reads_percent=$(echo "$unique_reads $total_frag" | awk ' { printf "%.2f", 100*$1/$2 } ')
    else
        # for scchip ended : duplicate number after PCR, RT and window
        unique_reads=$(echo "$uniquely_mapped_and_barcoded $pcr_duplicates" | awk ' { printf "%.2f", $1-$2 } ')
        unique_reads_percent=$(echo "$unique_reads $total_frag" | awk ' { printf "%.2f", 100*$1/$2 } ')
    fi

    ## Data for the barcode matching graph
    index_1_2_not_3=$(echo "$match_index_1_2 $match_barcode" | awk ' { printf "%.2f", $1-$2 } ')
    index_1_not_2_not_3=$(echo "$match_index_1 $index_1_2_not_3 $match_barcode" | awk ' { printf "%.2f", $1-$2-$3 } ')
    index_2_not_1_3=$(echo "$match_index_2 $match_index_1_2" | awk ' { printf "%.2f", $1-$2 } ')
    index_3_not_1_2=$(echo "$match_index_3 $match_barcode" | awk ' { printf "%.2f", $1-$2 } ')
    no_index_found=$(echo "$total_frag $match_barcode $index_1_2_not_3 $index_1_not_2_not_3 $index_2_not_1_3 $index_3_not_1_2" | awk ' { printf "%.2f", $1-$2-$3-$4-$5-$6 } ')
    uniquely_mapped_and_barcoded_percent=$(echo "$uniquely_mapped_and_barcoded $total_frag" | awk ' { printf "%.2f", 100*$1/$2 } ')

    # STAR 
    uniquely_mapped=`grep "Uniquely mapped reads number" star/${sample}Log.final.out | awk '{print $NF}'`
    uniquely_mapped_percent=`grep "Uniquely mapped reads %" star/${sample}Log.final.out | awk '{print $NF}' | sed -e 's/%//'`
    multimapped=$(grep -e "Number of reads mapped to multiple loci " star/${sample}Log.final.out | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
    multimapped_toomany=$(grep -e "Number of reads mapped to too many loci " star/${sample}Log.final.out | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')    

    ## Data for mapping - STAR
    total_mapped=$(echo "$uniquely_mapped $multimapped $multimapped_toomany" | awk ' { printf "%.2f", $1+$2+$3 } ')
    unmapped_count=$(echo "$total_frag $total_mapped" | awk ' { printf "%.2f", $1-$2 } ')
    total_unmapped_percent=$(echo "$unmapped_mismatches_percent $unmapped_tooshort_percent $unmapped_other_percent" | awk ' { printf "%.2f", $1+$2+$3 } ')
    uniquely_mapped_unbarcoded=$(echo "$uniquely_mapped $uniquely_mapped_and_barcoded" | awk ' { printf "%.2f", $1-$2 } ')
    multimapped=$(echo "$multimapped $multimapped_toomany" | awk ' { printf "%.2f", $1+$2 } ')
    unmapped=$unmapped_count

    ## Data for cell thresholds
    # total cells 
    nbCell=$(wc -l < cellThresholds/${sample}_rmDup.txt) #Barcodes found = 19133
    # nb cells with more than 1000 reads
    nbCellminReads=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.txt | awk -v limit=$minReads '$1>=limit && NR>1{c++} END{print c+0}')

	FRiP=$(grep "$sample" frip/${sample}_FRiP.tsv | awk '{print $2}')

    # Median reads per cell with more than 1000 reads
    if (( $nbCellminReads>1 ))
    then 
        awk -v limit=$minReads '$1>=limit && NR>1 {print $1}' cellThresholds/${sample}_rmDup.txt | sort -n > list_nbReads_overminReads
        mod=$(($nbCellminReads%2))
        if (( $mod == 0 ))
        then
            # get the first number that cut in two the number of read range
            line_first=$(( $nbCellminReads/2 ))
            first_num=$(sed "${line_first}q;d" list_nbReads_overminReads)
            # get the second number that cut in two the number of read range
            line_sec=$(( $line_first+1 ))
            sec_num=$(sed "${line_sec}q;d" list_nbReads_overminReads)
            
            #median=$( echo "scale=0; (($first_num+$sec_num)/2)" | bc -l )
            median=$(echo "$first_num $sec_num" | awk ' { printf "%.1f", ($1+$2)/2 } ')
        else
            #median=$( echo "scale=0; (($nbcell+1)/2)" | bc -l )
            line_median=$(( $nbCellminReads/2 ))
            median=$(sed "${line_median}q;d" list_nbReads_overminReads)
        fi
    else
        # if there is one cell take the first line == number of reads within this cell
        # if there is no cell, it will write 0
        awk -v limit=$minReads '$1>=limit && NR>1 {print $1}' cellThresholds/${sample}_rmDup.txt | sort -n > list_nbReads_overminReads
        median=$(cat list_nbReads_overminReads)
    fi

    # Barcode indexes summary table => only for indrop 
    if [[ $protocol =~ "indrop" ]]
    then
        echo "${sample},$sname,$match_barcode,$index_1_2_not_3,$index_1_not_2_not_3,$index_2_not_1_3,$index_3_not_1_2,$no_index_found" >> scChIPseq_barcode.csv
    fi

    # Duplicates summary table
    if [[ $protocol == "scchip" ]]
    then
        echo "${sample},$sname,$unique_reads,$window_dup,$rt_duplicates,$pcr_duplicates,$uniquely_mapped_unbarcoded,$multimapped,$unmapped" >> scChIPseq_alignments.csv
    else
        echo "${sample},$sname,$unique_reads,$pcr_duplicates,$uniquely_mapped_unbarcoded,$multimapped,$unmapped" >> scChIPseq_alignments.csv
    fi
    
    ## Summary table
    echo -e "${sample},$sname,$total_frag, $uniquely_mapped_percent, $uniquely_mapped_and_barcoded_percent, $unique_reads_percent, $nbCellminReads, $median, $FRiP" >> scChIPseq_table.csv

done

