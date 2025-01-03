#!/bin/bash

function usage {
    echo -e "usage : stats2multiqc.sh -s SAMPLE_PLAN -p PROTOCOL [-th]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo
    echo "stat2multiqc.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -s SAMPLE_PLAN"
    echo "   -p PROTOCOL"
    echo "   [-t]: min reads per cell"
    echo "   [-h]: help"
    exit;
}

while getopts "s:p:t:h" OPT
do
    case $OPT in
        s) splan=$OPTARG;;
        d) protocol=$OPTARG;;
        t) minReads=$OPTARG;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if  [[ -z $splan ]]; then
    usage
    exit
fi

all_samples=$(awk -F, '{print $1}' $splan | uniq )
n_header=0

for sample in $all_samples
do
                                                                                                                                                                                                          
    ## sample name
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan | uniq)
    header="Sample_id,Sample_name"
    output="${sample},${sname}"

    nb_frag=0
    nb_frag_barcoded=0
    if [[ -d barcodes ]]; then
        echo "for scChIP, scCutIndrop & scCut10x"
        echo "in if barcodes/"
        for batches in $(ls barcodes/${sample}*_addbarcodes.log)
        do
            echo $batches
            nb_frag_part=$(awk  '$0~"Total"{print $NF}' $batches)
            nb_barcoded_part=$(awk  '$0~"with barcodes"{print $NF}' $batches)
            nb_frag=$(( $nb_frag + $nb_frag_part ))
            nb_frag_barcoded=$(( $nb_frag_barcoded + $nb_barcoded_part ))
        done
        nb_reads=$(( $nb_frag * 2 ))
        nb_reads_barcoded=$(( $nb_frag_barcoded * 2 ))
        perc_barcoded=$(echo "${nb_reads_barcoded} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_frag,Number_of_reads,Number_barcoded_reads,Percent_barcoded"
        output+=",${nb_frag},${nb_reads},${nb_reads_barcoded},${perc_barcoded}"
    else
        echo "only for plate protocol"
        nb_reads=$(grep "raw total sequences" stats/${sample}.stats | awk '{print $5}')
        nb_frag=$(( $nb_reads / 2 ))
        nb_reads_barcoded=$nb_reads
        perc_barcoded=$(echo "${nb_reads_barcoded} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_frag,Number_of_reads,Number_barcoded_reads,Percent_barcoded"
        output+=",${nb_frag},${nb_reads},${nb_reads_barcoded},${perc_barcoded}"
    fi

    ## Mapped
    nb_reads_mapped=$(grep "reads mapped:" stats/${sample}.stats | awk '{print $4}')
    nb_paired_mapped=$(grep "reads mapped and paired:" stats/${sample}.stats | awk '{print $6}')
    nb_single_mapped=$(( $nb_reads_mapped - $nb_paired_mapped ))

    nb_paired_filter=$(grep "with itself and mate mapped" stats/${sample}_filtered.flagstats | awk '{print $1}')
    nb_single_filter=$(grep "singletons" stats/${sample}_filtered.flagstats | awk '{print $1}')
    nb_reads_filter=$(( $nb_paired_filter + $nb_single_filter ))

    perc_mapped=$(echo "${nb_reads_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_filter=$(echo "${nb_reads_filter} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    header+=",Number_of_aligned_reads,Percent_of_aligned_reads,Number_reads_after_filt,Percent_reads_after_filt"
    output+=",${nb_reads_mapped},${perc_mapped},${nb_reads_filter},${perc_filter}"

    ## Duplicates
    nb_reads_dups=$(grep "primary duplicates" stats/${sample}_markdup.flagstats | awk '{print $1}')
    perc_dups=$(echo "${nb_reads_dups} ${nb_reads_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    header+=",Number_of_duplicates_reads,Percent_of_duplicates"
    output+=",${nb_reads_dups},${perc_dups}"

    ## Filtering stats
    nb_not_barcoded=$(echo "${nb_reads} ${nb_reads_barcoded}" | awk ' { printf "%.*f",0,$1-$2 } ')
    nb_not_aligned=$(echo "${nb_reads_barcoded} ${nb_reads_mapped}" | awk ' { printf "%.*f",0,$1-$2 } ')
    nb_mapq_filters=$(echo "${nb_reads_mapped} ${nb_reads_filter} ${nb_reads_dups}" | awk ' { printf "%.*f",0,$1-($2+$3) } ')
    echo -e "Not barcoded\t$nb_not_barcoded" > ${sample}_filteringstats.mqc
    echo -e "Not aligned\t$nb_not_aligned" >> ${sample}_filteringstats.mqc
    echo -e "Duplicates\t$nb_reads_dups" >> ${sample}_filteringstats.mqc
    echo -e "Mapping quality\t$nb_mapq_filters" >> ${sample}_filteringstats.mqc
    echo -e "Final reads\t$nb_reads_filter" >> ${sample}_filteringstats.mqc

    ## Data for the barcode matching graph
    ## FRAG ::::::::
    #match_index_1=$(grep -e "## Number of matched indexes 1:" barcode/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    #match_index_2=$(grep -e "## Number of matched indexes 2:" barcde/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    #match_index_1_2=$(grep -e "## Number of matched indexes 1 and 2:" barcode/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    #match_index_3=$(grep -e "## Number of matched indexes 3:" barcode/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    #match_barcode=$(grep -e "## Number of matched barcodes:" barcode/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    #index_1_2_not_3=$(echo "$match_index_1_2 $match_barcode" | awk ' { printf "%.2f", $1-$2 } ')
    #index_1_not_2_not_3=$(echo "$match_index_1 $index_1_2_not_3 $match_barcode" | awk ' { printf "%.2f", $1-$2-$3 } ')
    #index_2_not_1_3=$(echo "$match_index_2 $match_index_1_2" | awk ' { printf "%.2f", $1-$2 } ')
    #index_3_not_1_2=$(echo "$match_index_3 $match_barcode" | awk ' { printf "%.2f", $1-$2 } ')
    #no_index_found=$(echo "$total_frag $match_barcode $index_1_2_not_3 $index_1_not_2_not_3 $index_2_not_1_3 $index_3_not_1_2" | awk ' { printf "%.2f", $1-$2-$3-$4-$5-$6 } ')

    ## Data for cell thresholds
    # total cells 
    #nbCell=$(wc -l < cellThresholds/${sample}_rmDup.txt) #Barcodes found = 19133
    # nb cells with more than 1000 reads
    #nbCellminReads=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.txt | awk -v limit=$minReads '$1>=limit && NR>1{c++} END{print c+0}')

    if [[ -e frip/${sample}_FRiP.tsv ]]
    then
	FRiP=$(grep "$sample" frip/${sample}_FRiP.tsv | awk '{print $2}')
	header+=",Fraction_of_reads_in_peaks"
	output+=",${frip}"
    fi
    #peakSizes=$(cut -f2 -d: peakSizes/${sample}_macs2_peaks.size_mqc.tsv)

    # Median reads per cell with more than 1000 reads
    countsfiles=$(ls finalBarcodeCounts/${sample}_final_barcodes_counts.txt)
    if [[ -e "${countsfiles[0]}" ]]
    then
	nbCell=$(wc -l ${countsfiles[0]} | awk '{print $1}')
	nbCellminReads=$( awk -v limit=$minReads '$1>=limit{c++} END{print c}' ${countsfiles[0]})
	header+=",Cell_number,Cell_number_minReads"
	output+=",${nbCell},${nbCellminReads}"
	if (( $nbCellminReads>1 ))
	then
        #int(i/2) = indice du milieu de tableau
	    median=$(sort -k1,1n ${countsfiles[0]} | awk -v limit=$minReads '$1>=limit{a[i++]=$1} END { print a[int(i/2)] }')
	    header+=",Median_reads_per_cell"
	    output+=",${median}"
	fi
    fi

    if [ $n_header == 0 ]; then
        echo -e $header > general_stats.mqc
        n_header=1
    fi
    echo -e $output >> general_stats.mqc
done

