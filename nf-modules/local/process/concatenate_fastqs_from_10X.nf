concatenate_fastqs_from_10X()
{ 
        fastq_dir=$1
        output_dir=$2
        sample_name=$3


        local log=$4/concatenate_fastqs_from_10X.log
	echo -e "Concatenating Fastqs From 10X together..."
        echo -e "Logs: $log"
        echo

        mkdir -p $output_dir

        echo "10X Fastq directory $fastq_dir" > $log
        echo "Concatenated Fastq directory $output_dir" >> $log
  
        # Remove any already existing FASTQ files
        cmd="rm -f ${output_dir}/${sample_name}.R1.fastq.gz"
        exec_cmd ${cmd} >> ${log} 2>&1
        cmd="rm -f ${output_dir}/${sample_name}.R2.fastq.gz"
        exec_cmd ${cmd} >> ${log} 2>&1
        cmd="rm -f ${output_dir}/${sample_name}.R3.fastq.gz"
        exec_cmd ${cmd} >> ${log} 2>&1
 
        for i in ${fastq_dir}/*_R1_*.fastq.gz
        do
                cmd="cat ${i} >> ${output_dir}/${sample_name}.R1.fastq.gz"
                exec_cmd ${cmd} >> ${log} 2>&1
        done

        for i in ${fastq_dir}/*_R2_*.fastq.gz
        do
                cmd="cat ${i} >> ${output_dir}/${sample_name}.R2.fastq.gz"
                exec_cmd ${cmd} >> ${log} 2>&1
        done
 
        for i in ${fastq_dir}*_R3_*.fastq.gz
        do
                cmd="cat ${i} >> ${output_dir}/${sample_name}.R3.fastq.gz"
                exec_cmd ${cmd} >> ${log} 2>&1
        done

	echo "Finished Concatenated Fastq" >> $log
	echo "" >> $log

}