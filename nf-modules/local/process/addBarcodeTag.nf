/*
 * Add barcode tag into bam files
 */

process addBarcodeTag {
  tag "$meta.id"
  label 'samtools'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bcList)

  output:
  tuple val(meta), path("*_flagged.bam"), emit: bam

  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  # Remove secondary aligned reads (256 <=> "not primary alignment") & If R1 is unmapped or multimapped (NH != 1), tag R1 & R2 with flag "4" <=> "unmapped" & "chr" = '*'
  samtools view -F 256 ${bam} | awk -v OFS='\t' 'NR%2==1{if(\$12==\"NH:i:1\"){mapped=1;print \$0} else{mapped=0;\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0}} NR%2==0{if(mapped==1){print \$0} else{\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0} }' > ${prefix}.sam

  # If read is mapped R1 & unmapped R2 -> set R2 position as '2147483647'
  cat ${prefix}.sam | awk -v OFS='\t' 'NR%2==1{print \$0} NR%2==0{if(\$3==\"*\"){\$4=2147483647;print \$0} else{print \$0} }' > ${prefix}_2.sam
  # Remove comments from the header that produce bugs in the count phase
  samtools view -H ${bam} | sed '/^@CO/ d' > ${prefix}_header.sam
  cat ${prefix}_2.sam >> ${prefix}_header.sam && mv ${prefix}_header.sam ${prefix}.sam && samtools view -b -@ ${task.cpus} ${prefix}.sam > ${prefix}_unique.bam
  rm -f ${prefix}_2.sam ${prefix}.sam
  
  # Keeping R1 aligned + R2 start as tag 'XS' (Switch from Paired End Bam to Single End Bam)
  samtools view ${prefix}_unique.bam | awk '{OFS = \"\t\" ; if(NR%2==1 && !(\$3==\"*\")) {R1=\$0} else if(NR%2==1){R1=0}; if(NR%2==0 && !(R1==0)){tagR2Seq=\"XD:Z:\"\$10; tagR2Pos=\"XS:i:\"\$4;print R1,tagR2Pos,tagR2Seq}}' > ${prefix}_unique.sam
  
  # Sort and join on read names reads barcoded and reads mapped to genome (barcode as tag 'XB') 
  sort -T ${params.tmpDir} --parallel=${task.cpus} -k1,1 ${prefix}_unique.sam > ${prefix}_unique_sorted.sam
  # filter out unbarcoded OR unmapped reads 
  join -1 1  -2 1  ${prefix}_unique_sorted.sam <(awk -v OFS=\"\t\" '{print \$1,\"XB:Z:\"\$2}' ${bcList}) > ${prefix}_flagged.sam
  sed -i 's/ /\t/g' ${prefix}_flagged.sam
  
  #Remove comments from the header that produce bugs in the count phase
  samtools view -H ${prefix}_unique.bam | sed '/^@CO/ d' > ${prefix}_header.sam
  cat ${prefix}_flagged.sam >> ${prefix}_header.sam && mv ${prefix}_header.sam ${prefix}_flagged.sam && samtools view -@ ${task.cpus} -b ${prefix}_flagged.sam > ${prefix}_flagged.bam
  
  #Cleaning
  rm -f ${prefix}_unique.bam ${prefix}_flagged.sam ${prefix}_unique_sorted.sam
  """
}