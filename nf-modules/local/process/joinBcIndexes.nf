/*
 * Bowtie2 index barcode alignment
 */

process joinBcIndexes {
  tag "$meta.id"
  label 'bowtie2'
  label 'highCpu'
  label 'medMem'

  input:
  tuple val(meta), path(readsMatchingSorted) // *_ReadsMatchingSorted.txt
  tuple val(meta), path(count_index) // *_count_index.txt

  output:
  // correctly barcoded reads
  tuple val(meta), path("*_read_barcodes.txt"), emit: results
  // mqc of counts
  path("*_bowtie2.log"), emit: logs
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  #Join indexes 1 & 2 together (inner join)
  join -t\$' ' -1 1 -2 1 ${prefix}_indexB_ReadsMatchingSorted.txt ${prefix}_indexC_ReadsMatchingSorted.txt > tmp
  
  #Count matched index 1 & 2
  echo \$(wc -l tmp) | cut -d' ' -f1 > count_index_1_2
  
  #Join indexes (1 & 2) & 3 together to recompose full barcode (inner join)
  join -t\$' ' -1 1 -2 1 tmp ${prefix}_indexD_ReadsMatchingSorted.txt > final
  
  #Reformat & count matched index (1 & 2 & 3) <=> barcode
  awk '{print substr(\$1,1)\"\tBC\"substr(\$2,2)substr(\$3,2)substr(\$4,2);count++} ;END{print count > \"count_index_1_2_3\"}' final > ${prefix}_read_barcodes.txt
  
  ##Write logs
  n_index_1=\$(cat ${prefix}_indexB_count_index.txt)
  n_index_2=\$(cat ${prefix}_indexC_count_index.txt)
  n_index_3=\$(cat ${prefix}_indexD_count_index.txt)
  n_index_1_2=\$(cat count_index_1_2)
  n_index_1_2_3=\$(cat count_index_1_2_3)

  ## logs
  echo "## Number of matched indexes 1: \$n_index_1" > ${prefix}_bowtie2.log
  echo "## Number of matched indexes 2: \$n_index_2" >> ${prefix}_bowtie2.log
  echo "## Number of matched indexes 1 and 2: \$n_index_1_2" >> ${prefix}_bowtie2.log
  echo "## Number of matched indexes 3: \$n_index_3" >> ${prefix}_bowtie2.log
  echo "## Number of matched barcodes: \$n_index_1_2_3" >> ${prefix}_bowtie2.log

  rm count_index_1_2 count_index_1_2_3 tmp final
  """
}