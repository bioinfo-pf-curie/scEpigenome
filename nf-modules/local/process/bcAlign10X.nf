/*
 * Bowtie2 index barcode alignment
 */

process bcAlign10X {
  tag "$meta.id"
  label 'bowtie2'
  label 'highCpu'
  label 'medMem'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path ("*_read_barcodes.txt"), emit: bcNames
  path ("*Bowtie2.log"), emit: logs
  tuple val(meta), path ("*_bowtie2.log"), emit: counts
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  ##Extract 1 barcode from reads : 1 - 16 = index 1
  gzip -cd  ${reads} | awk -v start_index_1=${params.barcode10X_start} -v size_index=${params.barcode10X_len} 'NR%4==1{print \">\"substr(\$1,2)}; NR%4==2{print \$0}' > ${prefix}_indexes_1_Reads.fasta

  #Map INDEXES 1 against Index1 library
  bowtie2 -x ${params.barcodes10X_bwt2} \
          -f ${prefix}_indexes_1_Reads.fasta \
          -p ${task.cpus} \
          ${args}  > ${prefix}Bowtie2.sam 2> ${prefix}Bowtie2.log

  #Keep only reads that were matched by a unique index 1 + counting matched index1
  awk -v prefix=${prefix} '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > prefix\"_count_index.txt\"}' ${prefix}Bowtie2.sam > ${prefix}ReadsMatching.txt
  
  ##Sort indexes by read name: 
  sort -T ${params.tmpDir} --parallel=${task.cpus} -k1,1 ${prefix}ReadsMatching.txt > ${prefix}_ReadsMatchingSorted.txt

  awk '{print substr(\$1,1)\"\tBC\"substr(\$2,2)}' ${prefix}_ReadsMatchingSorted.txt > ${prefix}_read_barcodes.txt

  ##Write logs
  n_index_1=\$(cat ${prefix}_count_index.txt)
  echo "## Number of matched barcodes: \$n_index_1" >> ${prefix}_bowtie2.log

  #delete useless files
  rm ${prefix}ReadsMatching.txt ${prefix}Bowtie2.sam ${prefix}_indexes_1_Reads.fasta ${prefix}_ReadsMatchingSorted.txt
  ## version
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  """
}
