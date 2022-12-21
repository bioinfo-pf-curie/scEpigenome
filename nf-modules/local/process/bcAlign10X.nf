/*
 * Bowtie2 index barcode alignment
 */

process bcAlign10X {
  tag "$meta.id"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads)
  val (bwt2)

  output:
  tuple val(meta), path ("*_ReadsMatchingSorted.txt"), emit: results
  tuple val(meta), path ("*_count_index.txt"), emit: counts
  tuple val(meta), path ("*_read_barcodes.txt"), emit: bcNames
  path ("*Bowtie2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  ##Extract 1 barcode from reads : 1 - 16 = index 1
  gzip -cd  ${reads} | awk -v start_index_1=${params.barcode10X_start} -v size_index=${params.barcode10X_len} 'NR%4==1{print \">\"substr(\$1,2)}; NR%4==2{print \$0}' > ${prefix}_indexes_1_Reads.fasta

  #Map INDEXES 1 against Index1 library
  bowtie2 -x ${bwt2} \
          -f ${prefix}_indexes_1_Reads.fasta \
          -p ${task.cpus} \
          ${args}  > ${prefix}Bowtie2.sam 2> ${prefix}Bowtie2.log

  #Keep only reads that were matched by a unique index 1 + counting matched index1
  awk '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > \"/count_index.txt\"}' ${prefix}Bowtie2.sam > ${prefix}ReadsMatching.txt
  
  ##Sort indexes by read name: 
  sort -T ${params.tmpDir} --parallel=${task.cpus} -k1,1 ${prefix}ReadsMatching.txt > ${prefix}_ReadsMatchingSorted.txt

  awk '{print substr(\$1,1)\"\tBC\"substr(\$2,2)}' ${prefix}_ReadsMatchingSorted.txt > ${prefix}_read_barcodes.txt

  #delete useless files
  rm ${prefix}ReadsMatching.txt ${prefix}Bowtie2.sam ${prefix}Reads.fasta
  ## version
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  """
}
