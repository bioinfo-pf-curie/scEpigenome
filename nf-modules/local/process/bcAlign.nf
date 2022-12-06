/*
 * Bowtie2 index barcode alignment
 */

process bcAlign {
  tag "$meta.id - ${index}"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads), val(index), path(bwt2Idx)

  output:
  tuple val(meta), path ("*${index}_ReadsMatchingSorted.txt"), emit: results
  tuple val(meta), path ("*${index}_count_index.txt"), emit: counts  
  path ("*Bowtie2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def size = params.barcodes[ index ].size
  def base = params.barcodes[ index ].base
  def prefix = task.ext.prefix ?: "${meta.id}"
  def oprefix = "${prefix}_${index}"
  if(params.darkCycleDesign == false) {
    start = params.barcodes[ index ].start_nodarkcycles
  } else {
    start = params.barcodes[ index ].start_darkcycles
  }
  """
  ##Extract three indexes from reads 
  # darkCycles design (==the first 4 bases are not read during the sequencing, the index begin at pos 1): 1 - 16 = index 1 ; 21 - 36 = index 2; 41 - 56 = index 3
  # not darkcycles design: 5 - 20 = index 1 ; 25 - 40 = index 2; 45 - 60 = index 3
  # => the start change but not the length
  gzip -cd  ${reads} | awk -v start_index_1=${start} -v size_index=${size} 'NR%4==1{print ">"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > ${oprefix}Reads.fasta

  #Map indexes (-f) against Index libraries (-x)
  bowtie2 \
    -x ${bwt2Idx}/${base} \
    -f ${oprefix}Reads.fasta \
    -N 1 -L 8 --rdg 0,7 --rfg 0,7 --mp 7,7 \
    --ignore-quals --score-min L,0,-1 -t \
    --no-unal --no-hd \
    -p ${task.cpus} > ${oprefix}Bowtie2.sam 2> ${oprefix}Bowtie2.log
  #Keep only reads that were matched by a unique index 1 + counting matched index1
  awk '/XS/{next} \$2!=4{print \$1,\$3}' ${oprefix}Bowtie2.sam > ${oprefix}ReadsMatching.txt 
  wc -l < ${oprefix}ReadsMatching.txt  > ${oprefix}_count_index.txt
  
  ##Sort indexes by read name: 
  sort -T ${params.tmpDir} --parallel=${task.cpus} -k1,1 ${oprefix}ReadsMatching.txt > ${oprefix}_ReadsMatchingSorted.txt 

  #delete useless files
  rm ${oprefix}ReadsMatching.txt ${oprefix}Bowtie2.sam ${oprefix}Reads.fasta
  
  ## version
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  """
}
