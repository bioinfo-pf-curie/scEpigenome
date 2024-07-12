process createBatches {
  label 'unix'
  label 'minCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(reads)
  val(batchSize)

  output:
  tuple val(meta),  path("*.fastq.gz"), emit: reads
  path('versions.txt'), emit: versions

  script:
  def bsizeOpts = batchSize ? "-L ${batchSize}" : ""
  """
  # merge all R1 fastq files per batchSize
  ls -1 $reads/*R1.fastq.gz | xargs $bsizeOpts echo | awk '{print "zcat " \$0 " > batch_"NR".R1.fastq"}' > create_batches.sh && bash create_batches.sh
  # merge all R2 fastq files per batchSize
  ls -1 $reads/*R2.fastq.gz | xargs $bsizeOpts echo | awk '{print "zcat " \$0 " > batch_"NR".R2.fastq"}' > create_batches.sh && bash create_batches.sh
  gzip *.fastq

  echo "gzip "\$(gzip --version | awk 'NR==1{print \$NF}') > versions.txt  
  """
}