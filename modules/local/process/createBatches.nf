process createBatches {
  label 'seqkit'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(reads)
  val(batchSize)

  output:
  tuple val(meta),  path("*.fastq.gz"), emit: reads

  script:
  def bsizeOpts = batchSize ? "-L ${batchSize}" : ""
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  # merge all R1 fastq files per batchSize
  ls -1 "$reads"/*R1.fastq.gz | parallel -N $batchSize 'zcat {} | pigz -p ${task.cpus} -c  > '${prefix}'_batch{#}.R1.fastq.gz'
  # merge all R2 fastq files per batchSize
  ls -1 "$reads"/*R2.fastq.gz | parallel -N $batchSize 'zcat {} | pigz -p ${task.cpus} -c  > '${prefix}'_batch{#}.R2.fastq.gz'
  """
}