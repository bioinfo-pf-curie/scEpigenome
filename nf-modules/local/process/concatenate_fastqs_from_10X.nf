/*
 * Concatenate 4 samples per reads 
 */

process concatenate_fastqs_from_10X{ 
  tag "$meta.id"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  path(reads)

  output:
  tuple val(meta), path ("*_R2.fastq.gz"), emit: barcodeRead
  tuple val(meta), path ("*_R1.fastq.gz"), path ("*_R3.fastq.gz"), emit: dnaRead
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  for sample in ${prefix}_S*_R1_*.fastq.gz
  do
        cat ${sample} >> ${prefix}_R1.fastq.gz
  done

  for sample in ${prefix}_S*_R2_*.fastq.gz
  do
        cat ${sample} >> ${prefix}_R2.fastq.gz
  done

  for sample in ${prefix}_S*_R3_*.fastq.gz
  do
        cat ${sample} >> ${prefix}_R3.fastq.gz
  done
  """
}