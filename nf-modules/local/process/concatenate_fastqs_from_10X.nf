/*
 * Concatenate 4 samples per reads 
 */

process concatenate_fastqs_from_10X{ 
  tag "$meta.id"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(s1r1), path(s1r2), path(s1r3)
  tuple val(meta), path(s2r1), path(s2r2), path(s2r3)
  tuple val(meta), path(s3r1), path(s3r2), path(s3r3)
  tuple val(meta), path(s4r1), path(s4r2), path(s4r3)

  output:
  tuple val(meta), path ("*_R2.fastq.gz"), emit: barcodeRead
  tuple val(meta), path ("*_R1.fastq.gz"), path ("*_R3.fastq.gz"), emit: dnaRead
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  for sample in ${prefix}_S*_R1_*.fastq.gz
  do
        cat \$sample >> ${prefix}_R1.fastq.gz
  done

  for sample in ${prefix}_S*_R2_*.fastq.gz
  do
        cat \$sample >> ${prefix}_R2.fastq.gz
  done

  for sample in ${prefix}_S*_R3_*.fastq.gz
  do
        cat \$sample >> ${prefix}_R3.fastq.gz
  done
  """
}