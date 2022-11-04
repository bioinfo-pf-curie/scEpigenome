/*
 * Bowtie2 index barcode alignment
 */

process bcSubset {
  tag "$meta.id"
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
  """
  """
}