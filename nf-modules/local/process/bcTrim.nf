/*
 * Cutadapt trim R2 barcode part
 */

process bcTrim {
  tag "$meta.id"
  label 'cutadapt'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(r1), path(r2)

  output:
  tuple val(meta), file("*.R1.fastq.gz"), file("*_trimmed.R2.fastq.gz"), emit: reads
  tuple val(meta), file("*_trimmedR2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  # Trim linker + barcode from R2 reads for genome aligning	
  cutadapt -u ${params.barcode_linker_length} --cores=${task.cpus} ${r1} -o ${prefix}_trimmed.R2.fastq > ${prefix}_trimmedR2.log
  gzip ${prefix}_trimmed.R2.fastq 
  cutadapt --version &> versions.txt
  """
}
