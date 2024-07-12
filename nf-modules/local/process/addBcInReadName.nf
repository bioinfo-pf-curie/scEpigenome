process addBcInReadName {
  label 'seqkit'
  label 'lowCpu'
  label 'medMem'
  tag "$meta.id"

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("barcodedFastq/"), emit: reads

  script:
  """
  mkdir barcodedFastq/
  for fastq in ${reads}/*.fastq.gz
  do
  base=\$(echo \$fastq | grep -o .R[1,2].fastq.gz)
  prefix=\$(basename \$fastq \$base)
  seqkit replace -p " " -r '_'\$prefix' '  \$fastq > "barcodedFastq/"\$prefix\$base
  gzip "barcodedFastq/"\$prefix\$base
  done

  echo \$(seqkit version 2>&1) > versions.txt
  """
}