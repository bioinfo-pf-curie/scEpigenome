// in each read name, add the barcode info at the end
process seqkitReplace {
  label 'seqkit'
  label 'lowCpu'
  label 'medMem'
  tag "$meta.id"

  input:
  tuple val(meta), path(dir)

  output:
  tuple val(meta), path("barcodedFastq/"), emit: reads
  path('versions.txt'), emit: versions

  script:
  """
  mkdir barcodedFastq/

  n1=\$(ls ${dir}/*R1*.fastq.gz | wc -l)
  n2=\$(ls ${dir}/*R2*.fastq.gz | wc -l)
  if [[ \$n1 != \$n2 ]] 
  then
    echo "Number of R1 file is different from R2 files. Please, be sure that the "R1/R2" tags are used in the file names"
    exit -1
  fi 

  for fastq in ${dir}/*R1*.fastq.gz
  do
  prefix=\$(basename \$fastq | sed -e 's/.fastq.gz//')
  base=\$(echo \$prefix | sed -e 's/.R[1,2].*\$//' | sed -e 's/_/-/g')
  seqkit replace -p " " -r '_'\$base' ' \$fastq > "barcodedFastq/"\$base"_R1.fastq"
  gzip "barcodedFastq/"\$base"_R1.fastq"
  done

  for fastq in ${dir}/*R2*.fastq.gz
  do
  prefix=\$(basename \$fastq | sed -e 's/.fastq.gz//')
  base=\$(echo \$prefix | sed -e 's/.R[1,2].*\$//' | sed -e 's/_/-/g')
  seqkit replace -p " " -r '_'\$base' ' \$fastq > "barcodedFastq/"\$base"_R2.fastq"
  gzip "barcodedFastq/"\$base"_R2.fastq"
  done

  echo \$(seqkit version 2>&1) > versions.txt
  """
}