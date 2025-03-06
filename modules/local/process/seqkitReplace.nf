// in each read name, add the barcode info at the end
process seqkitReplace {
  label 'seqkit'
  label 'medCpu'
  label 'medMem'
  tag "$meta.id"

  input:
  tuple val(meta), path(dir)
  path(sampleDescitpion)

  output:
  tuple val(meta), path("barcodedFastq/"), emit: reads
  path('versions.txt'), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  sampleDes = sampleDescitpion ? "${sampleDescitpion}" : "sampleDescitpion.txt"
  """
  mkdir barcodedFastq/

  n1=\$(ls ${dir}/*R1*.fastq.gz | wc -l)
  n2=\$(ls ${dir}/*R2*.fastq.gz | wc -l)
  if [[ \$n1 != \$n2 ]] 
  then
    echo "Number of R1 file is different from R2 files. Please, be sure that the "R1/R2" tags are used in the file names"
    exit -1
  fi 

  if [ -f ${sampleDes} ]; then
    for fastq in ${dir}/*R1*.fastq.gz
    do
    # Extract prefix
    prefix=\$(basename \$fastq | sed -e 's/.fastq.gz//')
    base=\$(echo \$prefix | sed -e 's/.R[1,2].*\$//')
    # Get prefix corresponding bioname in the 2nd column of the sample descritption
    # no _ is accepted in the bioname because it is used as field separator in read name !
    bioname=\$(grep \$base"|" ${sampleDes} | cut -f2 -d"|" | sed -e 's/_/--/g' )
    seqkit replace -p " " -r '_'\$bioname' ' \$fastq | pgzip -p ${task.cpus} -c > "barcodedFastq/"\$base"_R1.fastq.gz"
    done

    for fastq in ${dir}/*R2*.fastq.gz
    do
    # Extract prefix
    prefix=\$(basename \$fastq | sed -e 's/.fastq.gz//')
    base=\$(echo \$prefix | sed -e 's/.R[1,2].*\$//')
    # Get prefix corresponding bioname in the 2nd column of the sample descritption
    # no _ is accepted in the bioname because it is used as field separator in read name !
    bioname=\$(grep \$base"|" ${sampleDes} | cut -f2 -d"|" | sed -e 's/_/--/g' )
    seqkit replace -p " " -r '_'\$bioname' ' \$fastq | pgzip -p ${task.cpus} -c > "barcodedFastq/"\$base"_R2.fastq.gz"
    done

  else
    for fastq in ${dir}/*R1*.fastq.gz
    do
    prefix=\$(basename \$fastq | sed -e 's/.fastq.gz//')
    base=\$(echo \$prefix | sed -e 's/.R[1,2].*\$//' | sed -e 's/_/-/g')
    seqkit replace -p " " -r '_'\$base' ' \$fastq | pgzip -p ${task.cpus} -c > "barcodedFastq/"\$base"_R1.fastq.gz"
    done

    for fastq in ${dir}/*R2*.fastq.gz
    do
    prefix=\$(basename \$fastq | sed -e 's/.fastq.gz//')
    base=\$(echo \$prefix | sed -e 's/.R[1,2].*\$//' | sed -e 's/_/-/g')
    seqkit replace -p " " -r '_'\$base' ' \$fastq | pgzip -p ${task.cpus} -c > "barcodedFastq/"\$base"_R2.fastq.gz"
    done
  fi

  echo \$(seqkit version 2>&1) > versions.txt
  """
}