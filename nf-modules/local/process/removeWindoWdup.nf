/* 
 *  RT deduplaication Workflow
 */

process removeWindowDup {
  tag "$meta.id"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(bam)
  
  output:
  tuple val(meta), path("*_rmDup.bam"), emit: bam
  path("*_rmDup.log"), emit: logs

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  ## Index BAM file
  samtools index ${bam}

  # window param
  if [ ! -z ${params.distDup} ]; then
	  rmDup.py -v -i ${bam} -o ${prefix}_rmDup.bam -d ${params.distDup} > ${prefix}_rmDup.log
  else
	  rmDup.py -v -i ${bam} -o ${prefix}_rmDup.bam > ${prefix}_rmDup.log
  fi
  """
}