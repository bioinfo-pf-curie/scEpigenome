/* 
 *  removeBlackRegions Workflow
 */

process removeBlackRegions {
  tag "$meta.id"
  label 'bedtools'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(bam)
  path (blackListBed)
  
  output:
  path ("versions.txt"), emit: versions
  tuple val(meta), path("*BlackReg.bam"), emit: bam
  tuple val(meta), path("*BlackReg.bam.bai"), emit: bai
  tuple val(meta), path("*_rmDup.sam"), emit: sam

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  # Removing encode black regions 
  if [[ "${params.removeBlackRegions}" == "true" ]]
  then
    samtools index ${bam}
    bedtools intersect -v -abam ${bam} -b ${blackListBed} > ${prefix}_rmDup_rmBlackReg.bam
    samtools index ${prefix}_rmDup_rmBlackReg.bam
  else
    cp ${bam} ${prefix}_rmDup_withBlackReg.bam
    samtools index ${prefix}_rmDup_withBlackReg.bam
  fi

  samtools view ${prefix}*BlackReg.bam > ${prefix}_rmDup.sam

  bedtools --version &> versions.txt
  """
}