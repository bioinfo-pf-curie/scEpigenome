/* 
 *  Weighted histogram to highlight empty droplets 
 */

process gtfToTSSBed {
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  
  input:
  path(gtf)

  output:
  path("*_TSS*.bed"), emit: bed

  script:
  """
  # Add +/-window arround genes start sites
  create_transcript_annotation.sh $gtf ${params.tssWindow}
  """
}