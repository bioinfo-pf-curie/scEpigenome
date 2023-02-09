/*
 * Macs2 - peak calling
 */

process macs2{
  tag "$meta.id"
  label 'macs2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai), path(controlBam), path(controlBai)
  val(effGenomeSize)
  path(peakCountHeader)

  output:
  path("*.xls"), emit: outputXls
  tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peaks
  path("*_macs2_peaks.size_mqc.tsv"), emit: mqc_generalStat_peaksize
  path("*_macs2_peaks.count_mqc.tsv"), emit: mqc // macs2 module 
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def ctrl = controlBam ? "-c ${controlBam}" : ''
  def outputSuffix = (args.contains('--broad')) ? "broadPeak" : "narrowPeak"
  """
  echo \$(macs2 --version 2>&1) &> versions.txt
  macs2 callpeak \\
    ${args} \\
    -t ${bam} \\
    -n ${prefix}_macs2 \\
    -g $effGenomeSize \\

  cat ${prefix}_macs2_peaks.${outputSuffix} | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peakCountHeader - > ${prefix}_macs2_peaks.count_mqc.tsv
  grep ^chr[0-9]  ${prefix}_macs2_peaks.xls | awk '{ total += \$4 } END { print "average peak size :" total/NR }' > ${prefix}_macs2_peaks.size_mqc.tsv
  """
}


