/*
 * Add barcode tag into bam files
 */

process Macs2peakCalling {
  tag "$meta.id"
  label 'samtools'
  label 'highCpu'
  label 'highMem'

  input:
  bam

  output:
  broad_peaks_bed
  sharp_peaks_bed
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  # pseudo bulk
  macs2 callpeak {params.peak_calling_option} -t {bam} -n {output} 2> {log}
  # args = --nomodel --extsize 200 --keep-dup all --broad -f BAM -p 0.05 --max-gap 1000 
  cut -f1-3 {output}_peaks.*Peak | bedtools merge -d {params.peak_merging_option} -i /dev/stdin > {output} 2>> {log}
  # peak_merging_option = 5000
  bedtools sort 
  """
}