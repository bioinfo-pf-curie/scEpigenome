/*
 * FRiP - Fraction of reads in Peaks
 */

process frip{
  label 'bedtools'
  label 'medCpu'
  label 'medMem'
  tag("${meta.id}")

  input:
  tuple val(meta), path(bam), path(stats), path(peaks)
  path(fripScoreHeader)

  output:
  path("*_FRiP.tsv"), emit: fripTsv
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo "BEDtools"\$(intersectBed 2>&1 | grep "Version" | cut -f2 -d:) > versions.txt
  READS_IN_PEAKS=\$(intersectBed -a ${bam} -b ${peaks} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  peak_type=\$(echo ${peaks} | cut -f2 -d.)
  name=\$(basename ${peaks} .bed)
  grep 'mapped (' $stats | awk -v a="\$READS_IN_PEAKS" -v peakType="\$peak_type" '{printf "${prefix}_"peakType"\\t%.2f\\n", a/\$1}' | cat $fripScoreHeader - > "\$name"_FRiP.tsv
  """
}

