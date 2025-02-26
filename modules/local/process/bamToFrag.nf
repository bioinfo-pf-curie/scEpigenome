/*
 * Bam to fragment files files
 */

process bamToFrag {
  tag "${meta.id}"
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai)
  
  output:
  tuple val(meta), path("*tsv.gz"), path("*tbi"), emit: tsv
  path("*.log"), emit: log
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  bamToFrag.py ${args} --input ${bam} 2> ${prefix}_bam2frag.log | sort -T ${params.tmpDir} --parallel=${task.cpus} -k1,1V -k2,2n > ${prefix}.fragments.tsv
  bgzip -@ ${task.cpus} ${prefix}.fragments.tsv
  tabix -p bed ${prefix}.fragments.tsv.gz

  echo \$(python --version 2>&1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}
