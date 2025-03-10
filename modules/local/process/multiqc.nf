/*
 * MultiQC for RNA-seq report
 * External parameters :
 */

process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'lowMem'

  input:
  val customRunName
  path splan
  path metadata
  path multiqcConfig
  path ('barcodes/*')
  path ('barcodeCounts/*')
  path ("weightedHist/*")
  path ("stats/*")
  path ("duplicates/*")
  //path("cellThresholds/*")
  //path("removeWindowDup/*")
  path ('peaks/macs2/*')
  path ('peaks/frip/*')
  path ('peaks/qc/*')
  //path ('peaks/size/*')
  path('deeptools/*')
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path warnings

  output:
  path splan, emit: splan
  path "*report.html", emit: report
  path "*_data", emit: data

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : "--title \"${params.protocol}\""
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename scepi_report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  modulesList = "-m custom_content -m star -m bowtie2 -m deeptools -m macs2 -m homer -m samtools"
  warn = warnings.name == 'warnings.txt' ? "--warn warnings.txt" : ""
  minReads = params.protocol == "scchip_indrop" ? 1000 : params.protocol == "sccuttag_indrop" ? 500 : 100

  """
  stat2mqc.sh -s ${splan} -p ${params.protocol} -t ${minReads}
  mqc_header.py --splan ${splan} --name "scEpigenomic" --version ${workflow.manifest.version} ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}
