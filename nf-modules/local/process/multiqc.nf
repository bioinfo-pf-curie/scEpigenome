/*
 * MultiQC for RNA-seq report
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
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
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path warnings
   //Modules
  path ('star/*')
  path ('index/*')
  //Logs
  //path ("bowtie2/*")
  //path("removeRtPcr/*")
  //path("cellThresholds/*")
  //path("rmDup/*")
  // Weighted histogram
  path ('countUMI/*')

  output:
  path splan, emit: splan
  path "*report.html", emit: report
  path "*_data", emit: data

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_scchip_report" : "--filename scchip_report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  modulesList = "-m custom_content -m star"
  warn = warnings.name == 'warnings.txt' ? "--warn warnings.txt" : ""
  """
  mqc_header.py --splan ${splan} --name "scChIP-seq" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} ${warn} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}
