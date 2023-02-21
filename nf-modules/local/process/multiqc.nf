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
  //path warnings
  //Modules
  path ('star/*')
  path ('index/*')
  //Logs
  path ("bowtie2/*")
  path("allDup/*")
  path("cellThresholds/*")
  path("removeWindowDup/*")
  // Weighted histogram
  path ('countUMI/*')
  // macs2 module 
  path ('peakCounts/*')
  // general stat
  path('frip/*')
  // Homer module
  /*path ('peakAnnot/*')
  // deeptools module
  path('deeptoolsReadDistrib/*')
  // general stat
  path('peakSizes/*')*/

  output:
  path splan, emit: splan
  path "*report.html", emit: report
  path "*_data", emit: data

  script:
  splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  rtitle = customRunName ? "--title \"${params.protocol}\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  minReadsPerCellmqc = params.minReadsPerCellmqc ? "--minReadsPerCellmqc ${params.minReadsPerCellmqc}" : ""
  modulesList = "-m custom_content -m star -m bowtie2 -m deeptools -m macs2 -m homer"
  //warn = warnings.name == 'warnings.txt' ? "--warn warnings.txt" : ""
  """
  stat2mqc.sh ${splan} ${params.minReadsPerCellmqc} ${params.protocol}
  mqc_header.py --splan ${splan} --name ${params.protocol} --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts}  > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}
