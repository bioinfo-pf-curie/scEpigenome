include { seqkitReplace } from '../process/seqkitReplace'
include { createBatches } from '../process/createBatches'

// Set the meta.chunk value in case of multiple sequencing lanes
def setMetaChunk(row){
  def map = []
  row[1].eachWithIndex() { file, i ->
    meta = row[0].clone()
    meta.chunk = i+1
    meta.part = row[1].size()
    map += [meta, file]
  }
  return map
}

workflow sccuttagPlateFlow{

  take:
  reads
  batchSize

  main:
  chVersions = Channel.empty()

  seqkitReplace(
    reads
  )
  chVersions = chVersions.mix(seqkitReplace.out.versions)

  createBatches(
   seqkitReplace.out.reads,
   batchSize
  )
  chVersions = chVersions.mix(createBatches.out.versions)

  // group by read pairs
  chPairedFastq = createBatches.out.reads
    .map { meta, fastq -> fastq }
    .flatten()
    .collate(2)

  // Add the batch information
  chFastqBatches = createBatches.out.reads
    .map { meta, fastq -> meta }
    .combine(chPairedFastq)
    .map { meta, r1, r2 -> [meta, [r1,r2]] }
    .groupTuple()
    .flatMap { it -> setMetaChunk(it) }
    .collate(2)

  emit:
  versions = chVersions 
  reads = chFastqBatches
}