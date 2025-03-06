include { seqkitReplace } from '../process/seqkitReplace'
include { createBatches } from '../process/createBatches'

// Create batches of paired-end data
def splitByPairs(row){
  def map = []
  int chunk_nb = 1
  int part = row[1].size()/2
  for (i=0; i<row[1].size(); i+=2) { 
    meta = row[0].clone()
    meta.chunk = chunk_nb
    meta.part = part
    r1 = row[1][i]
    r2 = row[1][i+1]
    map += [meta, [r1,r2]]
    chunk_nb+=1
  }
  return map
}

workflow scepigenomePlateFlow{

  take:
  reads
  batchSize
  sampleDescitpion

  main:
  chVersions = Channel.empty()

  seqkitReplace(
    reads,
    sampleDescitpion
  )
  chVersions = chVersions.mix(seqkitReplace.out.versions)

  createBatches(
   seqkitReplace.out.reads,
   batchSize
  )

  // group by read pairs
  chPairedFastq = createBatches.out.reads
    .flatMap { it -> splitByPairs(it) }
    .collate(2)

  emit:
  versions = chVersions 
  reads = chPairedFastq
}