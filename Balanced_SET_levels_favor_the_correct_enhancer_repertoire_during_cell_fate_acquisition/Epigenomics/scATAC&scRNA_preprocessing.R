library(ArchR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(plyranges)

#set Number of threads
addArchRThreads(threads = 30)

#Reference genome

addArchRGenome("mm10")


