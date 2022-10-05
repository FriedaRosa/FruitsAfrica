library(ape)
trees <- read.nexus("TREES.nex")

length(trees)
trees[[1]]

for (i in 1:length(trees)){
  write.nexus(trees[[i]],file=paste0("tree",i))
}
