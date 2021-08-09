library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/fimoService/batchMode/fimoBatchTools.R")
# bach1 promoter: chr21:29,298,789-29,299,657

if(interactive()){
   targetGene <- "BACH1"
   fimoRegionsFile <- "~/github/bigFimo/inst/extdata/BACH1.01.fimoRegions-00008.RData"
   fimo.threshold <- 1e-6
} else {
    args <- commandArgs(trailingOnly=TRUE)
    printf("arg count: %d", length(args))
    stopifnot(length(args) == 3)
    targetGene <- args[1]
    fimoRegionsFile <- args[2]
    fimo.threshold <- as.numeric(args[3])
    }

tbl.regions <- get(load(fimoRegionsFile))

printf("%s, fimo region count: %d  threshold: %20.10f", targetGene, nrow(tbl.regions), fimo.threshold)

if(nrow(tbl.regions) > 0){
   meme.file <- "~/github/bigFimo/jaspar2018-hocomocoCoreA.meme"
   tbl.fimo <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
   } # if section has gh+atac hits

printf("fimo hits: %d", nrow(tbl.fimo))
date.string <- gsub(" +", "-", date())
section.chrom <- tbl.regions$chrom[1]
section.start <- min(tbl.regions$start)
section.end <- max(tbl.regions$end)
filename <- sprintf("fimo.%s.%s.%d.%d.%s.RData", targetGene, section.chrom, section.start, section.end, date.string)
if(!file.exists(targetGene))
   dir.create(targetGene)
path <- file.path(targetGene, filename)
save(tbl.fimo, file=path)

