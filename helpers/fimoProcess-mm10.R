library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/fimoService/batchMode/fimoBatchTools.R")
# bach1 promoter: chr21:29,298,789-29,299,657

if(interactive()){
   targetGene <- "Nfe2l2"
   fimoRegionsFile <- "~/github/bigFimo/inst/extdata/Nfe2l2.31.fimoRegions-00013.RData"
   fimo.threshold <- 1e-6
} else {
    args <- commandArgs(trailingOnly=TRUE)
    printf("arg count: %d", length(args))
    stopifnot(length(args) == 3)
    targetGene <- args[1]
    fimoRegionsFile <- args[2]
    fimo.threshold <- as.numeric(args[3])
    }

printf("--- staring ~/github/bigFimo/R/fimoProcess.R")
printf("    fimoRegionsFile: %s", fimoRegionsFile)
printf("    file.exists: %s", file.exists(fimoRegionsFile))
stopifnot(file.exists(fimoRegionsFile))

tbl.regions <- get(load(fimoRegionsFile))
printf("     %s, fimo region count: %d  threshold: %20.10f", targetGene, nrow(tbl.regions), fimo.threshold)

if(nrow(tbl.regions) > 0){
   meme.file <- "~/github/bigFimo/jaspar2018-hocomocov10-mouse.meme"
   tbl.fimo <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="mm10", pwmFile=meme.file)
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

