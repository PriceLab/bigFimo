library(RUnit)
#----------------------------------------------------------------------------------------------------
if(Sys.info()[["nodename"]] %in% c("hagfish.local", "khaleesi.systemsbiology.net")){
    source("~/github/fimoService/batchMode/fimoBatchTools.R")
} else { # assume this is a properly configured docker
    source("/usr/local/scripts/fimoBatchTools.R")
    }

args <- commandArgs(trailingOnly=TRUE)
printf("arg count: %d", length(args))
stopifnot(length(args) >= 3)
targetGene <- args[1]
fimoRegionsFile <- args[2]
fimo.threshold <- as.numeric(args[3])
meme.file <- "~/github/bigFimo/jaspar2022-human.meme"
if(length(args) == 4)
   meme.file <- args[4]

printf("--- starting ~/github/bigFimo/R/fimoProcess.R")
printf("    fimoRegionsFile: %s", fimoRegionsFile)
printf("    file.exists: %s", file.exists(fimoRegionsFile))
stopifnot(file.exists(fimoRegionsFile))

tbl.regions <- get(load(fimoRegionsFile))
printf("     %s, fimo region count: %d  threshold: %20.10f", targetGene, nrow(tbl.regions), fimo.threshold)

if(nrow(tbl.regions) > 0){
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

