library(ghdb)
library(GenomicRanges)
library(RUnit)
ghdb <- GeneHancerDB()
#----------------------------------------------------------------------------------------------------
runFimoProcessTests <- function()
{
    test_brandLab.combine.atac.and.gh.RCOR1()
    test_brandLab.combine.atac.and.gh.BACH1()

} # runFimoProcessTests
#----------------------------------------------------------------------------------------------------
brandLab.combine.atac.and.gh <- function(targetGene, maxGap.between.atac.and.gh=5000)
{
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")

      # reduce wide span by using only gh reports from two sources
    tbl.gh <- subset(tbl.gh, elite)

    if(!grepl("chr", tbl.gh$chrom[1]))
        tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)

    gh.chrom <- tbl.gh$chrom[1]
    gh.start <- min(tbl.gh$start)
    gh.end   <- max(tbl.gh$end)
    tbl.gh$width <- 1 + tbl.gh$end - tbl.gh$start
    round(sum(tbl.gh$width)/1000000, digits=3)  # 0.256M
    nrow(tbl.gh)
    gr.gh <- reduce(GRanges(tbl.gh))
    length(gr.gh)

    tbl.atacMerged <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
    shoulder <- 100
    tbl.atac <- subset(tbl.atacMerged, chrom==gh.chrom & start >= (gh.start-shoulder) & end <= (gh.end+shoulder))
    dim(tbl.atac)

    tbl.atac.expanded <- tbl.atac
    shoulder <- 500
    tbl.atac.expanded$start <- tbl.atac.expanded$start - shoulder
    tbl.atac.expanded$end <- tbl.atac.expanded$end + shoulder

    gr.atac <- reduce(GRanges(tbl.atac.expanded))

      # only interestd in atac within 5kb of gh
    tbl.ov <- as.data.frame(findOverlaps(gr.atac, gr.gh, maxgap=maxGap.between.atac.and.gh))

    gr.atac.near.gh <- reduce(gr.atac[tbl.ov[,1],])
    gr.atac.gh <- reduce(c(gr.atac.near.gh, gr.gh))

    tbl.atac.gh <- as.data.frame(gr.atac.gh)
    tbl.atac.gh$width <- 1 + tbl.atac.gh$end - tbl.atac.gh$start
    round(sum(tbl.atac.gh$width)/1000000, digits=3)
    nrow(tbl.atac.gh)
    invisible(tbl.atac.gh)

} # brandLab.combine.atac.and.gh
#----------------------------------------------------------------------------------------------------
test_brandLab.combine.atac.and.gh.RCOR1 <- function()
{
    message(sprintf("--- test_brandLab.cobmine.atac.and.gh.RCOR1"))

    tbl.rcor1 <- brandLab.combine.atac.and.gh("RCOR1")
    checkEquals(dim(tbl.rcor1), c(52, 5))
    size <- sum(with(tbl.rcor1, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    checkEqualsNumeric(size.mb, 0.16, tolerance=0.01)

       # maxGap 0

    tbl.rcor1 <- brandLab.combine.atac.and.gh("RCOR1", maxGap.between.atac.and.gh=0)
    dim(tbl.rcor1)
    checkEquals(dim(tbl.rcor1), c(39, 5))

    size <- sum(with(tbl.rcor1, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    size.mb

    checkEqualsNumeric(size.mb, 0.15, tolerance=0.01)

    tbl.rcor1 <- brandLab.combine.atac.and.gh("RCOR1", maxGap.between.atac.and.gh=10000)
    dim(tbl.rcor1)
    checkEquals(dim(tbl.rcor1), c(52, 5))
    size <- sum(with(tbl.rcor1, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    size.mb
    checkEqualsNumeric(size.mb, 0.16, tolerance=0.01)

} # test_brandLab.combine.atac.and.gh.RCOR1
#----------------------------------------------------------------------------------------------------
# empirical, n=1 support for maxGap.between.atac.and.gh=5000, seen in the "viz" block below.
# try with maxGap values 0, 100, 3000, 5000, 10000.  5000 seems to generously pick up all the relevant atac
test_brandLab.combine.atac.and.gh.BACH1 <- function()
{
    message(sprintf("--- test_brandLab.cobmine.atac.and.gh.BACH1"))

    targetGene <- "BACH1"

    tbl <- brandLab.combine.atac.and.gh(targetGene, 50)
    dim(tbl)
    checkEquals(dim(tbl), c(31, 5))
    size <- sum(with(tbl, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    size.mb
    checkEqualsNumeric(size.mb, 0.10, tolerance=0.01)

    viz <- FALSE
    if(viz){
       igv <- start.igv(targetGene)
       zoomOut(igv);       zoomOut(igv);
       tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
       if(!grepl("chr", tbl.gh$chrom[1]))
           tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)

       tbl.track <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
       #tbl.track$combinedscore <- asinh(tbl.track$combinedscore)
       track <- DataFrameQuantitativeTrack("gh", tbl.track, autoscale=TRUE, color="brown")
       displayTrack(igv, track)

       tbl.atacMerged <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
       tbl.atac <- subset(tbl.atacMerged, chrom==tbl.gh$chrom[1] & start > min(tbl.gh$start) & end < max(tbl.gh$end))
       dim(tbl.atac)

       track <- DataFrameAnnotationTrack("atac", tbl.atac, col="blue")
       displayTrack(igv, track)

       track <- DataFrameAnnotationTrack("atac+gh", tbl, col="darkGreen")
       displayTrack(igv, track)
       }

} # test_brandLab.combine.atac.and.gh.BACH1
#----------------------------------------------------------------------------------------------------
source("~/github/fimoService/batchMode/fimoBatchTools.R")
# bach1 promoter: chr21:29,298,789-29,299,657

if(interactive()){
   targetGene <- "BACH1"
   section.chrom <- "chr21"
   section.start <- 29304270    # much smaller than the region identified by gh+atac
   section.end   <- 29320980
   section.size <- 1 + section.end - section.start
   fimo.threshold <- 1e-6
} else {
    args <- commandArgs(trailingOnly=TRUE)
    printf("arg count: %d", length(args))
    stopifnot(length(args) == 5)
    targetGene <- args[1]
    section.chrom <- args[2]
    section.start <- as.numeric(args[3])
    section.end <- as.numeric(args[4])
    section.size <- 1 + section.end - section.start
    fimo.threshold <- as.numeric(args[5])
    }

tbl.gh.atac <- brandLab.combine.atac.and.gh(targetGene, 5000)
tbl.oi <- subset(tbl.gh.atac, start >= section.start & end <= section.end)
printf("nrow(tbl.oi): %d", nrow(tbl.oi))
tbl.fimo <- data.frame()

if(nrow(tbl.oi) > 0){
   colnames(tbl.oi)[1] <- "chrom"
   tbl.oi$chrom <- as.character(tbl.oi$chrom)
   meme.file <- "~/github/bigFimo/jaspar2018-hocomocoCoreA.meme"
   printf("%s:%d-%d  (%d)", section.chrom, section.start, section.end, section.size);
   printf("fimo threshold: %20.10f", fimo.threshold)

   tbl.fimo <- fimoBatch(tbl.oi, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
   } # if section has gh+atac hits

printf("fimo hits: %d", nrow(tbl.fimo))
date.string <- gsub(" +", "-", date())
filename <- sprintf("fimo.%s.%s.%d.%d.%s.RData", targetGene, section.chrom, section.start, section.end, date.string)
if(!file.exists(targetGene))
   dir.create(targetGene)
path <- file.path(targetGene, filename)
save(tbl.fimo, file=path)

