library(ghdb)
library(GenomicRanges)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_brandLab.combine.atac.and.gh.RCOR1()
    test_brandLab.combine.atac.and.gh.BACH1()

} # runTests
#----------------------------------------------------------------------------------------------------
brandLab.combine.atac.and.gh <- function(targetGene, maxGap.between.atac.and.gh=5000)
{
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")

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
    checkEquals(dim(tbl.rcor1), c(131, 5))
    size <- sum(with(tbl.rcor1, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    checkEqualsNumeric(size.mb, 0.37, tolerance=0.01)

       # maxGap 0

    tbl.rcor1 <- brandLab.combine.atac.and.gh("RCOR1", maxGap.between.atac.and.gh=0)
    dim(tbl.rcor1)
    checkEquals(dim(tbl.rcor1), c(81, 5))

    size <- sum(with(tbl.rcor1, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    size.mb

    checkEqualsNumeric(size.mb, 0.29, tolerance=0.01)

    tbl.rcor1 <- brandLab.combine.atac.and.gh("RCOR1", maxGap.between.atac.and.gh=10000)
    dim(tbl.rcor1)
    checkEquals(dim(tbl.rcor1), c(158, 5))
    size <- sum(with(tbl.rcor1, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    size.mb
    checkEqualsNumeric(size.mb, 0.41, tolerance=0.01)

} # test_brandLab.combine.atac.and.gh.RCOR1
#----------------------------------------------------------------------------------------------------
# empirical, n=1 support for maxGap.between.atac.and.gh=5000, seen in the "viz" block below.
# try with maxGap values 0, 100, 3000, 5000, 10000.  5000 seems to generously pick up all the relevant atac
test_brandLab.combine.atac.and.gh.BACH1 <- function()
{
    message(sprintf("--- test_brandLab.cobmine.atac.and.gh.BACH1"))

    targetGene <- "BACH1"

    tbl <- brandLab.combine.atac.and.gh(targetGene, 5000)
    dim(tbl)
    checkEquals(dim(tbl), c(194, 5))
    size <- sum(with(tbl, 1 + end - start))
    size.mb <- round(size/1000000, digits=2)
    checkEqualsNumeric(size.mb, 0.39, tolerance=0.01)

    viz <- FALSE
    if(viz){
       igv <- start.igv(targetGene)
       zoomOut(igv)
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
