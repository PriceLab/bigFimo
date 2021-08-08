library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/bigFimo/R/BrandLabBigFimo.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_specifiedRegionCtor()
    test_calculateRegionsForFimo_small()
    test_calculateRegionsForFimo_medium()

} # runTests
#---------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.atac.and.gh <- 5000

       # first, without an explicit genomic region

    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=NA, start=NA, end=NA)

    checkEquals(is(bf), "BrandLabBigFimo")
    tbl.gh <- bf$get.tbl.gh()
    checkTrue(nrow(tbl.gh) > 100)
    checkEquals(unique(tbl.gh$gene), targetGene)

       # with an explicit genomic region

    chrom <- "chr21"
    start <- 29297661
    end <- 29300191
    gh.elite.only <- TRUE
    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)

    checkEquals(is(bf), "BrandLabBigFimo")
    tbl.gh <- bf$get.tbl.gh()
    checkTrue(nrow(tbl.gh) < 40)

} # test_ctor
#---------------------------------------------------------------------------------------------------
test_specifiedRegionCtor <- function()
{
    message(sprintf("--- test_specifiedRegionCtor"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- TRUE
    maxGap.between.atac.and.gh <- 5000
    chrom <- "chr21"
    start <- 29297661
    end <- 29300191
    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)

    checkEquals(is(bf), "BrandLabBigFimo")
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 1)
    checkTrue(tbl.gh$elite)
    checkTrue(start >= tbl.gh$start & start <= tbl.gh$end)
    checkTrue(end >= tbl.gh$start & end <= tbl.gh$end)

} # test_specifiedRegionCtor
#---------------------------------------------------------------------------------------------------
test_calculateRegionsForFimo_small <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_small"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- TRUE
    maxGap.between.atac.and.gh <- 5000
       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 29195612
    end   <- 29207473
    end - start

    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 2)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.atac <- bf$get.tbl.gh.atac()
    checkEquals(dim(tbl.gh.atac), c(2, 5))

} # test_calculateRegionsForFimo_small
#---------------------------------------------------------------------------------------------------
test_calculateRegionsForFimo_medium<- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_medium"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- TRUE
    maxGap.between.atac.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 29075018
    end   <- 29330888
    end - start

    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 22)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.atac <- bf$get.tbl.gh.atac()
    checkEquals(dim(tbl.gh.atac), c(21, 5))

} # test_calculateRegionsForFimo_medium
#---------------------------------------------------------------------------------------------------
viz <- function()
{
   igv <- start.igv("BACH1")
   zoomOut(igv);    zoomOut(igv);
   library(ghdb)
   ghdb <- GeneHancerDB()
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, "BACH1", tissues="all")
   dim(tbl.gh)
   tbl.gh.strong <- subset(tbl.gh, elite)
   dim(tbl.gh.strong)

   track <- DataFrameQuantitativeTrack("gh.strong", tbl.gh.strong[, c("chrom", "start", "end", "combinedscore")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)

   dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions"
   full.path <- file.path(dir, "tbl.atacMerged.RData")
   stopifnot(file.exists(full.path))
   tbl.atacMerged <- get(load(full.path))

   roi <- getGenomicRegion(igv)
   tbl.atac.sub <- subset(tbl.atacMerged, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   dim(tbl.atac.sub)

   track <- DataFrameAnnotationTrack("atac", tbl.atac.sub[, c("chrom", "start", "end")],
                                     color="blue", trackHeight=24)
   displayTrack(igv, track)

     # good test region here:   "chr21:29,269,562-29,295,745"
   showGenomicRegion(igv, "chr21:29,269,562-29,295,745")

} # viz
#---------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()

