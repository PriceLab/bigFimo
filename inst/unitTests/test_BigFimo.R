library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/bigFimo/R/BigFimo.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_specifiedRegionCtor()

    test_calculateRegionsForFimo_small()
    test_calculateRegionsForFimo_small_brainFootprints()

    test_calculateRegionsForFimo_small_brainFootprints()
    
    test_calculateRegionsForFimo_medium()
    test_calculateRegionsForFimo_maximal()
    test_includeOnlyGeneHancerIntersectingOC()

    test_createFimoTables()
    test_runMany()

} # runTests
#---------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

       # first, without an explicit genomic region

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=NA, start=NA, end=NA)

    checkEquals(is(bf), "BigFimo")
    tbl.gh <- bf$get.tbl.gh()
    checkTrue(nrow(tbl.gh) > 100)
    checkEquals(unique(tbl.gh$gene), targetGene)

       # with an explicit genomic region

    chrom <- "chr21"
    start <- 29297661
    end <- 29300191
    gh.elite.only <- TRUE
    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)

    checkEquals(is(bf), "BigFimo")
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
    maxGap.between.oc.and.gh <- 5000
    chrom <- "chr21"
    start <- 29297661
    end <- 29300191
    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)

    checkEquals(is(bf), "BigFimo")
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
    maxGap.between.oc.and.gh <- 5000
       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 29195612
    end   <- 29207473
    end - start

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 2)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(3, 5))

} # test_calculateRegionsForFimo_small
#---------------------------------------------------------------------------------------------------
test_calculateRegionsForFimo_small_brainFootprints <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_small_brainFootprints"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000
       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 29195612
    end   <- 29207473
    start <- 29178575
    end   <- 29214651

    end - start

    bf <-  BigFimo$new(targetGene,
                       project="PriceLabBrainFootprints",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)

    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 10)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(160, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("gh.fp",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.oc.fp", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

    return(TRUE)


} # test_calculateRegionsForFimo_small_brainFootprints
#---------------------------------------------------------------------------------------------------
test_calculateRegionsForFimo_medium<- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_medium"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- TRUE
    maxGap.between.oc.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 29075018
    end   <- 29330888
    end - start

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 22)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(62, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("new.gh",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.atac", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }


} # test_calculateRegionsForFimo_medium
#---------------------------------------------------------------------------------------------------
test_includeOnlyGeneHancerIntersectingOC <- function()
{
    message(sprintf("--- test_includeOnlyGeneHancerIntersectingOC"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- TRUE
    maxGap.between.oc.and.gh <- 0

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
    end - start

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 2)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(3, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("new.gh",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.oc", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

} # test_includeOnlyGeneHancerIntersectingOC
#----------------------------------------------------------------------------------------------------
# maximal in these ways:  all genehancer regions for BACH1, elite and not
test_calculateRegionsForFimo_maximal <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_maximal"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 0

       # this region was discovered using viz function below
    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 144)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(114, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("new.gh",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.oc", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

} # test_calculateRegionsForFimo_maximal
#---------------------------------------------------------------------------------------------------
test_createFimoTables <- function()
{
    message(sprintf("--- test_createFimoTables"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
    end - start

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 11)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(16, 5))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 2)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

       # with three processes created, at least three fimo regions files are needed
       # with 16 rows in tbl.gh.oc, and 16/3 dividing with remainder 1, we
       # should get 4 files


    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount=3,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 11)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(16, 5))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 4)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))
    filenames.roi

} # test_createFimoTables
#---------------------------------------------------------------------------------------------------
test_runMany <- function()
{
    message(sprintf("--- test_runMany"))
    targetGene <- "BACH1"

    files <- list.files(path=targetGene, pattern="*.RData")
    if(length(files) > 0)
       unlink(file.path(targetGene, files))

    processCount <- 3    # four will actually created, to handle 16 gh.oc regions
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
    end - start

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    bf$calculateRegionsForFimo()
    filenames.roi <- bf$createFimoTables()
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

    bf$runMany()
    Sys.sleep(10)
    completed <- FALSE
    actual.processes.needed <- length(bf$getFimoRegionsFileList())

    while(!completed){
        file.count <- length(list.files(path=targetGene, pattern="^fimo.*"))
        completed <- (file.count == actual.processes.needed)
        if(!completed){
            printf("waiting for completion: %d/%d", file.count, actual.processes.needed)
            Sys.sleep(3)
        }
    } # while

    printf("complete %d/%d", actual.processes.needed, actual.processes.needed)

    result.files <- list.files(path=targetGene, pattern="^fimo.*")
    checkEquals(length(result.files), actual.processes.needed)
    tbls <- list()
    for(file in result.files){
        tbl <- get(load(file.path(targetGene, file)))
        tbls[[file]] <- tbl
        }
    tbl.fimo <- do.call(rbind, tbls)
    tbl.fimo <- tbl.fimo[order(tbl.fimo$start, decreasing=FALSE),]
    rownames(tbl.fimo) <- NULL
    checkEquals(dim(tbl.fimo), c(62, 9))
    checkTrue(min(tbl.fimo$start) >= start)
    checkTrue(max(tbl.fimo$end) <= end)

} # test_runMany
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

   track <- DataFrameQuantitativeTrack("gh.elite", tbl.gh.strong[, c("chrom", "start", "end", "combinedscore")],
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
