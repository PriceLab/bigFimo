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
    test_calculateRegionsForFimo_maximal()
    test_includeOnlyGeneHancerIntersectingAtac()

    test_createFimoTables()

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
    checkEquals(dim(tbl.gh.atac), c(3, 5))

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
    checkEquals(dim(tbl.gh.atac), c(62, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("new.gh",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.atac", tbl.gh.atac, color="red", trackHeight=23)
      displayTrack(igv, track)
      }


} # test_calculateRegionsForFimo_medium
#---------------------------------------------------------------------------------------------------
test_includeOnlyGeneHancerIntersectingAtac <- function()
{
    message(sprintf("--- test_includeOnlyGeneHancerIntersectingAtac"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- TRUE
    maxGap.between.atac.and.gh <- 0

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
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
    checkEquals(dim(tbl.gh.atac), c(3, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("new.gh",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.atac", tbl.gh.atac, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

} # test_includeOnlyGeneHancerIntersectingAtac
#----------------------------------------------------------------------------------------------------
# maximal in these ways:  all genehancer regions for BACH1, elite and not
test_calculateRegionsForFimo_maximal <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_maximal"))

    targetGene <- "BACH1"
    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.atac.and.gh <- 0

       # this region was discovered using viz function below
    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 144)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.atac <- bf$get.tbl.gh.atac()
    checkEquals(dim(tbl.gh.atac), c(114, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("new.gh",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.atac", tbl.gh.atac, color="red", trackHeight=23)
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
    maxGap.between.atac.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
    end - start

    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 11)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.atac <- bf$get.tbl.gh.atac()
    checkEquals(dim(tbl.gh.atac), c(16, 5))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 2)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

       # with three processes created, at least three fimo regions files are needed
       # with 16 rows in tbl.gh.atac, and 16/3 dividing with remainder 1, we
       # should get 4 files


    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount=3,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 11)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.atac <- bf$get.tbl.gh.atac()
    checkEquals(dim(tbl.gh.atac), c(16, 5))

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

    processCount <- 2
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.atac.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
    end - start

    bf <-  BrandLabBigFimo$new(targetGene,
                               processCount,
                               fimoThreshold,
                               gh.elite.only,
                               maxGap.between.atac.and.gh,
                               chrom=chrom, start=start, end=end)
    bf$calculateRegionsForFimo()
    filenames.roi <- bf$createFimoTables()
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

    bf$runMany()
    Sys.sleep(10)
    completed <- FALSE
    while(!completed){
        file.count <- length(list.files(path="BACH1", pattern="^fimo.*"))
        completed <- (file.count == processCount)
        if(!completed){
            printf("waiting for completion: %d/%d", file.count, processCount)
            Sys.sleep(3)
        }
    } # while

    printf("complete %d/%d", processCount, processCount)

    result.files <- list.files(path=targetGene, pattern="^fimo.*")
    checkEquals(length(result.files), processCount)

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

