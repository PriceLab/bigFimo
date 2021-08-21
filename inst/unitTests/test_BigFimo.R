library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/bigFimo/R/BigFimo.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_specifiedRegionCtor()
    test_createIndicesToDistributeTasks()

    # test_zbtb7a()

    test_calculateRegionsForFimo_small()
    test_calculateRegionsForFimo_small_brainFootprints()
    test_calculateRegionsForFimo_small_mayoATAC()
    test_calculateRegionsForFimo_wholeGene_mayoATAC()
    test_calculateRegionsForFimo_medium()
    test_calculateRegionsForFimo_maximal()
    test_includeOnlyGeneHancerIntersectingOC()

    test_createFimoTables_explicitRegion()
    #test_runMany()

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
test_zbtb7a <- function()
{
    targetGene <- "ZBTB7A"
    processCount <- 30
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000
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
    checkEquals(nrow(tbl.gh), 39)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(58, 5))

    filenames.roi <- bf$createFimoTables()
        # note: exactly 2 regions described in each file, for a balanced distribution
        # of tasks in runMany()
    checkEquals(filenames.roi[1],  "ZBTB7A.01.fimoRegions-00002.RData")
    checkEquals(filenames.roi[29], "ZBTB7A.29.fimoRegions-00002.RData")

    actual.processes.needed <- length(bf$getFimoRegionsFileList())
    checkEquals(actual.processes.needed, 29)   # one fewer than requested

    runMany <- FALSE   # will need to adjust the dimensions test of tbl.fimo
    if(runMany){
        bf$runMany()
        Sys.sleep(10)
        completed <- FALSE

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
        checkEquals(dim(tbl.fimo), c(487, 9))
        }

} # test_zbtb7a
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
test_calculateRegionsForFimo_small_mayoATAC <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_small_mayoATAC"))

    targetGene <- "PTK2B"
    processCount <- 2
    fimoThreshold <- 1e-4
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000
       # this region was discovered using viz function below
    chrom <- "chr8"
    start <- 27320895
    end   <- 27330083

    end - start

    bf <-  BigFimo$new(targetGene,
                       project="MayoATAC",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)

    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 3)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(3, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("gh.fp",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.oc.fp", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

    return(TRUE)


} # test_calculateRegionsForFimo_small_mayoATAC
#---------------------------------------------------------------------------------------------------
test_calculateRegionsForFimo_wholeGene_mayoATAC <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_wholeGene_mayoATAC"))

    targetGene <- "PTK2B"
    processCount <- 2
    fimoThreshold <- 1e-4
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000
    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BigFimo$new(targetGene,
                       project="MayoATAC",
                       processCount,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)

    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 37)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(28, 5))

    if(exists("igv")){
      track <- DataFrameQuantitativeTrack("gh.fp",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.oc.fp", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

    return(TRUE)

} # test_calculateRegionsForFimo_wholeGene_mayoATAC
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
test_createFimoTables_explicitRegion <- function()
{
    message(sprintf("--- test_createFimoTables_explicitRegion"))

    targetGene <- "BACH1"
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28900000
    end   <- 29120251
    end - start

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

       # with five processes created, at least six fimo regions files are needed

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount=5,
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
    checkEquals(length(filenames.roi), 6)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

} # test_createFimoTables_explicitRegion
#---------------------------------------------------------------------------------------------------
createIndices <- function(tbl.roi, processCount)
{
    numberOfGroups <- processCount
    if(nrow(tbl.roi) < numberOfGroups)   # not enough regions for the number of processes?
        numberOfGroups <- nrow(tbl.roi)
    remainder  <-  nrow(tbl.roi) %% numberOfGroups
    group.size <-  nrow(tbl.roi) %/% numberOfGroups
    if(remainder > group.size){
       group.size <- group.size + 1
       numberOfGroups <- nrow(tbl.roi) %/% group.size
       remainder <- nrow(tbl.roi) %% numberOfGroups
       }
    indices <- lapply(seq_len(numberOfGroups),
                      function(i) seq(from=(1 + (i-1)*group.size), length.out=group.size))
    indices.created <- sum(unlist(lapply(indices, length)))
    if(remainder > 0)
       indices[[numberOfGroups+1]] <- setdiff(seq_len(nrow(tbl.roi)), unlist(indices))

    indices

} # createIndices
#----------------------------------------------------------------------------------------------------
test_createIndicesToDistributeTasks <- function()
{
       # create a BigFimo object, but its details do not matter.  we
       # use it hear only to render possible a call to bf$createIndicesToDistributeTasks()

    targetGene <- "BACH1"
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount=3,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)

   indices <- bf$createIndicesToDistributeTasks(mtcars, 5)
   checkEquals(length(indices), 6)
   sizes <- unlist(lapply(indices, length))
   checkEquals(sizes, c(6,6,6,6,6,2))
   checkEquals(sum(sizes), nrow(mtcars))

   indices <- bf$createIndicesToDistributeTasks(mtcars, 6)
   checkEquals(length(indices), 7)
   sizes <- unlist(lapply(indices, length))
   checkEquals(sizes, c(5, 5, 5, 5, 5, 5, 2))
   checkEquals(sum(sizes), nrow(mtcars))

   indices <- bf$createIndicesToDistributeTasks(mtcars, 9)
   checkEquals(length(indices), 8)
   sizes <- unlist(lapply(indices, length))
   checkEquals(sizes, c(4, 4, 4, 4, 4, 4, 4, 4))
   checkEquals(sum(sizes), nrow(mtcars))

   indices <- bf$createIndicesToDistributeTasks(mtcars, 20)
   sizes <- unlist(lapply(indices, length))
   checkEquals(sum(sizes), 32)
   checkEquals(length(indices), 16)

   indices <- bf$createIndicesToDistributeTasks(mtcars, 19)
   sizes <- unlist(lapply(indices, length))
   checkEquals(sum(sizes), 32)
   checkEquals(length(indices), 16)

       # a 31 row table spread across 20 or fewer groups
   indices <- bf$createIndicesToDistributeTasks(mtcars[-1,], 20)
   sizes <- unlist(lapply(indices, length))
   checkEquals(sum(sizes), 31)
   checkEquals(length(indices), 16)

} # test_createIndicesToDistributeTasks
#----------------------------------------------------------------------------------------------------
test_createFimoTables_fullGene <- function()
{
    message(sprintf("--- test_createFimoTables_fullGene"))

    targetGene <- "ZBTB7A"
    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount=30,
                       fimoThreshold,
                       gh.elite.only,
                       maxGap.between.oc.and.gh,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 39)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(58, 5))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 31)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

       # with five processes created, at least six fimo regions files are needed

    bf <-  BigFimo$new(targetGene,
                       project="BrandLabErythropoiesis",
                       processCount=5,
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
    checkEquals(length(filenames.roi), 6)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

} # test_createFimoTables_fullGene
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
