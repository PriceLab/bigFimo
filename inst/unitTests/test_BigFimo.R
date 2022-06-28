library(RUnit)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(org.Mm.eg.db)
#----------------------------------------------------------------------------------------------------
library(BigFimo)
# source("~/github/bigFimo/R/BigFimo.R")
# source("~/github/endophenotypeExplorer/R/getExpressionMatrices.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_createIndicesToDistributeTasks()

    # test_zbtb7a()

    test_calculateRegionsForFimo_small()
    test_calculateRegionsForFimo_small_PTK2B()
    test_calculateRegionsForFimo_wholeGene_mayoATAC()
    test_calculateRegionsForFimo_medium()
    test_calculateRegionsForFimo_maximal()
    test_includeOnlyGeneHancerIntersectingOC()
    test_createFimoTables_explicitRegion()
    test_tspan14.noGeneHancer()
    test_mouseGene()
    test_b4galt3()

    #test_runMany()   # works, but takes a few minutes


} # runTests
#---------------------------------------------------------------------------------------------------
human.brain.tcx.boca.oc <- function()
{
    dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    f <- "tbl.boca.RData"
    full.path <- file.path(dir, f)
    checkTrue(file.exists(full.path))
    tbl.oc <- get(load(full.path))
    tbl.ctx <- subset(tbl.oc, tissue %in% c("DLPFC", "VLPFC"))
    gr.ctx <- reduce(GRanges(tbl.ctx))
    tbl.oc <- as.data.frame(gr.ctx)
    colnames(tbl.oc)[1] <- "chrom"
    tbl.oc$chrom <- as.character(tbl.oc$chrom)
    dim(tbl.oc)   # 125430

    invisible(tbl.oc)

} # human.brain.tcx.boca.oc
#---------------------------------------------------------------------------------------------------
human.brain.consensus.boca.oc <- function()
{
    dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    f <- "boca-hg38-consensus-ATAC.RData"
    full.path <- file.path(dir, f)
    checkTrue(file.exists(full.path))
    tbl.oc <- get(load(full.path))
    dim(tbl.oc)   # 125430

    invisible(tbl.oc)

} # human.brain.consensus.boca.oc
#---------------------------------------------------------------------------------------------------
human.brain.mayo.oc <- function()
{
   dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   file <- "mayoAllPeaks.merged.96064x4.RData"
   full.path <- file.path(dir, file)
   tbl.atac <- get(load(full.path))
   tbl.oc <- tbl.atac[, c("chrom", "start", "end")]

   invisible(tbl.oc)

} # human.brain.mayo.oc
#---------------------------------------------------------------------------------------------------
human.erythropoiesis.brand.oc <- function()
{
    dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions"
    full.path <- file.path(dir, "tbl.atacMerged.RData")
    checkTrue(file.exists(full.path))
    tbl.oc <- get(load(full.path))

    invisible(tbl.oc)

} # human.erythropoiesis.brand.oc
#---------------------------------------------------------------------------------------------------
mouse.liver.encode.oc <- function()
{
    data.dir <- "~/github/TrenaProjectMouseLC/inst/extdata/genomicRegions"
    filename <- "dnase-liver-8weeks-ENCSR836BNL.RData"
    full.path <- file.path(data.dir, filename)
    checkTrue(file.exists(full.path))
    tbl.oc <- get(load(full.path))

    size <- round(sum(width(GRanges(tbl.oc)))/1000000, digits=2)
    printf("mouse liver oc %s of %5.2fM bases, %d rows", filename, size, nrow(tbl.oc))

    invisible(tbl.oc)

} # mouse.liver.encode.oc
#---------------------------------------------------------------------------------------------------
mouse.brain.encode.oc <- function()
{
    data.dir <- "~/github/TrenaProjectMouseBrain/inst/extdata/genomicRegions"
    filename <- "brainTissueMaleAdult8weeks-ENCSR000COF.RData"
    full.path <- file.path(data.dir, filename)

    checkTrue(file.exists(full.path))
    tbl.oc <- get(load(full.path))
    dim(tbl.oc)
    size <- round(sum(width(GRanges(tbl.oc)))/1000000, digits=2)
    printf("mouse brain oc %s of %5.2fM bases, %d rows", filename, size, nrow(tbl.oc))
    invisible(tbl.oc)

} # mouse.brain.encode.oc
#----------------------------------------------------------------------------------------------------
test_getOpenChromatin <- function()
{
    message(sprintf("--- test_getOpenChromatin"))

       #----------------------------------------
       # consensus hg38 atac-seq from boca
       #----------------------------------------
    tbl.oc <- human.brain.consensus.boca.oc()
    checkTrue(all(colnames(tbl.oc)[1:3] == c("chrom", "start", "end")))
    checkTrue(class(tbl.oc[, "chrom"]) == "character")
    head(tbl.oc)
    checkEquals(length(unique(tbl.oc$chrom)), 24)
    checkTrue(all(grepl("chr", unique(tbl.oc$chrom))))
    checkEquals(nrow(tbl.oc), 125430)

       #----------------------------------------
       # all mayo ATAC.seq
       #----------------------------------------
    tbl.oc <- human.brain.mayo.oc()
    checkTrue(all(colnames(tbl.oc)[1:3] == c("chrom", "start", "end")))
    checkTrue(class(tbl.oc[, "chrom"]) == "character")
    head(tbl.oc)
    checkEquals(length(unique(tbl.oc$chrom)), 25)
    checkTrue(all(grepl("chr", unique(tbl.oc$chrom))))
    checkEquals(nrow(tbl.oc), 96064)

       #----------------------------------------
       # brand lab merged erythropoeis
       #----------------------------------------
    tbl.oc <- human.erythropoiesis.brand.oc()
    checkTrue(all(colnames(tbl.oc)[1:3] == c("chrom", "start", "end")))
    checkTrue(class(tbl.oc[, "chrom"]) == "character")
    head(tbl.oc)
    checkEquals(length(unique(tbl.oc$chrom)), 24)
    checkTrue(all(grepl("chr", unique(tbl.oc$chrom))))
    checkEquals(nrow(tbl.oc), 956694)

       #----------------------------------------
       # encode mouse brain
       #----------------------------------------

    tbl.oc <- mouse.brain.encode.oc()
    checkTrue(all(colnames(tbl.oc)[1:3] == c("chrom", "start", "end")))
    checkTrue(class(tbl.oc[, "chrom"]) == "character")
    checkEquals(length(unique(tbl.oc$chrom)), 21)   # 1-19, X, Y
    checkTrue(all(grepl("chr", unique(tbl.oc$chrom))))
    checkEquals(nrow(tbl.oc), 258549)

       #----------------------------------------
       # encode mouse liver
       #----------------------------------------

    tbl.oc <- mouse.liver.encode.oc()
    checkTrue(all(colnames(tbl.oc)[1:3] == c("chrom", "start", "end")))
    checkTrue(class(tbl.oc[, "chrom"]) == "character")
    checkEquals(length(unique(tbl.oc$chrom)), 21)   # 1-19, X, Y
    checkTrue(all(grepl("chr", unique(tbl.oc$chrom))))
    checkEquals(nrow(tbl.oc), 499693)

} # test_getOpenChromatin
#---------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

        #------------------------------------------------------------
        # BrandLab erythropoiesis merged (all condition) ATAC-seq
        # regions, thus w/o scores
        #------------------------------------------------------------

        # first, without an explicit genomic region
        # BigFimo intersect genehancer with tbl.oc

    targetGene <- "BACH1"

    if(file.exists(targetGene)){
       unlink(sprintf("%s/*", targetGene))
       unlink(targetGene, recursive=TRUE)
       }

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
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
    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-4,
                       gh.elite.only=TRUE,
                       maxGap.between.oc.and.gh=100,
                       chrom=chrom, start=start, end=end)

    checkEquals(is(bf), "BigFimo")
    tbl.gh <- bf$get.tbl.gh()
    checkTrue(nrow(tbl.gh) < 40)

} # test_ctor
#---------------------------------------------------------------------------------------------------
test_zbtb7a <- function()
{
    targetGene <- "ZBTB7A"

    if(file.exists(targetGene)){
       unlink(sprintf("%s/*", targetGene))
       unlink(targetGene, recursive=TRUE)
       }

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=30,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
                       chrom=NA, start=NA, end=NA)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 42)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkTrue(nrow(tbl.gh.oc) > 70)

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 25)
        # note: exactly 2 regions described in each file, for a balanced distribution
        # of tasks in runMany()
    checkEquals(filenames.roi[1],  "ZBTB7A.01.fimoRegions-00003.RData")
    checkEquals(filenames.roi[25], "ZBTB7A.25.fimoRegions-00003.RData")

    actual.processes.needed <- length(bf$getFimoRegionsFileList())
    checkEquals(actual.processes.needed, 25)   # one fewer than requested

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

    if(file.exists(targetGene)){
       unlink(sprintf("%s/*", targetGene))
       unlink(targetGene, recursive=TRUE)
       }

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-6,
                       gh.elite.only=TRUE,
                       maxGap.between.oc.and.gh=5000)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 32)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(73, 5))

} # test_calculateRegionsForFimo_small
#---------------------------------------------------------------------------------------------------
# compare the three hg38 oc sources:
#   human.erythropoiesis.brand.oc()
test_calculateRegionsForFimo_small_PTK2B <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_small_PTK2B"))

    targetGene <- "PTK2B"

    if(file.exists(targetGene)){
       unlink(sprintf("%s/*", targetGene))
       unlink(targetGene, recursive=TRUE)
       }

    chrom <- "chr8"
    start <- 27320895
    end   <- 27330083

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-4,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
                       chrom=chrom, start=start, end=end)

    tbl.gh <- bf$get.tbl.gh()
    dim(tbl.gh)
    checkEquals(nrow(tbl.gh), 4)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    dim(tbl.gh.oc)
    checkEquals(dim(tbl.gh.oc), c(6, 5))

    if(exists("igv")){
      displayGenomicRegion(igv, "chr8:27,309,174-27,347,212")
      track <- DataFrameQuantitativeTrack("gh.fp",
                                          tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                          autoscale=TRUE, color="darkgreen")
      displayTrack(igv, track)
      track <- DataFrameAnnotationTrack("gh.oc.fp", tbl.gh.oc, color="red", trackHeight=23)
      displayTrack(igv, track)
      }

    return(TRUE)

} # test_calculateRegionsForFimo_small_PTK2B
#---------------------------------------------------------------------------------------------------
# we currently have four human sources, three brain, one erythropoiesis.  how do they compare?
test_calculateRegionsForFimo_allHumanSources_NDUFS2 <- function()
{
   message(sprintf("--- test_calculateRegionsForFimo_allHumanSources_NDUFS2"))
   targetGene <- "NDUFS2"
   targetGene <- "BACH1"
   targetGene <- "PPOX"

   if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

   fimo.threshold <- 1e-5

   bf <-  BigFimo$new(targetGene,
                      tbl.oc=human.brain.mayo.oc(),
                      processCount=1,
                      fimoThreshold=fimo.threshold,
                      gh.elite.only=FALSE,
                      maxGap.between.oc.and.gh=5000,
                      ocExpansion=0)
   bf$calculateRegionsForFimo()
   tbl.regions.mayo <- bf$get.tbl.gh.oc()
   dim(tbl.regions.mayo)   # 17 5

   bf <-  BigFimo$new(targetGene,
                      tbl.oc=human.brain.consensus.boca.oc(),
                      processCount=1,
                      fimoThreshold=fimo.threshold,
                      gh.elite.only=FALSE,
                      maxGap.between.oc.and.gh=5000,
                      ocExpansion=0)
   bf$calculateRegionsForFimo()
   tbl.regions.boca <- bf$get.tbl.gh.oc()
   dim(tbl.regions.boca)  # 19 5

   bf <-  BigFimo$new(targetGene,
                      tbl.oc=human.erythropoiesis.brand.oc(),
                      processCount=1,
                      fimoThreshold=fimo.threshold,
                      gh.elite.only=FALSE,
                      maxGap.between.oc.and.gh=5000,
                      ocExpansion=0)
   bf$calculateRegionsForFimo()
   tbl.regions.brand <- bf$get.tbl.gh.oc()
   nrow(tbl.regions.brand)   # 53

   tbl.gh <- bf$get.tbl.gh()

   if(exists("igv")){
     showGenomicRegion(igv, targetGene)
     zoomOut(igv)
     zoomOut(igv)
     track <- DataFrameQuantitativeTrack("gh",
                                         tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                         autoscale=TRUE, color="darkgreen")
     displayTrack(igv, track)
     track <- DataFrameAnnotationTrack("mayo", tbl.regions.mayo, color="random")
     displayTrack(igv, track)
     track <- DataFrameAnnotationTrack("boca", tbl.regions.boca, color="random")
     displayTrack(igv, track)
     track <- DataFrameAnnotationTrack("brand", tbl.regions.brand, color="random")
     displayTrack(igv, track)
     }

    return(TRUE)

} # test_calculateRegionsForFimo_allHumanSources_NDUFS2
#----------------------------------------------------------------------------------------------------
test_calculateRegionsForFimo_wholeGene_mayoATAC <- function()
{
    message(sprintf("--- test_calculateRegionsForFimo_wholeGene_mayoATAC"))

    targetGene <- "PTK2B"

   if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }


    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-4,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000)

    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 40)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(81, 5))

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

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 29075018
    end   <- 29330888
    end - start

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-6,
                       gh.elite.only=TRUE,
                       maxGap.between.oc.and.gh=5000,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 23)
    checkTrue(all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(56, 5))

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

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28995780
    end   <- 29120251
    end - start

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-6,
                       gh.elite.only=TRUE,
                       maxGap.between.oc.and.gh=0,
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

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=2,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=0)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 142)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(102, 5))

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

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }


       # this region was discovered using viz function below
    chrom <- "chr21"
    start <- 28900000
    end   <- 29120251
    end - start

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=3,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 9)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(11, 5))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 4)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

       # with five processes created, at least six fimo regions files are needed

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=5,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
                       chrom=chrom, start=start, end=end)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 9)
    checkTrue(!all(tbl.gh$elite))

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(11, 5))

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

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=3,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
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
    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }


    fimoThreshold <- 1e-6
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

       # this region was discovered using viz function below
    chrom <- NA
    start <- NA
    end   <- NA

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=30,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
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
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=5,
                       fimoThreshold=1e-5,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
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

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

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
                       tbl.oc=human.erythropoiesis.brand.oc(),
                       processCount=3,
                       fimoThreshold=1e-6,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=5000,
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
    checkEquals(dim(tbl.fimo), c(78, 9))
    checkTrue(min(tbl.fimo$start) >= start)
    checkTrue(max(tbl.fimo$end) <= end)

} # test_runMany
#---------------------------------------------------------------------------------------------------
test_mouseGene <- function()
{
    message(sprintf("--- test_mouseGene"))

    targetGene <- "Nfe2l2"

    if(file.exists(targetGene)){
      unlink(sprintf("%s/*", targetGene))
      unlink(targetGene, recursive=TRUE)
      }

    tbl.geneInfo.mm10 <- get(load("~/github/TrenaProjectMM10/inst/extdata/geneInfoTable.RData"))
    tbl.targetGene <- subset(tbl.geneInfo.mm10, geneSymbol==targetGene)
    bf <-  BigFimo$new(targetGene=targetGene,
                       genome="mm10",
                       tbl.oc=mouse.liver.encode.oc(),
                       processCount=2,
                       fimoThreshold=1e-5,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=NA,
                       chrom=tbl.targetGene$chrom,
                       start=tbl.targetGene$start - 100,
                       end=tbl.targetGene$end + 100)
    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 0)

    tbl.regions <- bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(31,4))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 3)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))


    bf$runMany()
    # Sys.sleep(10)
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

} # test_mouseGene
#---------------------------------------------------------------------------------------------------
test_tspan14.noGeneHancer <- function()
{
    message(sprintf("--- test_tspan14.noGeneHancer"))

    fimoThreshold <- 1e-6     # 1e-4: 1654  1e-6: 32  1e-3  11,331
    processCount <- 3
    targetGene <- "TSPAN14"
    gh.elite.only <- FALSE
    maxGap.between.oc.and.gh <- 5000

    if(file.exists(targetGene)){
       unlink(sprintf("%s/*", targetGene))
       unlink(targetGene, recursive=TRUE)
       }

     # this region was discovered using viz function below

    chrom <- "chr10"
    start <- 80508193
    end   <- 80512482

        # create a fake table of open chromatin, the regions we want
        # fimo to search.   this fake table is 3 slightly overlapping
        # regions covering  start to end, with some shoulder
    landmarks <- seq(from=start, to=end, length.out=4)
    starts <- landmarks[1:3]
    ends <- landmarks[2:4]
    starts <- starts -100
    ends <- ends + 100

    tbl.oc.faux <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=tbl.oc.faux,
                       processCount=3,
                       fimoThreshold=fimoThreshold,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                   maxGap.between.oc.and.gh=5000,
                   chrom=chrom, start=start-100, end=end+100)
    bf$calculateRegionsForFimo()
    filenames.roi <- bf$createFimoTables()
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))


    bf$runMany()
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
       # with regions overlapping, there may be some duplicates.  eliminate them
    tbl.fimo <- unique(tbl.fimo[order(tbl.fimo$start, decreasing=FALSE),])
    rownames(tbl.fimo) <- NULL
    checkEquals(ncol(tbl.fimo), 9)
    checkTrue(min(tbl.fimo$start) >= (start - 1000))
    checkTrue(max(tbl.fimo$end) <= (end + 1000))

} # test_tspan14.noGeneHancer
#---------------------------------------------------------------------------------------------------
# B4GALT3 is in an eQTL downstream of ad GWAS variant rs4575098. a coding variant for this gene
# is used in recent oligogenic genomic risk score for LOAD.
# i add this test after some months not using BigFimo, and after making the package installable,
#
test_b4galt3 <- function()
{
    message(sprintf("--- test_b4galt3"))

    chrom <- "chr1"
    start <- 161184710
    end   <- 161186977

    span <- 1 + end - start
    print(span)
    targetGene <- "B4GALT3"

    fimo.pval.threshold <- 1e-6

      # BigFimo needs specified regions, as many as the processCount (or more)
      # which can be fed to different processes to work on
      # we start here without genomic regions created by ATAC-seq of GeneHancer
      # so gin some up.

    landmarks <- seq(from=start, to=end, length.out=4)
    starts <- landmarks[1:3]
    ends <- landmarks[2:4]
    starts <- starts - 10
    ends <- ends + 10

    tbl.oc.faux <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)

    if(file.exists(targetGene))
        unlink(targetGene, recursive=TRUE)

    bf <-  BigFimo$new(targetGene=targetGene,
                       tbl.oc=tbl.oc.faux,
                       processCount=3,
                       fimoThreshold=fimo.pval.threshold,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=0,
                       chrom=chrom, start=start, end=end)

    tbl.gh <- bf$get.tbl.gh()
    checkEquals(nrow(tbl.gh), 0)

    bf$calculateRegionsForFimo()
    tbl.gh.oc <- bf$get.tbl.gh.oc()
    checkEquals(dim(tbl.gh.oc), c(3, 3))

    filenames.roi <- bf$createFimoTables()
    checkEquals(length(filenames.roi), 3)
    checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

      #------------------------------------------------
      # now run fimobatch on those three region files
      #------------------------------------------------

    bf$runMany()
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
       # with regions overlapping, there may be some duplicates.  eliminate them
    tbl.fimo <- unique(tbl.fimo[order(tbl.fimo$start, decreasing=FALSE),])
    rownames(tbl.fimo) <- NULL
    expected.colnames <- c("chrom", "start", "end", "tf", "strand", "score", "p.value",
                           "matched_sequence", "motif_id")
    checkEquals(colnames(tbl.fimo), expected.colnames)
    checkTrue(all(tbl.fimo$start >= start - 10))
    checkTrue(all(tbl.fimo$end   <= end - 10))
    checkTrue(all(tbl.fimo$chrom == chrom))
    checkTrue(all(tbl.fimo$p.value <= fimo.pval.threshold))
    checkEquals(nrow(tbl.fimo), 33)

    viz <- FALSE
    if(viz){
       igv <- start.igv(targetGene, "hg38")
       zoomOut(igv)
       library(ghdb)
       ghdb <- GeneHancerDB()
       tbl.roi <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
       track <- DataFrameAnnotationTrack("roi", tbl.roi, color="darkgreen")
       displayTrack(igv, track)

       track <- DataFrameAnnotationTrack("regions.faux", tbl.oc.faux, color="random")
       displayTrack(igv, track)

       tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
       dim(tbl.gh)
       tbl.gh.strong <- subset(tbl.gh, elite)
       dim(tbl.gh.strong)
       track <- DataFrameQuantitativeTrack("gh",
                                           tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                           autoscale=TRUE, color="brown")
       displayTrack(igv, track)
       track <- DataFrameAnnotationTrack("fimo", tbl.fimo, color="blue")
       displayTrack(igv, track)
       }


} # test_b4galt3
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
