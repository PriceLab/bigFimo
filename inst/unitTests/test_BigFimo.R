library(RUnit)
library(MotifDb)
#----------------------------------------------------------------------------------------------------
library(BigFimo)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
# create a new meme file:
# library(MotifDb)
# motifs <- query(MotifDb, c("sapiens","jaspar2022"))
# length(motifs) # [1] 692
# export(motifs, con="~/github/bigFimo/jaspar2022-human.meme", format="meme")
# modiffy ~/github/bigFimo/helpers/fimoProcess-hg38.R:
# if(nrow(tbl.regions) > 0){
#    meme.file <- "~/github/bigFimo/jaspar2022-human.meme"
#    tbl.fimo <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
#    } #
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()

   test_subdivideRegion_1()
   test_subdivideRegion_2()
   test_subdivideRegion_3()
   test_subdivideRegion_10()
   test_subdivideBigRegion_30()

   test_singleSmallRegionUsage()
   test_singleHugeRegionUsage()
    
} # runTests
#---------------------------------------------------------------------------------------------------
# BigFimo originally (and optionally, still) depended on genehancer to provide genomic regions with
# reported transcriptional relations with the targetGene.  another option, which could be inovked at
# the same time, specifies known open chromatin regions though these options are still supported,
# and ctor tests can be found in old_test_BigFimo.R, our usage has evolved towards calling FIMO only
# on an entire contiguous genomic region (perhaps +/- 1M) around a gene of interest.  we make
# subsequent use of trena's FeatureTable to add in genehancer and open chromatin data, annotating
# corresponding fimo-identified regions with genehancer & open chromatin scores
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    targetGene <- "ZBTB7A"

    if(file.exists(targetGene)){
       unlink(targetGene, recursive=TRUE)
       }

       # hg38, targetGene's chomLoc plus some upstream
    chrom <- "chr19"
    start <- 4040801
    end <- 4077194

      #-----------------------------------------------------------------
      # partition the target region into 3 slightly overlapping parts
      # overlap here is 60, with duplicates removed in later processing
      #-----------------------------------------------------------------

    count <- 3    # number of processes, hence also number of genomic regions
    half.overlap <- 30
    landmarks <- seq(from=start, to=end, length.out=count+1)
    starts <- landmarks[1:count]
    ends <- landmarks[2:(count+1)]
    starts <- starts - half.overlap
    ends <- ends + half.overlap
    tbl.roi <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)

    motifs <- query(MotifDb, c("sapiens"), c("jaspar2022"))

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=tbl.roi,  # not actually open chromatin - this includes all bases
                       processCount=count,
                       fimoThreshold=1e-6,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=NA,
                       chrom=chrom, start=start, end=end,
                       motifs=motifs,
                       meme.file.path=file.path(targetGene, "motifs.meme"))

   filenames.roi <- bf$createFimoTables()
   bf$createMemeFile()

   checkEquals(is(bf), "BigFimo")
   tbl.gh <- bf$get.tbl.gh()
   checkEquals(nrow(tbl.gh), 0)

} # test_ctor
#---------------------------------------------------------------------------------------------------
# typical use proceeds by dividing a large genomic region into slightly overlapping
# subregions, running fimo independently on each region, then combining and uniquing
# the results.   test here our ability to reliably and reproducibly subdivide 
# a big region
test_subdivideRegion_1 <- function()
{
    message(sprintf("--- test_subdivideRegion_1"))

         #------------------------------------------------
         # just one small region, no subdivision to start
         #------------------------------------------------

    chrom <- "chr1"
    start <- 100
    end <- 200
    number.of.regions <- 1
    overlap <- 30
    half.overlap <- as.integer(overlap/2)

    tbl.roi <- subdivideGenomicRegion(chrom, start, end, 1, overlap)
    checkEquals(tbl.roi$chrom, chrom)
    checkEquals(tbl.roi$start, start - half.overlap)
    checkEquals(tbl.roi$end, end + half.overlap)
    checkEquals(tbl.roi$size, (2 * half.overlap) + 1 + end - start)

} # test_subdivideRegion_1
#----------------------------------------------------------------------------------------------------
test_subdivideRegion_2 <- function()
{
    message(sprintf("--- test_subdivideRegion_2"))

         #------------------------------------------------
         # just one small region, 2 subdivisions
         #------------------------------------------------

    chrom <- "chr1"
    start <- 100
    end <- 200
    number.of.regions <- 2
    overlap <- 30
    half.overlap <- as.integer(overlap/2)

    tbl.roi <-
        subdivideGenomicRegion(chrom, start, end, number.of.regions, overlap)

    checkEquals(nrow(tbl.roi), number.of.regions)
    checkEquals(tbl.roi$chrom, rep(chrom, number.of.regions))

        #------------------------------------------------------------
        # use reduce(GRanges): is the entire ranged plus shoulders, 
        # covered in one row?
        #------------------------------------------------------------

    tbl.reduced <- as.data.frame(reduce(GRanges(tbl.roi)))
    checkEquals(nrow(tbl.reduced), 1)
    expected.span <- (2 * half.overlap) + 1 + end - start
    checkEquals(tbl.reduced$width, expected.span)
    checkEquals(tbl.roi$size, c(81, 81))

} # test_subdivideRegion_2
#---------------------------------------------------------------------------------------------------
test_subdivideRegion_3 <- function()
{
    message(sprintf("--- test_subdivideRegion_3"))

         #------------------------------------------------
         # just one small region, 3 subdivisions
         #------------------------------------------------

    chrom <- "chr1"
    start <- 100
    end <- 200
    number.of.regions <- 3
    overlap <- 30
    half.overlap <- as.integer(overlap/2)

    tbl.roi <-
        subdivideGenomicRegion(chrom, start, end, number.of.regions, overlap)

    checkEquals(nrow(tbl.roi), number.of.regions)
    checkEquals(tbl.roi$chrom, rep(chrom, number.of.regions))

        #------------------------------------------------------------
        # use reduce(GRanges): is the entire ranged plus shoulders, 
        # covered in one row?
        #------------------------------------------------------------

    tbl.reduced <- as.data.frame(reduce(GRanges(tbl.roi)))
    checkEquals(nrow(tbl.reduced), 1)
    expected.span <- (2 * half.overlap) + 1 + end - start
    checkEquals(tbl.reduced$width, expected.span)
    checkEquals(tbl.roi$size, c(64, 65, 64))

} # test_subdivideRegion_3
#---------------------------------------------------------------------------------------------------
test_subdivideRegion_10 <- function()
{
    message(sprintf("--- test_subdivideRegion_10"))

    chrom <- "chr1"
    start <- 100
    end <- 200
    number.of.regions <- 10
    overlap <- 0
    half.overlap <- as.integer(overlap/2)

    tbl.roi <-
        subdivideGenomicRegion(chrom, start, end, number.of.regions, overlap)

    checkEquals(nrow(tbl.roi), number.of.regions)
    checkEquals(tbl.roi$chrom, rep(chrom, number.of.regions))

        #------------------------------------------------------------
        # use reduce(GRanges): is the entire ranged plus shoulders, 
        # covered in one row?
        #------------------------------------------------------------

    tbl.reduced <- as.data.frame(reduce(GRanges(tbl.roi)))
    checkEquals(nrow(tbl.reduced), 1)
    expected.span <- (2 * half.overlap) + 1 + end - start
    checkEquals(tbl.reduced$width, expected.span)
    checkEquals(tbl.roi$size, rep(11, number.of.regions))

} # test_subdivideRegion_10
#---------------------------------------------------------------------------------------------------
test_subdivideBigRegion_30 <- function()
{
    message(sprintf("--- test_subdivideBigRegion_30"))

    targetGene <- "ZBTB7A"   # not actually used here, just for explanation
    chrom <- "chr19"
    tss <- 4066900
    start <- tss - 1e6
    end   <- tss + 1e6
    number.of.regions <- 30
    overlap <- 100
    half.overlap <- as.integer(overlap/2)

    tbl.roi <-
        subdivideGenomicRegion(chrom, start, end, number.of.regions, overlap)

    checkEquals(nrow(tbl.roi), number.of.regions)
    checkEquals(tbl.roi$chrom, rep(chrom, number.of.regions))

        #------------------------------------------------------------
        # use reduce(GRanges): is the entire ranged plus shoulders, 
        # covered in one row?
        #------------------------------------------------------------

    tbl.reduced <- as.data.frame(reduce(GRanges(tbl.roi)))
    checkEquals(nrow(tbl.reduced), 1)
    expected.span <- (2 * half.overlap) + 1 + end - start
    checkEquals(tbl.reduced$width, expected.span)
    region.sizes <- tbl.roi$size
    checkTrue(all(region.sizes >= 66767 & region.sizes <= 66768))

} # test_subdivideBigRegion_30
#---------------------------------------------------------------------------------------------------
test_singleSmallRegionUsage <- function()
{
    message(sprintf("--- test_singleSmallRegionUsage"))

    targetGene <- "ZBTB7A"

    if(file.exists(targetGene)){
       unlink(targetGene, recursive=TRUE)
       }

       # hg38, targetGene's chomLoc plus some upstream
    chrom <- "chr19"
    start <- 4040801
    end   <- 4077194
    processCount <- 3    # number of processes, hence also number of genomic regions
    number.of.regions <- processCount

      #-----------------------------------------------------------------
      # partition the target region into 3 slightly overlapping parts
      # overlap here is 60, with duplicates removed in later processing
      #-----------------------------------------------------------------
    overlap <- 100
    half.overlap <- as.integer(overlap/2)
    tbl.roi <-
        subdivideGenomicRegion(chrom, start, end, number.of.regions, overlap)

       # quick sanity check: just one row, span slightly larger than requested
    tbl.reduced <- as.data.frame(reduce(GRanges(tbl.roi)))
    checkEqualsNumeric(tbl.reduced$width/(1+end-start), 1.0, tol=1e-2)
    checkTrue(tbl.reduced$width/(1+end-start) > 1)
    
    motifs <- query(MotifDb, c("sapiens"), c("jaspar2022"))
    meme.file <- file.path(targetGene, "motifs.meme")
    meme.file <- "motifs.meme"
    rtracklayer::export(motifs, con=meme.file, format="meme")

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=tbl.roi,  # not actually open chromatin - this includes all bases
                       processCount=processCount,
                       fimoThreshold=1e-6, 
                      use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=NA,
                       chrom=chrom, start=start, end=end,
                       motifs=motifs,
                       meme.file.path=meme.file)

   filenames.roi <- bf$createFimoTables()
   #bf$createMemeFile()
   checkTrue(file.exists(meme.file))

   checkEquals(length(filenames.roi), processCount)
      # make sure these RData files, which will be read by a script that directly
      # runs FIMO, match the regions calculated above
   for(i in seq_len(length(filenames.roi))){
       tbl.r <- get(load(file.path(targetGene, filenames.roi[i])))
       checkEquals(tbl.r, tbl.roi[i,])
       }
   bf$runMany()
   bf$waitForCompletion(sleepInterval=1)

   fimo.output.files.by.region <-
              list.files(path=targetGene,
                         pattern=sprintf("^fimo.%s.*RData", targetGene))

   checkEquals(length(fimo.output.files.by.region), processCount)
   tbl.fimo <- bf$combineResults()
   checkEquals(colnames(tbl.fimo), c("chrom", "start", "end", "tf", "strand", "score",
                                     "p.value", "matched_sequence", "motif_id"))
   checkTrue(nrow(tbl.fimo) > 400)
   checkTrue(nrow(tbl.fimo) < 500)

       #--------------------------------------------------
       # make sure that FIMO matches start and end very close to
       # the requested chromosomal coordinates
       #--------------------------------------------------
   checkTrue(abs(min(tbl.fimo$start) - start) < 100)
   checkTrue(abs(max(tbl.fimo$end) - end) < 100)

} # test_singleSmallRegionUsage
#----------------------------------------------------------------------------------------------------
test_singleHugeRegionUsage <- function()
{
    message(sprintf("--- test_singleHUGERegionUsage"))

    targetGene <- "NDUFS2"

    if(file.exists(targetGene)){
       unlink(targetGene, recursive=TRUE)
       }

       # hg38, targetGene's chomLoc plus some upstream
    chrom <- "chr1"
    start <- 160660231
    end   <- 161753478

      #-----------------------------------------------------------------
      # partition the target region into 3 slightly overlapping parts
      # overlap here is 60, with duplicates removed in later processing
      #-----------------------------------------------------------------

    processCount <- 20    # number of processes, hence also number of genomic regions
    number.of.regions <- processCount
    overlap <- 60
    half.overlap <- as.integer(overlap/2)

    tbl.roi <-
        subdivideGenomicRegion(chrom, start, end, number.of.regions, overlap)

    tbl.reduced <- as.data.frame(reduce(GRanges(tbl.roi)))
    checkEquals(nrow(tbl.reduced), 1)  # ensures we have continuous coverage
    checkEqualsNumeric(tbl.reduced$width/(1+end-start), 1.0, tol=1e-2)
      # tbl.roi should cover just a bit more than requested
    checkTrue(tbl.reduced$width/(1+end-start) > 1)
    checkTrue(tbl.reduced$width/(1+end-start) < 1.01)

    motifs <- query(MotifDb, c("sapiens"), c("jaspar2022"))
    meme.file <- "motifs.meme"
    rtracklayer::export(motifs, con=meme.file, format="meme")

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=tbl.roi,  # not actually open chromatin - this includes all bases
                       processCount=processCount,
                       fimoThreshold=1e-4,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=NA,
                       chrom=chrom, start=start, end=end,
                       motifs=motifs,
                       meme.file.path=meme.file)

   filenames.roi <- bf$createFimoTables()
   checkTrue(file.exists(meme.file))

   checkEquals(length(filenames.roi), processCount)
      # make sure these RData files, which will be read by a script that directly
      # runs FIMO, match the regions calculated above
   for(i in seq_len(length(filenames.roi))){
       tbl.r <- get(load(file.path(targetGene, filenames.roi[i])))
       checkEquals(tbl.r, tbl.roi[i,])
       }
   bf$runMany()
   bf$waitForCompletion(sleepInterval=10)

   fimo.output.files.by.region <-
              list.files(path=targetGene,
                         pattern=sprintf("^fimo.%s.*RData", targetGene))

   checkEquals(length(fimo.output.files.by.region), processCount)
   tbl.fimo <- bf$combineResults()
   checkEquals(colnames(tbl.fimo), c("chrom", "start", "end", "tf", "strand", "score",
                                     "p.value", "matched_sequence", "motif_id"))
   printf("tbl.fimo: %d nrows", nrow(tbl.fimo))
   checkTrue(nrow(tbl.fimo) > 200000)
   checkTrue(nrow(tbl.fimo) < 300000)

       #--------------------------------------------------
       # make sure that FIMO matches start and end very close to
       # the requested chromosomal coordinates
       #--------------------------------------------------
   checkTrue(abs(min(tbl.fimo$start) - start) < 100)
   checkTrue(abs(max(tbl.fimo$end) - end) < 100)

} # test_singleHugeRegionUsage
#----------------------------------------------------------------------------------------------------
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
