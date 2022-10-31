library(RUnit)
library(MotifDb)
#----------------------------------------------------------------------------------------------------
library(BigFimo)
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
    test_singleBigRegionUsage()

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
test_singleBigRegionUsage <- function()
{
    message(sprintf("--- test_singleBigRegionUsage"))

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
    meme.file <- file.path(targetGene, "motifs.meme")

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=tbl.roi,  # not actually open chromatin - this includes all bases
                       processCount=count,
                       fimoThreshold=1e-6,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=NA,
                       chrom=chrom, start=start, end=end,
                       motifs=motifs,
                       meme.file.path=meme.file)

   filenames.roi <- bf$createFimoTables()
   bf$createMemeFile()
   checkTrue(file.exists(meme.file))

   checkEquals(length(filenames.roi), 3)
      # make sure these RData files, which will be read by a script that directly
      # runs FIMO, match the regions calculated above
   for(i in seq_len(length(filenames.roi))){
       tbl.r <- get(load(file.path(targetGene, filenames.roi[i])))
       checkEquals(tbl.r, tbl.roi[i,])
       }
   bf$runMany()
   Sys.sleep(30)
   fimo.output.files.by.region <- list.files(path=targetGene, pattern=sprintf("^%s.*RData", targetGene))
   checkEquals(length(fimo.output.files.by.region), 3)
   tbl.fimo <- bf$combineResults()
   checkEquals(colnames(tbl.fimo), c("chrom", "start", "end", "tf", "strand", "score",
                                     "p.value", "matched_sequence", "motif_id"))
   checkTrue(nrow(tbl.fimo) > 400)
   checkTrue(nrow(tbl.fimo) < 500)

} # test_singleBigRegionUsage
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
