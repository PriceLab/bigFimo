library(R6)
library(ghdb)
source("~/github/fimoService/batchMode/fimoBatchTools.R")
BrandLabBigFimo = R6Class("BrandLabBigFimo",
#--------------------------------------------------------------------------------
private = list(targetGene=NULL,
               fimoThreshold=NULL,
               gh.elite.only=NULL,
               maxGap.between.atac.and.gh=NULL,
               tbl.gh=NULL,
               tbl.gh.atac=NULL,
               regionSize=NULL,
               processCount=NULL,
               chromosome=NULL,
               loc.start=NULL,
               loc.end=NULL,
               motifs=NULL
               ),
#--------------------------------------------------------------------------------
public = list(

    initialize = function(targetGene, processCount, fimoThreshold, gh.elite.only=TRUE,
                          maxGap.between.atac.and.gh=5000, chrom=NA, start=NA, end=NA){
         if(!file.exists(targetGene))
             dir.create(targetGene)
         private$targetGene  <- targetGene
         private$processCount <- processCount
         private$fimoThreshold <- fimoThreshold
         private$gh.elite.only <- gh.elite.only
         private$maxGap.between.atac.and.gh <- maxGap.between.atac.and.gh
         private$chromosome=chrom
         private$loc.start=start
         private$loc.end=end
         private$tbl.gh <- self$queryGeneHancer()

         if(!is.na(private$loc.start)){
            gr.explicitLoc <- GRanges(data.frame(chrom=private$tbl.gh$chrom[1],
                                                 start=private$loc.start,
                                                 end=private$loc.end))
            gr.gh  <- GRanges(private$tbl.gh)
            gh.regions <- subjectHits(findOverlaps(gr.explicitLoc, gr.gh))
            private$tbl.gh <- private$tbl.gh[gh.regions,]
         } else {
            private$chromosome <- private$tbl.gh$chrom[1]
            private$loc.start <- min(private$tbl.gh$start) - 1000
            private$loc.end <- max(private$tbl.gh$end) + 1000
            }
         },

    #------------------------------------------------------------------------
    queryGeneHancer = function(){
       suppressWarnings(
          db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
       if(length(db.access.test) == 0)
           stop("khaleesi postgres server unavailable")
       ghdb <- GeneHancerDB()
       tbl <- retrieveEnhancersFromDatabase(ghdb, private$targetGene, tissues="all")
       if(private$gh.elite.only)
          tbl <- subset(tbl, elite)
       return(tbl)
       },

    #------------------------------------------------------------------------
    get.tbl.gh = function(){
       if(colnames(private$tbl.gh)[1] == "seqnames"){
           colnames(private$tbl.gh)[1] <- "chrom"
           private$tbl.gh$chrom <- as.character(private$tbl.gh$chrom)
           }
        private$tbl.gh
        },

    #------------------------------------------------------------------------
    get.tbl.gh.atac = function(){
       if(colnames(private$tbl.gh.atac)[1] == "seqnames"){
           colnames(private$tbl.gh.atac)[1] <- "chrom"
           private$tbl.gh.atac$chrom <- as.character(private$tbl.gh.atac$chrom)
           }
       private$tbl.gh.atac
       },

    #------------------------------------------------------------------------
    createMemeFile = function(){
       meme.file <- "jaspar2018-hocomocoCore.meme"
       private$motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))
       printf("--- exporting %d motifs to %s", length(private$motifs), meme.file)
       export(private$motifs, con=meme.file, format="meme")
       },

    #------------------------------------------------------------------------
    setup = function(){
      },

    #------------------------------------------------------------------------
    calculateRegionsForFimo = function(){
       tbl.gh <- private$tbl.gh   # convenience
       private$regionSize <- with(tbl.gh, max(end) - min(start))
       printf("--- full genehancer region: %dk", round(private$regionSize/1000, digits=0))
       if(private$gh.elite.only)
          tbl.gh <- subset(tbl.gh, elite)
       if(!grepl("chr", tbl.gh$chrom[1]))
          tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)

       gr.gh <- reduce(GRanges(tbl.gh))

       dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions"
       full.path <- file.path(dir, "tbl.atacMerged.RData")
       stopifnot(file.exists(full.path))
       tbl.atacMerged <- get(load(full.path))

         # this shoulder artificially expands the atac hit regions
         # allowing for imprecision in those results.

       shoulder <- 100
       tbl.atac <- subset(tbl.atacMerged, chrom==private$chromosome &
                                          start >= (private$loc.start-shoulder) &
                                          end <=   (private$loc.end+shoulder))
       gr.atac <- GRanges(tbl.atac)

         # identify the atac regions which permissively "overlap" with gh,
         # where any atac region with in maxgap of a gh region is considerd
         # an overlap
       gr.ov <- findOverlaps(gr.atac, gr.gh,
                             maxgap=private$maxGap.between.atac.and.gh)
       gr.atac.near.gh <- reduce(gr.atac[subjectHits(gr.ov)])
       tbl.gh.atac <- as.data.frame(gr.atac.near.gh)

         # add a width column, for subsequent convenience
       tbl.gh.atac$width <- 1 + tbl.gh.atac$end - tbl.gh.atac$start
       private$tbl.gh.atac <- tbl.gh.atac
       },

    #------------------------------------------------------------------------
    runMany = function(){
           total.span <- 1 + private$loc.end - private$loc.start
           size <- as.numeric(round(1 + (total.span/private$processCount)))
           script <- "~/github/bigFimo/R/fimoProcess.R"
           printf("---- starting %d processes", private$processCount)
           for(i in seq_len(private$processCount)){
               start <- private$loc.start + ((i-1) * size)
               end <- start + size + 20
               cmd <- sprintf("Rscript %s %s %s %d %d %10.8f",
                              script,
                              private$targetGene, private$chromosome, start, end, private$fimoThreshold)
               printf("cmd: %s", cmd)
               system(cmd, wait=FALSE)
               } # for i
           }
    ) # public

) # class BigFimo
#--------------------------------------------------------------------------------

