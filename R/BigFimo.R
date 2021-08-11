library(R6)
library(RPostgreSQL)
library(ghdb)
source("~/github/fimoService/batchMode/fimoBatchTools.R")
BigFimo = R6Class("BigFimo",
#--------------------------------------------------------------------------------
private = list(targetGene=NULL,
               project=NULL,
               fimoThreshold=NULL,
               gh.elite.only=NULL,
               maxGap.between.oc.and.gh=NULL,
               ocExpansion=NULL,
               tbl.gh=NULL,
               tbl.gh.oc=NULL,
               regionSize=NULL,
               processCount=NULL,
               chromosome=NULL,
               loc.start=NULL,
               loc.end=NULL,
               motifs=NULL,
               fimoRegionsFileList=NULL
               ),
#--------------------------------------------------------------------------------
public = list(

    initialize = function(targetGene, project, processCount, fimoThreshold,
                          gh.elite.only=TRUE, maxGap.between.oc.and.gh=5000,
                          ocExpansion=100, chrom=NA, start=NA, end=NA){
        stopifnot(project %in% c("BrandLabErythropoiesis",
                                 "PriceLabBrainFootprints",
                                 "MayoATAC"))
         if(!file.exists(targetGene))
             dir.create(targetGene)
         private$targetGene  <- targetGene
         private$project <- project
         private$processCount <- processCount
         private$fimoThreshold <- fimoThreshold
         private$gh.elite.only <- gh.elite.only
         private$maxGap.between.oc.and.gh <- maxGap.between.oc.and.gh
         private$ocExpansion <- ocExpansion
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

       if(grepl("hagfish", Sys.info()[["nodename"]])){
          suppressWarnings(
             db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
          if(length(db.access.test) == 0)
              stop("khaleesi postgres server unavailable")
          } # if hagfish
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
    get.tbl.gh.oc = function(){
       if(colnames(private$tbl.gh.oc)[1] == "seqnames"){
           colnames(private$tbl.gh.oc)[1] <- "chrom"
           private$tbl.gh.oc$chrom <- as.character(private$tbl.gh.oc$chrom)
           }
       private$tbl.gh.oc
       },

    #------------------------------------------------------------------------
    createMemeFile = function(){
       meme.file <- "jaspar2018-hocomocoCore.meme"
       private$motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))
       printf("--- exporting %d motifs to %s", length(private$motifs), meme.file)
       export(private$motifs, con=meme.file, format="meme")
       },

    #------------------------------------------------------------------------
    calculateRegionsForFimo = function(){
       tbl.gh <- private$tbl.gh   # convenience
       private$regionSize <- with(tbl.gh, max(end) - min(start))
       printf("   full genehancer region: %dk", round(private$regionSize/1000, digits=0))
       if(private$gh.elite.only)
          tbl.gh <- subset(tbl.gh, elite)
       if(!grepl("chr", tbl.gh$chrom[1]))
          tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)

       gr.gh <- reduce(GRanges(tbl.gh))


         # this shoulder artificially expands the oc hit regions
         # allowing for imprecision in those results.


       shoulder <- private$ocExpansion

       if(private$project == "BrandLabErythropoiesis"){
          dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions"
          full.path <- file.path(dir, "tbl.atacMerged.RData")
          stopifnot(file.exists(full.path))
          tbl.atacMerged <- get(load(full.path))
          tbl.oc <- subset(tbl.atacMerged, chrom==private$chromosome &
                                           start >= (private$loc.start-shoulder) &
                                           end <=   (private$loc.end+shoulder))
          } # brandLabErythropoiesis
       if(private$project == "PriceLabBrainFootprints"){
          db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_16", host="khaleesi")
          # browser()
          query <- sprintf("select * from regions where chrom='%s' and start >= %d and endpos <= %d",
                           private$chromosome, private$loc.start, private$loc.end)

          tbl.oc <- dbGetQuery(db, query)
          dbDisconnect(db)
          if(nrow(tbl.oc) > 0){
             tbl.oc <- tbl.oc[, c("chrom", "start", "endpos")]
             colnames(tbl.oc) <- c("chrom", "start", "end")
             }
          } # PriceLabBrainFootprints

       if(private$project == "MayoATAC"){
          f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"
          stopifnot(file.exists(f))
          tbl.atac <- get(load(f))
          tbl.oc <- tbl.atac[, c("chrom", "start", "end", "max.score")]
          tbl.oc <- subset(tbl.oc, max.score >= 0.001)
          } # MayoATAC

       gr.oc <- GRanges(tbl.oc)

         # identify the oc regions which permissively "overlap" with gh,
         # where any oc region with in maxgap of a gh region is considerd
         # an overlap
       gr.ov <- findOverlaps(gr.oc, gr.gh, maxgap=private$maxGap.between.oc.and.gh)
       gr.oc.near.gh <- reduce(gr.oc[queryHits(gr.ov)])
       tbl.gh.oc <- as.data.frame(gr.oc.near.gh)

         # add a width column, for subsequent convenience
       tbl.gh.oc$width <- 1 + tbl.gh.oc$end - tbl.gh.oc$start
       private$tbl.gh.oc <- tbl.gh.oc
       },

    #------------------------------------------------------------------------
       # create one binary data.frame per process, splitting regions
       # equally among them
    createFimoTables = function(){
       tbl.roi <- self$get.tbl.gh.oc()
       n <- private$processCount
       group.size <-  nrow(tbl.roi) %/% n
       remainder  <-  nrow(tbl.roi) %% n
       indices <- lapply(seq_len(n), function(i) seq(from=(1 + (i-1)*group.size), length.out=group.size))
       filenames <- list()
       for(i in seq_len(length(indices))){
           elements.this.file <- indices[[i]]
           tbl.out <- tbl.roi[elements.this.file,]
           filename <- sprintf("%s.%02d.fimoRegions-%05d.RData", private$targetGene, i, nrow(tbl.out))
           dir <- private$targetGene
           full.path <- file.path(dir, filename)
           filenames[[i]] <- filename
           save(tbl.out, file=full.path)
           }
       leftover.indices <- setdiff(seq_len(nrow(tbl.roi)), unlist(indices))
       if(length(leftover.indices) > 0){
           filename <- sprintf("%s.%02d.fimoRegions-%05d.RData", private$targetGene, i+1, length(leftover.indices))
           tbl.out <- tbl.roi[leftover.indices,]
           dir <- private$targetGene
           full.path <- file.path(dir, filename)
           filenames[[i+1]] <- filename
           save(tbl.out, file=full.path)
           }
       private$fimoRegionsFileList <- unlist(filenames)
       private$fimoRegionsFileList
       },

    #------------------------------------------------------------------------
    getFimoRegionsFileList = function(){
        private$fimoRegionsFileList
        },

    #------------------------------------------------------------------------
    runMany = function(){
       script <- "~/github/bigFimo/R/fimoProcess.R"
       printf("---- starting %d processes", length(private$fimoRegionsFileList)) # private$processCount)
       for(fimoRegionsFile in private$fimoRegionsFileList){
           full.path <- file.path(private$targetGene, fimoRegionsFile)
           stopifnot(file.exists(full.path))
           cmd <- sprintf("Rscript %s %s %s %10.8f",
                          script,
                          private$targetGene, full.path, private$fimoThreshold)
           printf("cmd: %s", cmd)
           system(cmd, wait=FALSE)
           } # for i
       }
    ) # public

) # class BigFimo
#--------------------------------------------------------------------------------
