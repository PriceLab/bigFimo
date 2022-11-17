# BigFimo
#----------------------------------------------------------------------------------------------------
#' @import R6
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import GenomicRanges
#' @import ghdb
#'
#' @title BigFimo
#------------------------------------------------------------------------------------------------------------------------
#' @name BigFimo
#' @rdname BigFimo
#' @aliases BigFimo
#----------------------------------------------------------------------------------------------------
#' @description
#' An R6 Class to prepare for, and then run, multiple instances of FIMO.
#'
#' @export
#'
#----------------------------------------------------------------------------------------------------
if(Sys.info()[["nodename"]] %in% c("hagfish.local", "khaleesi.systemsbiology.net")){
    source("~/github/fimoService/batchMode/fimoBatchTools.R")
} else { # assume this is a properly configured docker
    source("/usr/local/scripts/fimoBatchTools.R")
    }

#----------------------------------------------------------------------------------------------------
BigFimo = R6Class("BigFimo",
    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   genome=NULL,
                   fimoThreshold=NULL,
                   use.genehancer=NULL,
                   gh.elite.only=NULL,
                   maxGap.between.oc.and.gh=NULL,
                   ocExpansion=NULL,
                   tbl.gh=NULL,
                   tbl.oc=NULL,
                   tbl.gh.oc=NULL,
                   regionSize=NULL,
                   processCount=NULL,
                   chromosome=NULL,
                   loc.start=NULL,
                   loc.end=NULL,
                   motifs=NULL,
                   meme.file.path=NULL,
                   fimoRegionsFileList=NULL
                   ),
#--------------------------------------------------------------------------------
    public = list(


      #' @description
      #' Create a new BigFimo
      #' @param targetGene  character, the gene of interest.
      #' @param genome character default "hg38"
      #' @param tbl.oc data.frame, specified areas of (additional) open chromatin in which to run fio
      #' @param processCount integer number of parallel processes
      #' @param fimoThreshold numeric, 1e-4 for example
      #' @param use.genehancer logical, use GeneHancer for open chromatin, default TRUE
      #' @param gh.elite.only logical, use only doubly-attested GeneHancer regions, default TRUE
      #' @param maxGap.between.oc.and.gh numeric, 5000 by default
      #' @param ocExpansion numeric, pad reported open chromatin by this amount, default 100
      #' @param chrom character optional fimo region, used rather than fimo & tbl.oc? default NA
      #' @param start numeric optional fimo region, used rather than fimo & tbl.oc? default NA
      #' @param end  numeric optional fimo region, used rather than fimo & tbl.oc? default NA
      #' @param motifs a list, default empty, or motifs obtained from querying MotifDb
      #' @param meme.file.path character, optionally points to a (possibly smaller, custom) meme file,
      #'    in contrast to frequently adequate "~/github/bigFimo/jaspar2022-human.meme"
      #' @return A new `BigFimo` object.
    initialize = function(targetGene, genome="hg38", tbl.oc, processCount, fimoThreshold,
                          use.genehancer=TRUE, gh.elite.only=TRUE,
                          maxGap.between.oc.and.gh=5000, ocExpansion=100,
                          chrom=NA, start=NA, end=NA, motifs=list(),
                          meme.file.path=NA){
         if(!file.exists(targetGene))
             dir.create(targetGene)
         private$targetGene  <- targetGene
         private$genome <- genome
         stopifnot(all(colnames(tbl.oc)[1:3] == c("chrom", "start", "end")))
         stopifnot(substr(colnames(tbl.oc)[1], 1, 3) == "chr")
         stopifnot(class(tbl.oc[, "chrom"]) == "character")
         private$tbl.oc <- tbl.oc
         private$processCount <- processCount
         private$fimoThreshold <- fimoThreshold
         private$use.genehancer <- use.genehancer
         private$gh.elite.only <- gh.elite.only
         private$maxGap.between.oc.and.gh <- maxGap.between.oc.and.gh
         private$ocExpansion <- ocExpansion
         private$chromosome=chrom
         private$loc.start=start
         private$loc.end=end
         private$motifs <- motifs
         private$meme.file.path <- meme.file.path

         private$tbl.gh <- self$queryGeneHancer()

         if(!is.na(private$loc.start)){
            if(private$use.genehancer){
                gr.explicitLoc <- GRanges(data.frame(chrom=private$tbl.gh$chrom[1],
                                                     start=private$loc.start,
                                                     end=private$loc.end))
              gr.gh  <- GRanges(private$tbl.gh)
              gh.regions <- subjectHits(findOverlaps(gr.explicitLoc, gr.gh))
              private$tbl.gh <- private$tbl.gh[gh.regions,]
              private$chromosome <- private$tbl.gh$chrom[1]
              private$loc.start <- min(private$tbl.gh$start) - 1000
              private$loc.end <- max(private$tbl.gh$end) + 1000
              } # use.genehancer
            if(!private$use.genehancer){
               tbl.ov <- as.data.frame(findOverlaps(GRanges(private$tbl.oc),
                                                    GRanges(seqnames=private$chromosome,
                                                            IRanges(start=private$loc.start,
                                                                    end=private$loc.end))))
               tbl.oc.sub <- private$tbl.oc[tbl.ov$queryHits,]
               private$tbl.oc <- tbl.oc.sub
               }
           } #
         }, # initialize

    #------------------------------------------------------------------------
      #' @description
      #' call the genehancer database on khaleesi
      #' @return data.frame
    queryGeneHancer = function(){

       if(!private$use.genehancer)
            return(data.frame())
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
      #' @description
      #' access to the previously extracted genehancer table
      #' @return data.frame
    get.tbl.gh = function(){
       if(nrow(private$tbl.gh) == 0) return(private$tbl.gh)
       if(colnames(private$tbl.gh)[1] == "seqnames"){
           colnames(private$tbl.gh)[1] <- "chrom"
           private$tbl.gh$chrom <- as.character(private$tbl.gh$chrom)
           }
        private$tbl.gh
        },

    #------------------------------------------------------------------------
      #' @description
      #' access to the previously calculated overlap of genehancer and open chromatin tables
      #' @return data.frame
    get.tbl.gh.oc = function(){
       if(nrow(private$tbl.gh.oc) == 0) return(data.frame())
       if(colnames(private$tbl.gh.oc)[1] == "seqnames"){
           colnames(private$tbl.gh.oc)[1] <- "chrom"
           private$tbl.gh.oc$chrom <- as.character(private$tbl.gh.oc$chrom)
           }
       private$tbl.gh.oc
       },

    #------------------------------------------------------------------------
      #' @description
      #' fimo needs a meme file with motifs of interest.  just jaspar2022 for now
      #' used to be (before june 2022)    jaspar2018 & hocomoco-core-A
      #  are our standard choice
      #' @return nothing
    createMemeFile = function(){
       if(is.null(private$meme.file.path)){
          private$motifs <- query(MotifDb, c("sapiens"), c("jaspar2022"))
          private$meme.file.path < "bigfimo.meme"  # our default, reused each time
          message(sprintf("--- exporting %d motifs to %s", length(private$motifs),
                          private$meme.file.path))
          export(private$motifs, con=private$meme.file.path, format="meme")
          } # if no meme.file.path provided
       },  # createMemeFile

    #------------------------------------------------------------------------
      #' @description
      #' find the overlap of oc with genehancer - but very generously
      #' any open chromatin within < maxGap.between.oc.and.gh (often 5000)
      #'  will be included.
      #' @return data.frame the overlap of open chromatin with genehancer
    calculateRegionsForFimo = function(){

       if(!private$use.genehancer){
          private$tbl.gh.oc <- private$tbl.oc
          return(private$tbl.gh.oc)
          }
       tbl.gh <- private$tbl.gh   # convenience
       private$regionSize <- with(tbl.gh, max(end) - min(start))
       message(sprintf("   full genehancer region: %dk", round(private$regionSize/1000, digits=0)))
       if(private$gh.elite.only)
          tbl.gh <- subset(tbl.gh, elite)
       if(!grepl("chr", tbl.gh$chrom[1]))
          tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)

       gr.gh <- reduce(GRanges(tbl.gh))

         # this shoulder artificially expands the oc hit regions
         # allowing for imprecision in those results.

       shoulder <- private$ocExpansion
       gr.oc <- GRanges(private$tbl.oc)
       gr.oc <- gr.oc + shoulder   # expands both up and downstream

         # identify the oc regions which permissively "overlap" with gh,
         # where any oc region within maxgap of a gh region is considerd
         # an overlap
       gr.ov <- findOverlaps(gr.oc, gr.gh, maxgap=private$maxGap.between.oc.and.gh)
       gr.oc.near.gh <- reduce(gr.oc[queryHits(gr.ov)])
       tbl.gh.oc <- as.data.frame(gr.oc.near.gh)

         # add a width column, for subsequent convenience
       tbl.gh.oc$width <- 1 + tbl.gh.oc$end - tbl.gh.oc$start
       private$tbl.gh.oc <- tbl.gh.oc
       tbl.gh.oc
       },

    #------------------------------------------------------------------------
      #' @description
      #' we want to partition the fimo regions equally across the approximate
      #' number of suggested processes.  this function returns a list of indices
      #' each a vector. sometimes an extra process is needed, sometimes fewer
      #' @param tbl data.frame, specifies regions we want to distribute across
      #'             multiple processes
      #' @param suggestedProcessCount integer, may be adjusted up or down
      #' @return list of integer vectors
    createIndicesToDistributeTasks = function(tbl, suggestedProcessCount){

       numberOfGroups <- suggestedProcessCount
       if(nrow(tbl) < numberOfGroups)   # not enough regions for the number of processes?
           numberOfGroups <- nrow(tbl)
       remainder  <-  nrow(tbl) %% numberOfGroups
       group.size <-  nrow(tbl) %/% numberOfGroups
       if(remainder > group.size){
           group.size <- group.size + 1
           numberOfGroups <- nrow(tbl) %/% group.size
           remainder <- nrow(tbl) %% numberOfGroups
           }
       indices <- lapply(seq_len(numberOfGroups),
                         function(i) seq(from=(1 + (i-1)*group.size), length.out=group.size))
       indices.created <- sum(unlist(lapply(indices, length)))
       if(remainder > 0)
           indices[[numberOfGroups+1]] <- setdiff(seq_len(nrow(tbl)), unlist(indices))

       indices
       },

    #------------------------------------------------------------------------
      #' @description
      #' create one binary data.frame per process, splitting regions equally among them
      #' @return character vector, the list of files containing the data.frames
      #
    createFimoTables = function(){
       if(private$use.genehancer)
          tbl.roi <- self$get.tbl.gh.oc()
       else
          tbl.roi <- private$tbl.oc
       indices <- self$createIndicesToDistributeTasks(tbl.roi, private$processCount)
       filenames <- list()
       for(i in seq_len(length(indices))){
           elements.this.file <- indices[[i]]
           tbl.out <- tbl.roi[elements.this.file,]
           tbl.out$start <- as.integer(tbl.out$start)
           tbl.out$end <- as.integer(tbl.out$end)
           tbl.out$chrom <- as.character(tbl.out$chrom)
           filename <- sprintf("%s.%02d.fimoRegions-%05d.RData", private$targetGene, i, nrow(tbl.out))
           dir <- private$targetGene
           full.path <- file.path(dir, filename)
           filenames[[i]] <- filename
           save(tbl.out, file=full.path)
           }
       private$fimoRegionsFileList <- unlist(filenames)
       return(private$fimoRegionsFileList)
       },

    #------------------------------------------------------------------------
      #' @description
      #' the files listed here specify regions in which fimo will search
      #' @return character vector
    getFimoRegionsFileList = function(){
        private$fimoRegionsFileList
        },

    #------------------------------------------------------------------------
      #' @description
      #' initiates one fimo processes per region file, each previously created
      #' @return nothing
    runMany = function(){
       hg38.script <- system.file(package="BigFimo", "helpers/fimoProcess-hg38.R")
       mm10.script <- system.file(package="BigFimo", "helpers/fimoProcess-mm10.R")
       script <-  switch(private$genome,
                         hg38=hg38.script,
                         mm10=mm10.script)
       message(sprintf("---- starting %d processes", length(private$fimoRegionsFileList)))
       for(fimoRegionsFile in private$fimoRegionsFileList){
           full.path <- file.path(private$targetGene, fimoRegionsFile)
           stopifnot(file.exists(full.path))
           cmd <- sprintf("Rscript %s %s %s %10.8f %s",
                          script,
                          private$targetGene,
                          full.path,
                          private$fimoThreshold,
                          private$meme.file.path
                          )
           message(sprintf("cmd: %s", cmd))
           system(cmd, wait=FALSE)
           } # for i
         message(sprintf("--- leaving runMany"))
         },
    #------------------------------------------------------------------------
      #' @description
      #' the output file count tells us when all runMany processes have completed
      #' @param sleepInterval numeric the number of seconds to wait before next file count check
      #' @return logical
    waitForCompletion = function(sleepInterval=60){
       done <- FALSE
       while(!done){
          fimo.output.files.by.region <-
              list.files(path=private$targetGene,
                         pattern=sprintf("^fimo.%s.*RData", private$targetGene))
          file.count <- length(fimo.output.files.by.region)
          message(sprintf("--- processes now complete: %d/%d",
                          file.count, private$processCount))
          if(file.count == private$processCount){
             message(sprintf("%d BigFimo processes complete", file.count))
             done <- TRUE
           } else {
               Sys.sleep(sleepInterval)
               }
          } # while !done
       return(TRUE)
       }, # waitForCompletion

    #------------------------------------------------------------------------
      #' @description
      #' runMany creates one output fimo results file for each region, which
      #' we combine into a single possibly large data.frame here
      #' @return data.frame
    combineResults = function(){
       fimo.output.files <- list.files(path=private$targetGene,
                                       pattern=sprintf("^fimo.%s.*RData", private$targetGene))
       stopifnot(length(fimo.output.files) == private$processCount)
       tbls <- list()
       f <- 0
       for(file in fimo.output.files){
          f <- f + 1
          tbl.fimo <- get(load(file.path(private$targetGene, file)))
          tbls[[f]] <- tbl.fimo
          } # for
       tbl.fimo <- do.call(rbind, tbls)
       new.order <- order(tbl.fimo$start, decreasing=FALSE)
       tbl.fimo <- tbl.fimo[new.order,]
       tbl.fimo <- unique(tbl.fimo)
       rownames(tbl.fimo) <- NULL
       invisible(tbl.fimo)
       } # combineResults
    #------------------------------------------------------------------------
    ) # public

) # class BigFimo
#--------------------------------------------------------------------------------
#' given chrom, start, end, and a number of regions, subdivide with overlaps
#'
#' @param chrom character chromosome
#' @param start numeric of region
#' @param end numeric of region
#' @param count numeric the number to calculate
#' @param overlap numeric of neighboring regions, also used to extend start and end
#'
#' @return a data.frame (chrom, start, end, size) of approximately equal count regions
#'
#' @export
#'
#' @aliases subdivideGenomicRegion
#' @rdname subdivideGenomicRegion
#'
subdivideGenomicRegion <- function(chrom, start, end, count, overlap)
{
    
    half.overlap <- as.integer(overlap/2)  # put half backwards, half forward
    landmarks <- seq(from=start, to=end, length.out=count+1)
    starts <- round(landmarks[1:count])
    ends <- round(landmarks[2:(count+1)])
    starts <- starts - half.overlap
    ends <- ends + half.overlap
    tbl.roi <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)
    tbl.roi$size <- with(tbl.roi, 1 + end - start)

    tbl.roi

} # subdivideGenomicRegion
#----------------------------------------------------------------------------------------------------

