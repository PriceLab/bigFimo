library(R6)
library(ghdb)
source("~/github/fimoService/batchMode/fimoBatchTools.R")
BigFimo = R6Class("DemoApp",
#--------------------------------------------------------------------------------
private = list(targetGene=NULL,
               fimoThreshold=NULL,
               tbl.gh=NULL,
               processCount=NULL,
               chromosome=NULL,
               loc.start=NULL,
               loc.end=NULL,
               motifs=NULL
               ),
#--------------------------------------------------------------------------------
public = list(
     initialize = function(targetGene, fimoThreshold, processCount, chrom=NA, start=NA, end=NA){
         if(!file.exists(targetGene))
             dir.create(targetGene)
         private$targetGene=targetGene
         private$fimoThreshold=fimoThreshold
         private$processCount=processCount
         if(is.na(chrom)){
            ghdb <- GeneHancerDB()
            private$tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
            private$chromosome <- private$tbl.gh$chrom[1]
            private$loc.start <- min(private$tbl.gh$start) - 1000
            private$loc.end <- max(private$tbl.gh$end) + 1000
         } else {
            private$chromosome <- chrom
            private$loc.start <- start
            private$loc.end <- end
            }
         meme.file <- "jaspar2018-hocomocoCore.meme"
         private$motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))
         export(private$motifs, con=meme.file, format="meme")
         },
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

