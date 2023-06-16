library(BigFimo)
library(MotifDb)

printf <- function(...) print(noquote(sprintf(...)))

if(interactive()){   # for testing only
    args <- c("pilraMicro.csv",
              "PILRA",
              "5",
              "1e-5",
              "human-jaspar2022-hocomoco-core-A.meme",
              "./")
   }else{
      args <- commandArgs(trailingOnly=TRUE)
      }

stopifnot(length(args) == 6)
chromLocFile <- args[1]
stopifnot(file.exists(chromLocFile))
targetGene <- args[2]
processCount <- as.numeric(args[3])
fimo.pval.threshold <- as.numeric(args[4])
meme.file.path <- args[5]
output.directory <- args[6]

if(file.exists(targetGene)) {
    unlink(targetGene, recursive=TRUE)
    }            


tbl.roi <- read.table(chromLocFile, sep="\t", as.is=TRUE)
if(colnames(tbl.roi)[1] != "chrom")
    colnames(tbl.roi) <- c("chrom", "start", "end")

stopifnot(file.exists(meme.file.path))
gh.elite.only <- FALSE
maxGap.between.oc.and.gh <- 0

bf <-  BigFimo$new(targetGene=targetGene,
                   tbl.oc=tbl.roi,
                   processCount=processCount,
                   fimoThreshold=fimo.pval.threshold,
                   use.genehancer=FALSE,
                   gh.elite.only=FALSE,
                   maxGap.between.oc.and.gh=0,
                   chrom=tbl.roi$chrom[1],
                   start=min(tbl.roi$start),
                   end=max(tbl.roi$end),
                   meme.file.path=meme.file.path)

filenames.roi <- bf$createFimoTables()
length(filenames.roi)
stopifnot(file.exists(meme.file.path))

bf$runMany()
bf$waitForCompletion(sleepInterval=10)

fimo.output.files.by.region <-
           list.files(path=targetGene,
                      pattern=sprintf("^fimo.%s.*RData", targetGene))
tbl.fimo <- bf$combineResults()
printf("tbl.fimo: %d nrows", nrow(tbl.fimo))
out.filename <- sprintf("tbl.fimo.%s.%s.RData", targetGene, chromLocFile)
out.filepath <- file.path(output.directory, out.filename)
printf("saving %d hits to %s", nrow(tbl.fimo), out.filepath)
save(tbl.fimo, file=out.filepath)

