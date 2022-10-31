library(BigFimo)
library(MotifDb)
library(ghdb)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)
chrom <- args[1]
start <- as.numeric(args[2])
end <- as.numeric(args[3])
targetGene <- args[4]

if(file.exists(targetGene)) {
    message(sprintf("%s directory already exists, exiting", targetGene))
    quit("no")
    }            

span <- 1 + end - start
printf("running bigFimo across %5.2fk", span/1000)

landmarks <- seq(from=start, to=end, length.out=31)
starts <- landmarks[1:30]
ends <- landmarks[2:31]
starts <- starts - 30
ends <- ends + 30

tbl.oc.faux <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)

processCount <- 30
fimo.pval.threshold <- 1e-3
gh.elite.only <- FALSE
maxGap.between.oc.and.gh <- 0
motifs <- query(MotifDb, c("sapiens"), c("jaspar2022"))
printf("--- using %d motifs", length(motifs))

bf <-  BigFimo$new(targetGene=targetGene,
                   tbl.oc=tbl.oc.faux,
                   processCount=processCount,
                   fimoThreshold=fimo.pval.threshold,
                   use.genehancer=FALSE,
                   gh.elite.only=FALSE,
                   maxGap.between.oc.and.gh=0,
                   chrom=chrom, start=start, end=end,
                   motifs=motifs)

bf$createMemeFile()
bf$calculateRegionsForFimo()
tbl.gh.oc <- bf$get.tbl.gh.oc()

filenames.roi <- bf$createFimoTables()
length(filenames.roi)
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

# result.files <- list.files(path=targetGene, pattern="^fimo.*")
# tbls <- list()
# for(file in result.files){
#     tbl <- get(load(file.path(targetGene, file)))
#     tbls[[file]] <- tbl
#     }
# tbl.fimo <- do.call(rbind, tbls)
# tbl.fimo <- tbl.fimo[order(tbl.fimo$start, decreasing=FALSE),]
# rownames(tbl.fimo) <- NULL

