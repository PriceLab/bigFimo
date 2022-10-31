library(GenomicRanges)
source("~/github/bigFimo/R/BigFimo.R")

fimoThreshold <- 1e-3     # 1e-4: 1654  1e-6: 32  1e-3  11,331
processCount <- 40
targetGene <- "TSPAN14"

if(file.exists(targetGene)){
   unlink(targetGene, recursive=TRUE)
   }


gh.elite.only <- FALSE
maxGap.between.oc.and.gh <- 5000

   # this region was discovered using viz function below

chrom <- "chr10"
start <- 79902142
end   <- 81008099

landmarks <- seq(from=start, to=end, length.out=processCount)
starts <- landmarks[1:(processCount-1)]
ends <- landmarks[2:processCount]
starts <- starts - 100
ends <- ends + 100

  #
tbl.oc.faux <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)
stopifnot(width(reduce(GRanges(tbl.oc.faux))) == 1106158)

bf <-  BigFimo$new(targetGene,
                   tbl.oc=tbl.oc.faux,
                   processCount=processCount,
                   fimoThreshold=fimoThreshold,
                   use.genehancer=FALSE,
                   gh.elite.only=FALSE,
                   maxGap.between.oc.and.gh=5000,
                   chrom=chrom, start=start-100, end=end+100)
bf$calculateRegionsForFimo()
filenames.roi <- bf$createFimoTables()
checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))
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
printf("tbl.fimo aggregated: %d rows", nrow(tbl.fimo))
checkEquals(ncol(tbl.fimo), 9)
checkTrue(min(tbl.fimo$start) >= (start - 1000))
checkTrue(max(tbl.fimo$end) <= (end + 1000))
filename <- sprintf("%s/%s-%s-%d-%d_%f.RData", targetGene, targetGene, chrom, start, end, fimoThreshold)
save(tbl.fimo, file=filename)
