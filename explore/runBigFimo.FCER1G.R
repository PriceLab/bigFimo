library(BigFimo)
targetGene <- "FCER1G"

# files <- list.files(path=targetGene, pattern="*")
# if(length(files) > 0)
#    unlink(file.path(targetGene, files))

# library(ghdb)
# ghdb <- GeneHancerDB()
# tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
# start <- min(tbl.gh$start)
# end   <- max(tbl.gh$end)
# span <- 1 + end - start
# printf("gh span for %s: %5.2fk", targetGene, span/1000)  # 808k
# FCER1G
chrom <- "chr1"
start <-  160726012
end <- 161534001

span <- 1 + end - start
printf("running bigFimo across %5.2fk", span/1000)

landmarks <- seq(from=start, to=end, length.out=41)
starts <- landmarks[1:40]
ends <- landmarks[2:41]
starts <- starts - 100
ends <- ends + 100

tbl.oc.faux <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)

processCount <- 30
fimo.pval.threshold <- 1e-3
gh.elite.only <- FALSE
maxGap.between.oc.and.gh <- 0

bf <-  BigFimo$new(targetGene=targetGene,
                   tbl.oc=tbl.oc.faux,
                   processCount=processCount,
                   fimoThreshold=fimo.pval.threshold,
                   use.genehancer=FALSE,
                   gh.elite.only=FALSE,
                   maxGap.between.oc.and.gh=0,
                   chrom=chrom, start=start, end=end)


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

