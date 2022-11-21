library(BigFimo)
library(MotifDb)

printf <- function(...) print(noquote(sprintf(...)))

if(interactive()){   # for testing only
    args <- c("chr19",
              "4040801",
              "4041801",
              "ZBTB7A",
              "1",
              "1e-5",
              "motifs.meme",
              "./")
    args <- c("chr1",
              "161200670",
              "161204670",
              "NDUFS2",
              "1",
              "1e-5",
              "human-jaspar2018-hocomoco-swissregulon.meme",
              "./")
   }else{
      args <- commandArgs(trailingOnly=TRUE)
      }

stopifnot(length(args) == 8)
chrom <- args[1]
start <- as.numeric(args[2])
end <- as.numeric(args[3])
targetGene <- args[4]
processCount <- as.numeric(args[5])
fimo.pval.threshold <- as.numeric(args[6])
meme.file.path <- args[7]
output.directory <- args[8]

if(file.exists(targetGene)) {
    unlink(targetGene, recursive=TRUE)
    }            

span <- 1 + end - start
printf("running bigFimo across %5.2fk", span/1000)
tbl.roi <- subdivideGenomicRegion(chrom, start, end, 1, overlap=60)
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
                   chrom=chrom, start=start, end=end,
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
out.filename <- sprintf("tbl.fimo.%s.%s:%d-%d.RData", targetGene, chrom, start, end)
out.filepath <- file.path(output.directory, out.filename)
printf("saving %d hits to %s", nrow(tbl.fimo), out.filepath)
save(tbl.fimo, file=out.filepath)

