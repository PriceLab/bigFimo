library(BigFimo)
library(MotifDb)
library(yaml)

printf <- function(...) print(noquote(sprintf(...)))  # basic utility

yaml.input.directory <- "/usr/local/data/inputs"
if(interactive())
    yaml.input.directory <- "./"

yaml.file.path <- file.path(yaml.input.directory, "bigFimo.yaml")
printf("looking for yaml configuration file here: %s", yaml.file.path)
stopifnot(file.exists(yaml.file.path))

args <- yaml.load(readLines(yaml.file.path))


chrom <- args$chrom
start <- args$start
end <- args$end
targetGene <- args$targetGene
processCount <- args$processCount
fimo.pval.threshold <- as.numeric(args$fimo.pval.threshold)
meme.file.path <- args$meme.file.path
output.directory <- args$output.directory

printf("%30s: %s", "chrom", chrom)
printf("%30s: %d", "start", start)
printf("%30s: %d", "end", end)
printf("%30s: %s", "targetGene", targetGene)
printf("%30s: %d", "processCount", processCount)
printf("%30s: %f", "fimo.pval.threshold", fimo.pval.threshold)
printf("%30s: %s", "meme.file.path", meme.file.path)
printf("%30s: %s", "output.directory", output.directory)

stopifnot(file.exists(meme.file.path))
stopifnot(file.exists(output.directory))

if(file.exists(targetGene)) {
    unlink(targetGene, recursive=TRUE)
    }            

span <- 1 + end - start
printf("running bigFimo across %5.2fk", span/1000)
tbl.roi <- subdivideGenomicRegion(chrom, start, end, processCount, overlap=60)
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

