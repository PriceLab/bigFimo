args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
targetGene <- args[1]

source("~/github/bigFimo/R/BigFimo.R")

files <- list.files(path=targetGene, pattern="*.RData")

if(length(files) > 0)
   unlink(file.path(targetGene, files))

processCount <- 30
fimoThreshold <- 1e-3
gh.elite.only <- FALSE
maxGap.between.atac.and.gh <- 5000

   # this region was discovered using viz function below
chrom <- NA
start <- NA
end   <- NA

bf <-  BigFimo$new(targetGene,
                   project="MayoATAC",
                   processCount,
                   fimoThreshold,
                   gh.elite.only,
                   maxGap.between.atac.and.gh,
                   chrom=chrom, start=start, end=end)
bf$calculateRegionsForFimo()
filenames.roi <- bf$createFimoTables()
checkTrue(all(file.exists(file.path(targetGene, filenames.roi))))

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
