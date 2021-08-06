library(RUnit)
#----------------------------------------------------------------------------------------------------
source("~/github/bigFimo/R/BigFimo.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_bach1.promoter()
    test_bach1.with.empty.section()

} # runTests
#---------------------------------------------------------------------------------------------------
test_bach1.promoter <- function()
{
    message(sprintf("--- test_bach1.promoter"))

       # delete any results files left over from previous run
    files <- list.files(path="BACH1", pattern="*.RData")
    length(files)
    if(length(files) > 0)
       unlink(file.path("BACH1", files))

    fimoThreshold <- 1e-6
    processCount <- 2

    chrom <- "chr21"
      # a 16.7kb region
    start <- 29304270
    end   <- 29320980
    fimoThreshold <- 1e-6

    runner <- BigFimo$new("BACH1", fimoThreshold, processCount, chrom, start, end)

    runner$runMany()

    Sys.sleep(10)
    completed <- FALSE
    while(!completed){
       file.count <- length(list.files(path="BACH1", pattern="*.RData"))
       completed <- (file.count == processCount)
       if(!completed){
          printf("waiting for completion: %d/%d", file.count, processCount)
          Sys.sleep(3)
          }
       } # while
     printf("complete %d/%d", processCount, processCount)

    files <- list.files(path="BACH1", pattern="*.RData")
    tbls <- list()
    for(i in seq_len(length(files))){
        tbls[[i]] <- get(load(file.path("BACH1", files[i])))
        }

    tbl <- do.call(rbind, tbls)
    sort.order <- order(tbl$start, decreasing=FALSE)
    tbl <- unique(tbl[sort.order,])
    checkEquals(dim(tbl), c(17, 9))

} # test_bach1.promoter
#---------------------------------------------------------------------------------------------------
# if an empty region is encountered, in batch/runMany mode, that process never completes.
# reproduce that here, figure it out
test_bach1.with.empty.section <- function()
{
    message(sprintf("--- test_bach1.with.empty.section"))

    files <- list.files(path="BACH1", pattern="*.RData")
    length(files)
    if(length(files) > 0)
       unlink(file.path("BACH1", files))

    processCount <- 4

    chrom <- "chr21"
      # a 16.7kb region
    start <- 29268880
    end   <- 29282977
    fimoThreshold <- 1e-6

    runner <- BigFimo$new("BACH1", fimoThreshold, processCount, chrom, start, end)

    runner$runMany()

    Sys.sleep(10)
    completed <- FALSE
    while(!completed){
       file.count <- length(list.files(path="BACH1", pattern="*.RData"))
       completed <- (file.count == processCount)
       if(!completed){
          printf("waiting for completion: %d/%d", file.count, processCount)
          Sys.sleep(3)
          }
       } # while
     printf("complete %d/%d", processCount, processCount)

    files <- list.files(path="BACH1", pattern="*.RData")
    tbls <- list()
    for(i in seq_len(length(files))){
        tbls[[i]] <- get(load(file.path("BACH1", files[i])))
        }
      # make sure there were some empty sections
    checkTrue(0 %in% unlist(lapply(tbls, nrow)))
    tbl <- do.call(rbind, tbls)

    sort.order <- order(tbl$start, decreasing=FALSE)
    tbl <- unique(tbl[sort.order,])
    checkEquals(dim(tbl), c(10, 9))

} # test_bach1.with.empty.section
#----------------------------------------------------------------------------------------------------
test_bach1.geneHancer <- function()
{
    message(sprintf("--- test_bach1.geneHancer"))

       # delete any results files left over from previous run
    files <- list.files(path="BACH1", pattern="*.RData")
    if(length(files) > 0)
       unlink(file.path("BACH1", files))

    fimoThreshold <- 1e-3
    processCount <- 30

    runner <- BigFimo$new("BACH1", fimoThreshold, processCount)
    runner$runMany()

    Sys.sleep(10)
    completed <- FALSE
    while(!completed){
       file.count <- length(list.files(path="BACH1", pattern="*.RData"))
       completed <- (file.count == processCount)
       if(!completed){
          printf("waiting for completion: %d/%d", file.count, processCount)
          Sys.sleep(3)
          }
       } # while
     printf("completed %d/%d",  processCount,processCount)

  #  files <- list.files(path="BACH1", pattern="*.RData")
  #  tbls <- list()
  #  for(i in seq_len(length(files))){
  #      tbls[[i]] <- get(load(file.path("BACH1", files[i])))
  #      }

  #  tbl <- do.call(rbind, tbls)
  #  sort.order <- order(tbl$start, decreasing=FALSE)
  #  tbl <- unique(tbl[sort.order,])
  #  checkEquals(dim(tbl), c(2021, 9))

} # test_bach1.geneHancer
#---------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()

