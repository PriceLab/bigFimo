source("~/github/bigFimo/R/FimoParallel.R")

fimoThreshold <- 1e-6
processCount <- 3

runner <- FimoParallel$new("BACH1", fimoThreshold, processCount, "chr21", 29305400, 29320060)
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
checkEquals(dim(tbl), c(161, 9))
