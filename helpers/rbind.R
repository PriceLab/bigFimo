args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
targetGene <- args[1]

fimo.files <- list.files(pattern="^fimo")
tbls <- list()
for(f in fimo.files){
    tbl <- get(load(f))
    tbls[[f]] <- tbl
    }
length(tbls)
lapply(head(tbls), dim)
tbl.fimo <- do.call(rbind, tbls)
rownames(tbl.fimo) <- NULL
dim(tbl.fimo)
new.order <- order(tbl.fimo$start, decreasing=FALSE)
head(new.order)
tbl.fimo <- tbl.fimo[new.order,]
head(tbl.fimo)
filename <- sprintf("tbl.fimo.%s.RData", targetGene)
save(tbl.fimo, file=filename)
