args <- commandArgs(trailingOnly=TRUE)
printf("arg count: %d", length(args))
stopifnot(length(args) == 5)

source("~/github/fimoService/batchMode/fimoBatchTools.R")
# bach1 promoter: chr21:29,298,789-29,299,657
targetGene <- args[1]
chrom <- args[2]
start <- as.numeric(args[3])
end <- as.numeric(args[4])
size <- 1 + end - start

fimo.threshold <- as.numeric(args[5])

tbl.oi <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
meme.file <- "~/fimoParallel/jaspar2018-hocomocoCoreA.meme"
printf("%s:%d-%d  (%d)", chrom, start, end, size);
printf("fimo threshold: %20.10f", fimo.threshold)

tbl.fimo <- fimoBatch(tbl.oi, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
printf("fimo hits: %d", nrow(tbl.fimo))
filename <- sprintf("fimo.%s.%s.%d.%d.RData", targetGene, chrom, start, end)
path <- file.path(targetGene, filename)
save(tbl.fimo, file=path)
