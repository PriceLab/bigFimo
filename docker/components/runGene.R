    targetGene <- "ZBTB7A"

    if(file.exists(targetGene)){
       unlink(targetGene, recursive=TRUE)
       }

       # hg38, targetGene's chomLoc plus some upstream
    chrom <- "chr19"
    start <- 4040801
    end <- 4077194

      #-----------------------------------------------------------------
      # partition the target region into 3 slightly overlapping parts
      # overlap here is 60, with duplicates removed in later processing
      #-----------------------------------------------------------------

    processCount <- 3    # number of processes, hence also number of genomic regions
    half.overlap <- 30
    landmarks <- seq(from=start, to=end, length.out=processCount+1)
    starts <- landmarks[1:processCount]
    ends <- landmarks[2:(processCount+1)]
    starts <- starts - half.overlap
    ends <- ends + half.overlap
    tbl.roi <- data.frame(chrom=chrom, start=starts, end=ends, stringsAsFactors=FALSE)

    motifs <- query(MotifDb, c("sapiens"), c("jaspar2022"))
    meme.file <- file.path(targetGene, "motifs.meme")

    bf <-  BigFimo$new(targetGene,
                       tbl.oc=tbl.roi,  # not actually open chromatin - this includes all bases
                       processCount=processCount,
                       fimoThreshold=1e-6,
                       use.genehancer=FALSE,
                       gh.elite.only=FALSE,
                       maxGap.between.oc.and.gh=NA,
                       chrom=chrom, start=start, end=end,
                       motifs=motifs,
                       meme.file.path=meme.file)

   filenames.roi <- bf$createFimoTables()
   bf$createMemeFile()
   checkTrue(file.exists(meme.file))

   checkEquals(length(filenames.roi), processCount)
      # make sure these RData files, which will be read by a script that directly
      # runs FIMO, match the regions calculated above
   for(i in seq_len(length(filenames.roi))){
       tbl.r <- get(load(file.path(targetGene, filenames.roi[i])))
       checkEquals(tbl.r, tbl.roi[i,])
       }
   bf$runMany()
   bf$waitForCompletion(sleepInterval=1)

   fimo.output.files.by.region <-
              list.files(path=targetGene,
                         pattern=sprintf("^fimo.%s.*RData", targetGene))

   checkEquals(length(fimo.output.files.by.region), processCount)
   tbl.fimo <- bf$combineResults()
