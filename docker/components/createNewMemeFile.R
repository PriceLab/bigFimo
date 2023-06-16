library(MotifDb)
motifs <- query(MotifDb, c("sapiens"), c("jaspar2022", "hocomocov11-core-A"))
length(motifs)
filename <- "human-jaspar2022-hocomoco-core-A.meme"
export(motifs, con=filename, format="meme")
