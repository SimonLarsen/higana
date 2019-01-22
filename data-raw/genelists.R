library(usethis)

files <- c(
  hg18 = "glist-hg18",
  hg19 = "glist-hg19",
  hg38 = "glist-hg38"
)

genes <- lapply(files, function(file) {
  read.table(file, header=FALSE, col.names=c("chr","start","end","name"))
})

usethis::use_data(genes, overwrite=TRUE)
