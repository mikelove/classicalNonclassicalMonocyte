scratch <- "fastq"

# see https://gist.github.com/mikelove/f539631f9e187a8931d34779436a1c01 for accession2url() definition

source("https://gist.githubusercontent.com/mikelove/f539631f9e187a8931d34779436a1c01/raw/6e6633aa5123358b70390ab738be1eef03a3df31/accession2url.R")

x <- scan("SRR_Acc_List.txt", what="char")

# problem w 3, 32

for (i in 1:length(x)) {
  print(paste("---",i,"---"))
  run <- x[i]
  file <- paste0(run,".fastq.gz")
  url <- file.path(accession2url(run), file)
  dest <- file.path(scratch, file)
  if (!file.exists(dest))
    download.file(url, dest)
}
