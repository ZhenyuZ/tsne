library(DESeq2)

load("count.rda")

ptm <- proc.time()
rld = varianceStabilizingTransformation(count)
proc.time() - ptm


2
   user  system elapsed
  7.540  27.260   6.953
4
   user  system elapsed
  7.307  28.274   6.956

8
 16.438  47.593  15.296

16
 29.228  96.303  23.343
