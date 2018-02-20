# get.hg38.segmean.r
# This is a modifed version GetGeneLevelCNA.r from
# https://github.com/ZhenyuZ/eqtl/blob/master/cnv/GetGeneLevelCNA.r
#
# Input GDC segmentation data, window size and calculate segment mean
# in genomic regions every window.size bases 

# split chromosome by window.size to create region
# input "chromosome" contains two columns in order: chr and len
# window.size is any interger > 1
make.region = function(chromosomes, window.size) {
  chromosome$nseg = ceiling(chromosome$len / window.size)
  region = data.frame(chr = unlist(apply(chromosome, 1, function(x) rep(x[1], x[3]))), 
            start = unlist(apply(chromosome, 1, function(x) (seq(as.numeric(x[3]))-1)*window.size + 1)), 
            end = unlist(apply(chromosome, 1, function(x) c((seq(as.numeric(x[3]))-1)*window.size, as.numeric(x[2]))[-1]))
    )
  region$name = with(region, paste0(chr, ":", start, "-", end))
  return(region)
}


# read segmentation file and convert to GRanges object (only autosomes)
seg2gr = function(file){
  seg = read.table(file, h=T, stringsAsFactors=F)
  # load into granges object
  gr = with(seg[seg$Chromosome %in% 1:22, ], GRanges( seqnames = paste0("chr", Chromosome),
                                                    ranges = IRanges(start = Start, end = End),
                                                    Segment_Mean = Segment_Mean))
  return(gr)
}

# calculate mean of segment_mean in gr by region 
# region is a data frame containing the following 4 columns in order: chr, start, end, name
get.segmean = function(gr, region) {
  window.mean = NULL
  for (chr in chromosome$chr) {
    my.gr = gr[seqnames(gr) == chr]
    my.region = region[region$chr == chr, ]
    window.mean = c(window.mean, apply(my.region, 1, 
                  function(x){y = restrict(my.gr, start  = as.numeric(x[2]), end = as.numeric(x[3])); 
                  weighted.mean(y$Segment_Mean, width(y))})  )
  }
  names(window.mean) = region$name
  return(window.mean)
}

# calculate extremum of segment_mean in gr by region 
# region is a data frame containing the following 4 columns in order: chr, start, end, name
get.extremum = function(gr, region) {
  window.extremum = NULL
  for (chr in chromosome$chr) {
    my.gr = gr[seqnames(gr) == chr]
    my.region = region[region$chr == chr, ]
    window.extremum = c(window.extremum, unlist(apply(my.region, 1, 
                  function(x){y = restrict(my.gr, start  = as.numeric(x[2]), end = as.numeric(x[3]))$Segment_Mean; 
                  if(length(y) > 0) {y[which.max(abs(y))]} else { NA }})))
  }
  names(window.extremum) = region$name
  return(window.extremum)
}



library(data.table)
library(GenomicRanges)

# load library and init
options(stringsAsFactors=F)
options(scipen=999)

# read file list
nocnv.files = fread("~/SCRATCH/cbs/data/all/nocnv", h=F)$V1
n = length(nocnv.files )

# collect data into list
segs = array(list(), n)
for (i in 1:n) {
  if(i %% 100 == 0) {
    cat(i, "out of", n, "done\n")
  }
  f = nocnv.files[i]
  segs[[i]] = seg2gr(paste0("~/SCRATCH/cbs/data/all/", f))
}

# annotate aliquot names
meta = fread("../meta/snp6.sample.txt", h=T)
meta$core.name = sapply(meta$filename, function(x) strsplit(x, "\\.")[[1]][1])
m = match(sapply(nocnv.files, function(x) strsplit(x, "\\.")[[1]][1]), meta$core.name)
meta = meta[m, ]
names(segs) = meta$aliquot_barcode

save(segs, file = "TCGA.nocnv.segment.rda")


# reconstruct autosomal chromosome information for GRCh38
chromosome = data.frame(chr=paste0("chr", 1:22), len=c(248956422, 242193529, 198295559, 190214555, 181538259, 
              170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 
              101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468))

# subset by window.size and get extreme cnvs
window.size = 5e6
region = make.region(chromosome, window.size)

extremum = matrix(0, nrow(region), n)
for (i in 1:n) {
  if(i %% 100 == 0) {
    cat(i, "out of", n, "done\n")
  }
  gr = segs[[i]]
  extremum[, i] = get.extremum(gr, region)
}



library(data.table)

cnv = fread("TCGA.cnv.window_1M.txt", skip=1)
region = cnv$V1
cnv = as.matrix(cnv[, 2:ncol(cnv), with=F])
header = readLines("TCGA.cnv.window_1M.txt", n=1)
aliquot = strsplit(header, "\t")[[1]]
colnames(cnv) = aliquot
rownames(cnv) = region
save(cnv, file="TCGA.cnv.1M.matrix.rda")


# aliquot selection
meta = fread("../meta/snp6.sample.txt", h=T)
aliquot2 = SelectAliquot(aliquot, my.type = c("01", "03"))
cnv = cnv[ , match(aliquot2, aliquot)]
meta = meta[match(aliquot2, meta$aliquot_barcode), ]


# cut off
cnv2 = cnv
cutoff = 0.1
cnv2[which(is.na(cnv2))] = 0
cnv2[which(cnv2 < -cutoff)] = -1
cnv2[which(cnv2 > cutoff)] = 1
cnv2[which(cnv2 < cutoff & cnv > -cutoff)] = 0

cnv2 = CompletenessFilter(cnv2, feature.missing = 1, sample.mising = 1, dup.removal = T)
w = duplicated(cnv2)
cnv2 = cnv2[-w, ]

library(Rtsne)



pca2 = pca(cnv2)
x = pca2$x
w = duplicated(x)
tsne2 = tsne(x, permute=1, use=200)

