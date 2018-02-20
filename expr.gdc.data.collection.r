library(data.table)
####################
# miRNA data aggregation
# ls *mirnas.quantification.txt > f
# while read -r f; do tail -n+2 $f | awk '{printf $3"\t"}; END{printf "\n"}' >> combine.tpm.txt; done <"f"
# while read -r f; do tail -n+2 $f | awk '{printf $2"\t"}; END{printf "\n"}' >> combine.count.txt; done <"f"
####################
setwd("~/SCRATCH/tsne/mirna/data")

count = fread("combine.count.txt", h=F)
tpm = fread("combine.tpm.txt", h=F)
files = fread("f", h=F)$V1

sample = fread(files[1], h=T)
mirna_id =sample$miRNA_ID
sample_id = sapply(files, function(x) strsplit(x, "\\.")[[1]][1])

count = t(as.matrix(count))
count = count[-nrow(count), ]
rownames(count) = mirna_id
colnames(count) = sample_id
write.table(count, "mirna.count.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(count, file="mirna.count.rda")

tpm = t(as.matrix(tpm))
tpm = tpm[-nrow(tpm), ]
rownames(tpm) = mirna_id
colnames(tpm) = sample_id
write.table(tpm, "mirna.tpm.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(tpm, file="mirna.tpm.rda")

######################
# autosomal transcript 
######################
# autosomal transcript 
transcript = fread("/mnt/SCRATCH/tsne/meta/mirbase.v21.transcript.txt")
autosomal.transcript = transcript[seqname %in% paste0("chr",1:22)]

m = match(autosomal.transcript$Name, rownames(count))
autosomal.count = count[m, ]
mirna.meta = fread("~/SCRATCH/tsne/meta/mirna.meta.txt", h=T)
m = match(colnames(count), mirna.meta$aligned_file_gdcid)
colnames(autosomal.count) = mirna.meta$barcode[m]

write.table(autosomal.count, "autosomal.mirna.count.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(autosomal.count, file="autosomal.mirna.count.rda")

m = match(autosomal.transcript$Name, rownames(tpm))
autosomal.tpm = tpm[m, ]
m = match(colnames(tpm), mirna.meta$aligned_file_gdcid)
colnames(autosomal.tpm) = mirna.meta$barcode[m]

write.table(autosomal.tpm, "autosomal.mirna.tpm.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(autosomal.tpm, file="autosomal.mirna.tpm.rda")



########################################################################################################################
# RNA UQFPKM data aggregation
# ls *PKM-UQ.txt.gz > f.uqfpkm
# while read -r f; do tail -n+2 $f | awk '{printf $3"\t"}; END{printf "\n"}' >> combine.uqfpkm.txt; done <"f.uqfpkm"
########################################################################################################################
library(data.table)
# miRNA first
setwd("~/SCRATCH/tsne/rna/data")

uqfpkm = fread("combine.uqfpkm.txt", h=F)
uqfpkm.files = fread("f.uqfpkm", h=F)$V1

sample = fread(paste("zcat", uqfpkm.files[1]), h=F)
rna_id =sample$V1
sample_id = sapply(uqfpkm.files, function(x) strsplit(x, "\\.")[[1]][1])

uqfpkm = t(as.matrix(uqfpkm))
uqfpkm = uqfpkm[-nrow(uqfpkm), ]
rownames(uqfpkm) = rna_id
colnames(uqfpkm) = sample_id
write.table(uqfpkm, "rna.uqfpkm.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(uqfpkm, file="rna.uqfpkm.rda")

######################
# Protein-coding genes 
######################
# protein coding genes
gene = fread("~/SCRATCH/tsne/meta/gencode.v22.genes.txt", h=T)
pc.gene = gene[gene_type == "protein_coding" & seqname %in% paste0("chr", 1:22)]
# remove genes with duplicated names
w = unique(c(which(duplicated(pc.gene$gene_name)), which(duplicated(pc.gene$gene_name, fromLast=T))))
pc.gene = pc.gene[-w, ]
m = match(pc.gene$gene_id, rownames(uqfpkm))
pc.uqfpkm = uqfpkm[m, ]
row.names(pc.uqfpkm) = pc.gene$gene_name

rna.meta = fread("~/SCRATCH/tsne/meta/rna.meta.txt", h=T)
m = match(colnames(pc.uqfpkm), rna.meta$htseq_count_id)
colnames(pc.uqfpkm) = rna.meta$barcode[m]

write.table(pc.uqfpkm, "pc.rna.uqfpkm.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(pc.uqfpkm, file="pc.rna.uqfpkm.rda")



########################################################################################################################
# RNA Count data aggregation
# ls *htseq.counts.gz > f.count
# while read -r f; do zcat $f | awk '{printf $2"\t"}; END{printf "\n"}' >> combine.count.txt; done <"f.count"
########################################################################################################################
library(data.table)
# miRNA first
setwd("~/SCRATCH/tsne/rna/data")

count = fread("combine.count.txt", h=F)
count.files = fread("f.count", h=F)$V1

sample = fread(paste("zcat", count.files[1]), h=F)
rna_id =sample$V1
sample_id = sapply(count.files, function(x) strsplit(x, "\\.")[[1]][1])

count = t(as.matrix(count))
count = count[-nrow(count), ]
rownames(count) = rna_id
colnames(count) = sample_id
write.table(count, "rna.count.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(count, file="rna.count.rda")

######################
# Protein-coding genes 
######################
# protein coding genes
gene = fread("~/SCRATCH/tsne/meta/gencode.v22.genes.txt", h=T)
pc.gene = gene[gene_type == "protein_coding" & seqname %in% paste0("chr", 1:22)]
# remove genes with duplicated names
w = unique(c(which(duplicated(pc.gene$gene_name)), which(duplicated(pc.gene$gene_name, fromLast=T))))
pc.gene = pc.gene[-w, ]
m = match(pc.gene$gene_id, rownames(count))
pc.count = count[m, ]
row.names(pc.count) = pc.gene$gene_name

rna.meta = fread("~/SCRATCH/tsne/meta/rna.meta.txt", h=T)
m = match(colnames(pc.count), rna.meta$htseq_count_id)
colnames(pc.count) = rna.meta$barcode[m]

write.table(pc.count, "pc.rna.count.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(pc.count, file="pc.rna.count.rda")
