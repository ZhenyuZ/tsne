library(data.table)
library(ggplot2)


############################################################
# RNA
############################################################

# load rna data
load("~/SCRATCH/tsne/rna/tcga.gaf.rsem.rda")
tcga = data
load("~/SCRATCH/tsne/rna/rna.uqfpkm.rda")
gdc = uqfpkm

# load metadata
aliquot = fread("disease_aliquot.txt", h=T)
gencode = fread("~/SCRATCH/tsne/meta/gencode.v22.genes.txt", h=T)

# exact genes names 
gdc.genes = rownames(gdc)
gdc.genes = toupper(gencode$gene_name[match(gdc.genes, gencode$gene_id)])

tcga.genes = rownames(tcga)
tcga.genes = toupper(sapply(tcga.genes, function(x) strsplit(x, "\\|")[[1]][1]))

# remove duplicated gene names, and find commons ones
gdc.uniq.genes = gdc.genes[-which(duplicated(gdc.genes) | duplicated(gdc.genes, fromLast=T))]
tcga.uniq.genes = tcga.genes[-which(duplicated(tcga.genes) | duplicated(tcga.genes, fromLast=T))]
common.genes = intersect(gdc.uniq.genes, tcga.uniq.genes)

# subset by common genes
gdc = gdc[match(common.genes, gdc.genes), ]
tcga = tcga[match(common.genes, tcga.genes), ]

# find common aliquot
gdc.aliquot = colnames(gdc)
rna.meta = fread("~/SCRATCH/tsne/meta/rna.meta.txt", h=T)
m = match(gdc.aliquot, rna.meta$htseq_count_id)
gdc.aliquot = rna.meta$barcode[m]
tcga.aliquot = colnames(tcga)
common.aliquot = intersect(gdc.aliquot, tcga.aliquot)

# subset by common aliquot
gdc = gdc[, match(common.aliquot, gdc.aliquot)]
tcga = tcga[, match(common.aliquot, tcga.aliquot)]

# save(gdc, tcga, file = "~/SCRATCH/tsne/rna/common/common.tcga_rsem.gdc_uqfpkm.rda")

# calculate spearman correlation
corr = array(0, ncol(gdc))
for (i in 1: ncol(gdc)) {
	if(i %% 1000 == 0) {
		cat(i, "out of", ncol(gdc), "samples completed\n")
	}
	corr[i] = cor(gdc[,i], tcga[,i], method = "spearman")
}
corr = data.frame(corr)
corr$project = aliquot$disease[match(common.aliquot, aliquot$aliquot)]
corr$aliquot = common.aliquot

write.table(corr, "~/SCRATCH/tsne/rna/common/tcga.gdc.spearman.corr.txt", col.names=T, row.names=F, sep="\t", quote=F)


png("~/Downloads/paper/tcga_rsem.gdc_uqfpkm.spearman.png", width=1024, height=800)
theme_set(theme_gray(base_size = 18))
p = ggplot(corr, aes(x=project, y=corr, fill=project))  + 
	geom_boxplot(alpha=0.8) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	labs(x = "Project", y="Spearman Correlation between TCGA RSEM and GDC UQ-FPKM values") + 
	theme(legend.position="none")
p
dev.off()


# average correlation 0.944
# sd correlation 0.007


############################################################
# miRNA
############################################################
library(data.table)


# load mirna data
load("~/SCRATCH/tsne/mirna/tcga.mirna.count.rda")
tcga = data
load("~/SCRATCH/tsne/mirna/mirna.count.rda")
gdc = count

# load metadata
aliquot = fread("disease_aliquot.txt", h=T)

# exact genes names, and find commons ones
gdc.genes = rownames(gdc)
tcga.genes = rownames(tcga)
common.genes = intersect(gdc.genes, tcga.genes)

# subset by common genes
gdc = gdc[match(common.genes, gdc.genes), ]
tcga = tcga[match(common.genes, tcga.genes), ]

# find common aliquot
gdc.aliquot = colnames(gdc)
mirna.meta = fread("~/SCRATCH/tsne/meta/mirna.meta.txt", h=T)
m = match(gdc.aliquot, mirna.meta$aligned_file_gdcid)
gdc.aliquot = mirna.meta$barcode[m]
tcga.aliquot = colnames(tcga)
common.aliquot = intersect(gdc.aliquot, tcga.aliquot)

# subset by common aliquot
gdc = gdc[, match(common.aliquot, gdc.aliquot)]
tcga = tcga[, match(common.aliquot, tcga.aliquot)]

# save(gdc, tcga, file = "~/SCRATCH/tsne/mirna/common/common.tcga_count.gdc_count.rda")

# transform to tpm
total = apply(gdc, 2, sum) / 1e6
gdc = t(t(gdc)/total)
total = apply(tcga, 2, sum) / 1e6
tcga = t(t(tcga)/total)

# calculate spearman correlation
corr = array(0, ncol(gdc))
for (i in 1: ncol(gdc)) {
	if(i %% 1000 == 0) {
		cat(i, "out of", ncol(gdc), "samples completed\n")
	}
	corr[i] = cor(gdc[,i], tcga[,i], method = "spearman")
}
corr = data.frame(corr)
corr$project = aliquot$disease[match(common.aliquot, aliquot$aliquot)]
corr$aliquot = common.aliquot

write.table(corr, "~/SCRATCH/tsne/mirna/common/tcga.gdc.spearman.corr.txt", col.names=T, row.names=F, sep="\t", quote=F)


png("tcga_count.gdc_count.spearman.png", width=1024, height=800)
theme_set(theme_gray(base_size = 18))
p = ggplot(corr, aes(x=project, y=corr, fill=project))  + 
	geom_boxplot(alpha=0.8) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	labs(x = "Project", y="Spearman Correlation between TCGA and GDC miRNA TPM values") + 
	theme(legend.position="none")
p
dev.off()





