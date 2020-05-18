library(data.table)
library(ggplot2)

setwd("~/SCRATCH/tsne")
# source("./tsne.util.r")

load("./rna/TCGA.rna.filter1.allgene.pca.rda")
pca.rna <- pca
rm(list=c("pca"))
patient.rna <- substr(rownames(pca.rna$x), 1, 12)

load("./mirna/TCGA.mirna.filter1.allgene.pca.rda")
pca.mirna <- pca2
rm(list=c("pca2"))
patient.mirna <- substr(rownames(pca.mirna$x), 1, 12)

load("./meth/TCGA.met.pca.perm100.rda")
pca.met <- pca
rm(list=c("pca"))
patient.met <- substr(rownames(pca.met$x), 1, 12)

load("./cnv/TCGA.cnv.1m.pca.rda")
pca.cnv <- pca2
rm(list=c("pca2"))
patient.cnv <- substr(rownames(pca.cnv$x), 1, 12)

# find common patient
patients <- Reduce(intersect, list(patient.met, patient.rna, patient.mirna, patient.cnv))

rna <- pca.rna$x[match(patients, patient.rna), ]
mirna <- pca.mirna$x[match(patients, patient.mirna), ]
met <- pca.met$x[match(patients, patient.met), ]
cnv <- pca.cnv$x[match(patients, patient.cnv), ]


# merge all PCs 
data <- cbind(rna, 
							mirna[match(patients, patient.mirna), ], 
							met[match(patients, patient.met), ], 
							cnv[match(patients, patient.cnv), ])
rownames(data) <- patients



# in order to conbine PCs from different dataset, we need to scale the PCs
# we will calculate sqrt(sum of sdev^2), and normalize PCs agaist this factor 
sum.var.rna <- sum(apply(rna, 2, sd)^2)
rna <- rna / sqrt(sum.var.rna)
colnames(rna) <- paste0("rna_", colnames(rna))


sum.var.mirna <- sum(apply(mirna, 2, sd)^2)
mirna <- mirna / sqrt(sum.var.mirna)
colnames(mirna) <- paste0("mirna_", colnames(mirna))

sum.var.met <- sum(apply(met, 2, sd)^2)
met <- met / sqrt(sum.var.met)
colnames(met) <- paste0("met_", colnames(met))

sum.var.cnv <- sum(apply(cnv, 2, sd)^2)
cnv <- cnv / sqrt(sum.var.cnv)
colnames(cnv) <- paste0("cnv_", colnames(cnv))

# merge all PCs 
data <- cbind(rna, mirna, met, cnv)
rownames(data) <- patients

# UMAP and save results
d <- data.table(patient = patients)

set.seed(1)
umap.rna <- umap(data[, which(grepl("^rna_PC", colnames(data)))])
d$umap.rna.x <- umap.rna$layout[, 1]
d$umap.rna.y <- umap.rna$layout[, 2]
save(umap.rna, file="umap.rna.rda")

set.seed(1)
umap.mirna <- umap(data[, which(grepl("^mirna_PC", colnames(data)))])
d$umap.mirna.x <- umap.mirna$layout[, 1]
d$umap.mirna.y <- umap.mirna$layout[, 2]
save(umap.mirna, file="umap.mirna.rda")

set.seed(1)
umap.met <- umap(data[, which(grepl("^met_PC", colnames(data)))])
d$umap.met.x <- umap.met$layout[, 1]
d$umap.met.y <- umap.met$layout[, 2]
save(umap.met, file="umap.met.rda")

set.seed(1)
umap.cnv <- umap(data[, which(grepl("^cnv_PC", colnames(data)))])
d$umap.cnv.x <- umap.cnv$layout[, 1]
d$umap.cnv.y <- umap.cnv$layout[, 2]
save(umap.cnv, file="umap.cnv.rda")

set.seed(1)
umap.all <- umap(data)
d$umap.all.x <- umap.all$layout[, 1]
d$umap.all.y <- umap.all$layout[, 2]
save(umap.all, file="umap.all.rda")

# add metadata
meta <- fread("all.tsne.data.20170724.txt")
m <- match(d$patient, meta$patient)
d <- cbind(d, meta[m, c("project", "death", "day", "age", "gender", "race")])
write.table(d, "tcga.umap.tsv", col.names=T, row.names=F, sep="\t", quote=F)





png("umap.rna.global.png", width=1280, height=800)
p = ggplot(d, aes(x=umap.rna.x, y=umap.rna.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()

png("umap.rna.local.png", width=1280, height=800)
p = ggplot(d[umap.rna.x > -20 & umap.rna.y > -10], aes(x=umap.rna.x, y=umap.rna.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


png("umap.met.global.png", width=1280, height=800)
p = ggplot(d, aes(x=umap.met.x, y=umap.met.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


png("umap.mirna.global.png", width=1280, height=800)
p = ggplot(d, aes(x=umap.mirna.x, y=umap.mirna.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


png("umap.cnv.global.png", width=1280, height=800)
p = ggplot(d, aes(x=umap.cnv.x, y=umap.cnv.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


png("umap.all.global.png", width=1280, height=800)
p = ggplot(d, aes(x=umap.all.x, y=umap.all.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()

png("umap.all.hnsc.png", width=1280, height=800)
p = ggplot(d[umap.all.x > -3.5 & umap.all.x < 2.8 & umap.all.y > -4.3 & umap.all.y < 2.4], aes(x=umap.all.x, y=umap.all.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()

png("umap.all.hnsc.only.png", width=1280, height=800)
p = ggplot(d[project == "HNSC"], aes(x=umap.all.x, y=umap.all.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()




ggplot(d[project == "HNSC"], aes(x=umap.all.x, y=umap.all.y, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8))
