library(data.table)
source("tsne.util.r")

load("TCGA.methylation.beta450.before.filter.rda")
met450 = met

load("TCGA.methylation.beta27.before.filter.rda")
# met27

# common probeset
probe450 = rownames(met450)
probe27 = rownames(met27)
common.probe = intersect(probe450, probe27)
met450 = met450[match(common.probe, probe450), ]
met27 = met27[match(common.probe, probe27), ]
met450.backup = met450
met27.backup = met27

# common sample
aliquot450 = colnames(met450)
aliquot27 = colnames(met27)
common.aliquot = intersect(aliquot450, aliquot27)
met450 = met450[, match(common.aliquot, aliquot450)]
met27 = met27[, match(common.aliquot, aliquot27)]
met450.backup2 = met450
met27.backup2 = met27

# combine
colnames(met450) = paste0(common.aliquot, "_450")
colnames(met27) = paste0(common.aliquot, "_27")
met = cbind(met450, met27)
met = CompletenessFilter(met)
met = ImputeNA(met)

# pca
pca = pca(met)
pc = data.frame(pca$x[, 1:3])
pc$exp = rownames(pc)
rownames(pc) = NULL
pc$aliquot = sapply(pc$exp, function(x) strsplit(x, "_")[[1]][1])
pc$platform = sapply(pc$exp, function(x) strsplit(x, "_")[[1]][2])

save(pc, file="met450.met27.before.correction.pc.rda")

# confirm the possibility to use linear transformation to convert data
png("~/Downloads/paper/met450.met27.before.correction.png", width=640, height=640)
p = ggplot(data=pc, aes(PC1, PC2, col=platform)) + 
	geom_point(size=4) + 
	theme(legend.position = c(0.9,0.9), axis.text=element_blank(), axis.ticks=element_blank()) + 
	geom_path(data=pc, aes(PC1, PC2, group=aliquot), alpha=0.5, show.legend = FALSE) 
p
dev.off()

# Trying to correct met27 data to met450
met450 = CompletenessFilter(met450.backup2)
met27 = CompletenessFilter(met27.backup2)
# common probeset
probe450 = rownames(met450)
probe27 = rownames(met27)
common.probe = intersect(probe450, probe27)
met450 = met450[match(common.probe, probe450), ]
met27 = met27[match(common.probe, probe27), ]
# common sample
aliquot450 = colnames(met450)
aliquot27 = colnames(met27)
common.aliquot = intersect(aliquot450, aliquot27)
met450 = met450[, match(common.aliquot, aliquot450)]
met27 = met27[, match(common.aliquot, aliquot27)]

met = cbind(met450, met27)
lms = apply(met, 1, function(x) lm(x[1:193] ~ x[194:386]))

# probes with 0 intersect
intercept.pvals = unlist(lapply(lms, function(x) coef(summary(x))[1, 4]))
w.intercept = which(intercept.pvals > 0.05/length(intercept.pvals))

# probes with 0 slope
slope.pvals = unlist(lapply(lms, function(x) coef(summary(x))[2, 4]))
w.slope= which(slope.pvals > 0.05/length(slope.pvals))

# common stable probes
best.probes = intersect(names(w.intercept), names(w.slope))
write.table(best.probes, "best.probes.txt", col.names=F, row.names=F, sep="\t", quote=F)

met450 = met450.backup[match(best.probes, rownames(met450.backup)), ]
met27 = met27.backup[match(best.probes, rownames(met27.backup)), ]


met450.backup3 = met450
met27.backup3 = met27

# common sample
aliquot450 = colnames(met450)
aliquot27 = colnames(met27)
common.aliquot = intersect(aliquot450, aliquot27)
met450 = met450[, match(common.aliquot, aliquot450)]
met27 = met27[, match(common.aliquot, aliquot27)]
colnames(met450) = paste0(common.aliquot, "_450")
colnames(met27) = paste0(common.aliquot, "_27")
met = cbind(met450, met27)
met = CompletenessFilter(met)
met = ImputeNA(met)
pca = pca(met)
pc = data.frame(pca$x[, 1:3])
pc$exp = rownames(pc)
rownames(pc) = NULL
pc$aliquot = sapply(pc$exp, function(x) strsplit(x, "_")[[1]][1])
pc$platform = sapply(pc$exp, function(x) strsplit(x, "_")[[1]][2])
# confirm the possibility to use linear transformation to convert data
png("met450.met27.after.correction.png", width=640, height=640)
p = ggplot(data=pc, aes(PC1, PC2, col=platform)) + 
	geom_point(size=4) + 
	theme(legend.position = c(0.9,0.9), axis.text=element_blank(), axis.ticks=element_blank()) + 
	geom_path(data=pc, aes(PC1, PC2, group=aliquot), alpha=0.5, show.legend = FALSE) 
p
dev.off()

# combine met27 and met450 data, and filter for best aliquot
met = cbind(met450.backup3, met27.backup3[, which(!colnames(met27.backup3) %in% common.aliquot)])
met = CompletenessFilter(met)
met = ImputeNA(met)
pca = pca(met)
tsne = tsne(pca$x, permute=10, use=50)
aliquot = fread("aliquot.txt")
pdata = annotate.tsne(tsne, aliquot)





