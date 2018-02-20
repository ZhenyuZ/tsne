library(data.table)
source("./tsne.util.r")

meta = fread("met.meta.txt")
load("met.rda")

# find patient with at least one non primary tumor sample
w = which(meta$sample_code %in% c("02", "05", "06","07","11"))
temp = table(meta$patient[which(meta$patient %in% meta$patient[w])])
patient = names(temp[temp>1])
w = which(meta$patient %in% patient)
met2 = met[, w]

# ptm <- proc.time()
met2 = CompletenessFilter(met2)
# proc.time() - ptm
#     user   system  elapsed
# 1492.042  118.861 1611.207

# NA imputation
ptm <- proc.time()
met2 = ImputeNA(met2, large.matrix=T)
proc.time() - ptm

save(met2, file="met2.nonsolidtumorpatient.rda")

# object too large for PCA
# filter by SD
std = apply(met2, 1, sd)
met2 = met2[which(std > 0.2), ]


# pca
ptm <- proc.time()
pca = pca(met2)
proc.time() - ptm
save(pca, file="pca.nonsolidtumorpatient.rda")

pc = data.frame(pca$x[, 1:2])
pc$patient = substr(rownames(pc), 1, 12)
pc$sample_type = substr(rownames(pc), 14, 15)
pc$project = meta$project[match(pc$patient, meta$patient)]
save(pc, file="pc.nonsolidtumorpatient.rda")

png("met2.nonsolidtumorpatient.tsne.perm100.png", width=1200, height=800)
ggplot(pc[which(pc$project %in% c("LUAD", "LUSC"))], aes(x=PC1, y=PC2, col=project, shape = sample_type, group=patient)) + 
	geom_point(size = 3) + geom_line(lty=2, alpha=0.5)

# t-SNE
library(Rtsne)
ptm <- proc.time()
tsne = tsne(pca$x, permute=100, use=50)
proc.time() - ptm
save(tsne, file="tsne.nonsolidtumorpatient.rda")


tsne$patient = substr(rownames(tsne), 1, 12)
tsne$sample_type = substr(rownames(tsne), 14, 15)
tsne$project = meta$project[match(tsne$patient, meta$patient)]
save(tsne, file="tsne.nonsolidtumorpatient.rda")

# t-SNE of randome sampleing
# pca3 = pca(met3)
# tsne3 = tsne(pca3$x)
# tsne3.pdata = annotate.tsne(tsne3, aliquot)
# tsne2.pdata = annotate.tsne(tsne2, aliquot)
# p3 = ggplot(tsne3.pdata, aes(x=feature_1, y=feature_2, col=project)) + geom_point()

pdata$sh = rep(c(0, 2, 4), 11)

png("met2.nonsolidtumorpatient.tsne.perm100.png", width=1200, height=800)
ggplot(tsne, aes(x=feature_1, y=feature_2, col=project, shape = sample_type, group=patient)) + 
	geom_point(size = 3) + geom_line(lty=2, alpha=0.5)
p
dev.off()

