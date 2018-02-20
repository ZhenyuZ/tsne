# Download data from firehose
# ls | grep gz > f
# while read -r line ; do  tar xvzf $line --strip=1; done <"f"
# ls | grep data.txt > f	
# while read -r line; do sed '2d' $line.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt | awk -F"\t" '{for(i=2;i<=NF;i=i+4){printf $i"\t"}; print ""}' > $line.beta.txt & done <"../../project"


library(data.table)
source("~/tsne.util.r")

disease = fread("../../project", h=F)$V1
aliquot = NULL
all = NULL

for(d in disease) {
	print(d)
	file = paste0(d, ".beta.txt")
	data = fread(file, h=T)
	barcode = names(data)
	if (d == "LAML"){
		barcode2keep = SelectAliquot(barcode, my.type="03")
	} else {
		barcode2keep = SelectAliquot(barcode, my.type="01")
	}
	met = as.matrix(data[, barcode2keep, with=F])
	all = cbind(all, met)
	my.aliquot = data.frame(disease = d, aliquot = barcode2keep)
	aliquot = rbind(aliquot, my.aliquot)
}

# add probesetnames
probes = unlist(fread("UCS.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", h=F, select=1, skip=2L))

met = all
rownames(met) = probes

save(met, file="TCGA.methylation.beta.before.filter.rda")

# data cleaning
# ptm <- proc.time()
met = CompletenessFilter(met)
# proc.time() - ptm
#     user   system  elapsed
# 1492.042  118.861 1611.207

# NA imputation
ptm <- proc.time()
met = ImputeNA(met, large.matrix=T)
proc.time() - ptm

save(met, file="TCGA.methylation.beta.after.filter.rda")

# object too large for PCA
# filter by SD
std = apply(met, 1, sd)
met2 = met[which(std > 0.2), ]
save(met2, file="TCGA.methylation.beta.sd0.2.rda")

# randome sampling
met3 = met[sample.int(nrow(met), nrow(met2)), ]
save(met3, file="TCGA.methylation.beta.sample.rda")
