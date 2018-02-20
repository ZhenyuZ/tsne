library(data.table)

######################################################
# working on level 4 pick tier 1 clinical data first
######################################################

# get file manifest with `find . -name "*clin.merged.picked.txt" > picked.files.txt`
# get file manifest with `find . -name "All_CDEs.txt" > All_CDEs.files.txt`
f = fread("picked.files.txt", h=F, col.names="pick")
f$project = sapply(basename(f$pick), function(x) strsplit(x, "\\.")[[1]][1])
f2 = fread("All_CDEs.files.txt", h=F, col.names="cde")
f2$project = gsub(".+gdac.broadinstitute.org_([A-Z]+).Clinical_Pick_Tier1.+", "\\1", f2$cde)
f = merge(f, f2)

n_project = nrow(f)
clin = array(data.frame(), n_project)
names(clin) = f$project

for(i in 1:n_project){
#	raw.pick = read.table(f$pick[i], stringsAsFactors=F, sep="\t", h=F)
	raw.cde = read.table(f$cde[i], stringsAsFactors=F, sep="\t", h=F, quote='"')

	cde = data.frame(t(raw.cde)[2:ncol(raw.cde), ])
	colnames(cde) = raw.cde$V1
	cde$project = f$project[i]
	clin[[i]] = cde
}

cnames = unique(unlist(lapply(clin, colnames)))
merge = bind_rows(clin)
write.table(merge, "raw.tcga.clin.full.20160128.txt", col.names=T, row.names=F, sep="\t", quote=F)
saveRDS(clin, "raw.tcga.clin.full.20160128.rds")

# removing singletons 
w = NULL
for(i in 1: ncol(merge)) {
	temp = merge[, i]
	if(length(unique(temp[!is.na(temp)])) < 2) w = c(w, i)
}
merge = merge[, -w] 
write.table(merge, "raw.tcga.clin.nonsingleton.20160128.txt", col.names=T, row.names=F, sep="\t", quote=F)





######################################################
# ignore below for the full clin read in
######################################################
# read merged clin file list 
file = fread("f", h=F)$V1
project = sapply(f, function(x) strsplit(x, "\\.")[[1]][1])
clin = array(data.table(), length(project))
names(clin) = project
for(i in 1:33){ f = file[i]; clin[[i]] = fread(f, h=F)}

# merge
field = NULL
for(i in 1:33) field = c(field, clin[[i]]$V1)
field = unique(field)

data =  data.table(property = field)
for(i in 1:33) {
	f = clin[[1]]$V1
	m = match(field, f)
	data = cbind(data, clin[[i]][, -"V1", with=F][m, ])
}
saveRDS(data, file="tcga.clin.raw.wide.rds")
write.table(data, "tcga.clin.raw.wide.txt", col.names=F, row.names=F, sep="\t", quote=F)

clin = t(data)
write.table(clin, "tcga.clin.raw.txt", col.names=F, row.names=F, sep="\t", quote=F)
clin = fread("tcga.clin.raw.txt")

clin[clin==""] = NA

meta = data.table(field)
meta$class = sapply(field, function(x) {y=strsplit(x, "\\.")[[1]]; paste(y[1:(length(y)-1)], collapse=".")})
meta$value = sapply(field, function(x) {y=strsplit(x, "\\.")[[1]]; y[length(y)]})
meta$level = sapply(meta$class, function(x) length(strsplit(x, "\\.")[[1]]))
meta$l1 = sapply(meta$class, function(x) strsplit(x, "\\.")[[1]][1])
meta$l2 = sapply(meta$class, function(x) strsplit(x, "\\.")[[1]][2])

save(clin, meta, file="clin.rda")





