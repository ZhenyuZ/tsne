library(data.table)

clin = read.table("tcga.clin.patient.coded.pruned.txt", h=T, sep="\t", stringsAsFactors=F)
class = sapply(clin, class)
load("TCGA.autosome.pc10.tsne2.perm1000.met_cnv_rna_mirna.rda")


projects = unique(clin$project)
variables = colnames(clin)[-1]
n_project = length(projects)
n_variable = length(variables)
pcnames = paste0("PC", 1:10)

# Merge PC and Clinical dataset
pc = met.pc
pc = data.frame(pc) %>% 
	mutate(aliquot = rownames(pc), 
			bcr_patient_barcode = substr(aliquot, 1, 12), 
			sample_code = substr(aliquot, 14 ,15)) 
data = data.table(merge(clin, pc))
data = data[(project == "LAML" & sample_code == "03") | 
			(project == "SKCM" & sample_code == "06") |
			(!project %in% c("LAML", "SKCM") &  sample_code == "01")]
data = data[!duplicated(data$bcr_patient_barcode)]

d = array(data.frame(), n_project)
names(d) = projects 

for(my.project in projects) {
	print(my.project)
	my.data = data[project == my.project]
	result = data.frame(matrix(NA, 10, n_variable))
	colnames(result) = variables
	rownames(result) = pcnames
	for(x in variables) {
		for(y in 1:10) {
			result[y, x] = GetCorP(my.data[[paste0("PC", y)]], my.data[[x]])
		}
	}
	d[[my.project]] = result
}

# repeat the code above for all data types
rna = d
met = d
mirna = d
cnv = d

# combne
d = array(data.frame(), n_project)
names(d) = projects 
for(my.project in projects) {
	result = rbind(rna[[my.project]], mirna[[my.project]], met[[my.project]], cnv[[my.project]])
	rownames(result) = c(paste0("rna.", pcnames), 
						paste0("mirna.", pcnames), 
						paste0("met.", pcnames), 
						paste0("cnv.", pcnames))
	d[[my.project]] = result
}
saveRDS(d, file="tcga.rna_mirna_met_cnv_pca.clin.correlation.rds")

# total non-NA p-values is 47080, which translate to 1e-6 cutoff
sum(unlist(lapply(d, function(x) sum(sum(!is.na(x))))))

# collect metrics with significant p-values
p = NULL
for(my.project in projects) {
	data = as.matrix(d[[my.project]])
	w1 = which(data < 1e-6)
	if(length(w1) > 0) {
		w2 = data.frame(which(data < 1e-6, arr.ind=T))
		w2$pval = data[w1]
		w2$project = my.project
		p = rbind(p, w2)
	}
}
p$pc = rownames(data)[p$row]
p$variable = colnames(data)[p$col]
p$type = sapply(p$pc, function(x) strsplit(x, "\\.")[[1]][1])
p$pc = sapply(p$pc, function(x) strsplit(x, "\\.")[[1]][2])
p = data.table(p[, c("project", "type", "pc", "variable", "pval")])
p = p[order(p$pval)]
p$pc = gsub("PC", "", p$pc)

write.table(p, "tcga.rna_mirna_met_cnv_pca.clin.significant.correlation.txt", col.names=T, row.names=F, sep="\t", quote=F)

GetCorP = function(x, y) {
	if(class(x) %in% c("integer", "logical", "numeric") & class(y) %in% c("integer", "logical", "numeric") & length(unique(y[!is.na(y)])) > 1) {
		l = lm(x ~ y)
		f = summary(l)$fstatistic
		p = pf(f[1],f[2],f[3],lower.tail=F)
	} else if(class(x) %in% c("integer", "logical", "numeric") & class(y) %in% c("character") & length(unique(y[!is.na(y)])) > 1) {
		fit = aov(x ~ y)
		p = summary(fit)[[1]][["Pr(>F)"]][1]
	} else {
		p = NA
	}
	return(p)
}







clin = read.table("tcga.clin.patient.coded.pruned.txt", h=T, sep="\t", stringsAsFactors=F)
class = sapply(clin, class)
w = which(class %in% c("integer", "numeric"))

for(i in w) {
	y =  clin[, i]
	if (length(unique(y[!is.na(y)])) == 1) cat(colnames(clin)[i], " : 1 enum\n")
	if (length(unique(y[!is.na(y)])) == 2) {
		cat(colnames(clin)[i], " : 2 enum\n")
		print(table(y))
		cat("\n")
	}
}









