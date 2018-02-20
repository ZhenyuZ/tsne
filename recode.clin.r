library(data.table)


clin = read.table("raw.tcga.clin.nonsingleton.20160128.txt", h=T, stringsAsFactors=T, sep="\t", quote='"')
clin$bcr_patient_barcode = toupper(clin$bcr_patient_barcode)
class = sapply(clin, class)


action = array("", ncol(clin))
w = which(grepl("_text", colnames(clin)))
action[w] = "delete"
w = which(apply(clin, 2, function(x) sum(!is.na(x))) < 50)
action[w] = "delete"

for(i in 211: ncol(clin)){
	if(action[i] != "delete"){
		default = "keep"
		cat("\n")
		cat(i, ":", names(clin)[i], ":\n")
		cat("  - class:", class[i], "\n")
		if( class[i] %in% c("integer", "numeric")) {
			print(summary(clin[, i]))
		} else if( class[i] == "factor") {
			t = table(clin[, i])
			cat(" - Number of enums:", length(t), "\n")
			if(length(t) > 20) {
				cat("more than 20 enums. \n")
				cat(" - Sample enums:\n")
				print(names(t[order(-t)])[1:5])
			}
			else {
				print(names(t[order(-t)]))
			}
		} else {
			cat("class type is", class[i], "\n")
		}	
		line = readline(prompt="Press [enter] to keep, [i] to invest, [d] to delete, [r] to recode: \n")
		if(line == "d") default = "delete"
		if(line == "r") default = "recode"
		if(line == "i") default = "invest"
		action[i] = default
	}
}




map = data.frame(variable=NULL, value=NULL, alias=NULL)
recode = which(action=="recode")
for(i in 1: length(recode)){
	w = recode[i]
	temp = names(table(clin[, w]))
	map = rbind(map, data.frame(variable = colnames(clin)[w], 
								value = names(table(clin[, w])), 
								alias = names(table(clin[, w]))))
}
write.table(map, "map.txt", col.names=T, row.names=F, sep="\t", quote=F)

map = data.frame(variable=NULL, value=NULL, alias=NULL)
recode = which(action=="recode" & class=="character")
for(i in 1: length(recode)){
	w = recode[i]
	temp = names(table(clin[, w]))
	map = rbind(map, data.frame(variable = colnames(clin)[w], 
								value = names(table(clin[, w])), 
								alias = ""))
}
write.table(map, "map.txt", col.names=T, row.names=F, sep="\t", quote=F)


# manual editing
map = fread("map.txt")
recode = data.frame(lapply(clin[, which(action=="recode")], as.character), stringsAsFactors=FALSE)
for(i in 1: nrow(map)) {
	data = recode[, map$variable[i]]
	data[which(data==map$value[i])] = map$alias[i]
	recode[, map$variable[i]] = data
}
colnames(recode) = paste0("recoded_", colnames(recode))
action[which(action=="recode")] = "to_be_recoded"
action = c(action, rep("recoded", ncol(recode)))
clin = cbind(clin, recode)

# remove variable with number non-NA values < 50 or second most enums < 20 
w = which(action %in% c("keep", "recoded"))
to_remove = NULL
for(i in w) {
	data = table(clin[, i])
	data = data[order(-data)]
	if(sum(data) < 50 | data[2] < 20) to_remove = c(to_remove, i)
}
to_remove = to_remove[-1]
action[to_remove] = "delete"

save(clin, action, file="tcga.clin.patient.coded.all.rda")

clin = clin[, which(action %in% c("keep", "recoded"))]
write.table(clin, "tcga.clin.patient.coded.pruned.txt", col.names=T, row.names=F, sep="\t", quote=F)
clin = read.table("tcga.clin.patient.coded.pruned.txt", h=T, stringsAsFactors=T, sep="\t", quote='"')
saveRDS(clin, "tcga.clin.patient.coded.pruned.rds")















# manual editing
map = fread("map.txt")
for(i in 1: nrow(map)) {
	data = clin[, map$variable[i]]
	data[which(data==map$value[i])] = map$alias[i]
	clin[, map$variable[i]] = data
}
write.table(clin, "tcga.clin.patient.coded.all.txt", col.names=T, row.names=F, sep="\t", quote=F)
clin = clin[, which(action != "ignore")]
write.table(clin, "tcga.clin.patient.coded.pruned.txt", col.names=T, row.names=F, sep="\t", quote=F)
saveRDS(clin, "tcga.clin.patient.coded.pruned.rds")

clin2 = fread("tcga.clin.patient.coded.pruned.txt")
temp = NULL
class = sapply(clin, class)
w = which(class %in% c("logical", "character"))
for(i in 1: length(w)) {
	data = unlist(clin[, w[i], with=F])
	if(table(data)[order(-table(data))][2] < 20) {
		temp = c(temp, w[i])
	}
}




for(i in 1:ncol(clin)){
	data = clin[[i]]
	temp = unique(data)
	temp = temp[!is.na(temp)]
	if(length(temp)==4 & sum(temp %in% c("0", "1", "FALSE", "TRUE")) ==4) {
		data[which(data=="TRUE")] = "1"
		data[which(data=="FALSE")] = "0"
		clin[[i]] = data
	}


for(i in 1:ncol(clin)){
	data = clin[[i]]
	temp = unique(data)
	temp = temp[!is.na(temp)]
	if(length(temp)==4 & sum(temp %in% c("0", "1", "FALSE", "TRUE")) ==4) {
	print(i)
	}
}











clin = read.table("tcga.clin.patient.full.txt")
clin = as.matrix(clin)
w = which(clin == "" | clin == "[Completed]" | clin == "[Discrepancy]")
clin[w] = NA

write.table(clin, "tcga.clin.patient.temp.txt", col.names=T, row.names=F, sep="\t", quote=F)

clin = fread("tcga.clin.patient.temp.txt")
class = sapply(clin, class)
names(clin) = tolower(names(clin))




action = array("", ncol(clin))
for(i in 34: ncol(clin)){
	cat(i, ":", names(clin)[i], ":\n")
	cat("  - class:", class[i], "\n")
	if( class[i] %in% c("integer", "numeric")) {
		table(clin[, i, with=F])
		cat("  - numeric: suggest to keep\n")
		default = "keep"
	} else if( class[i] == "logical") {
		print(table(clin[, i, with=F]))
		cat("  - logic: suggest to keep\n")
		default = "keep"		
	} else if( class[i] == "character") {
		t = table(clin[, i, with=F])
		if(length(t) > 20) {
			cat("more than 20 enums. suggest to ignore\n")
			default = "ignore" 
		}
		else {
			print(table(clin[, i, with=F]))
			cat("  - categorical: suggest to keep\n")
			default = "keep"
		}
	} else {
		cat("class type is", class[i], "\n")
	}	
	line = readline(prompt="Press [enter] to accept suggestion, [k] to keep, [i] to ignore, [r] to recode: \n")
	if(line == "k") default = "keep"
	if(line == "i") default = "ignore"
	if(line == "r") default = "recode"
	action[i] = default
}

clin = as.matrix(clin)
w = which(clin=="")
clin[w] = NA

recode = which(action=="recode" & class=="character")
for(i in 1: length(recode)){
	w = recode[i]
	temp = table(clin[, w])
	if(temp[order(-temp)][2] < 20) {
		action[w] = "ignore"
	}
}

recode = which(action=="recode" & class=="character")
for(i in 1: length(recode)){
	w = recode[i]
	temp = names(table(clin[, w]))
	if(length(temp)==4 & sum(temp %in% c("0", "1", "FALSE", "TRUE")) ==4) {
		data = clin[, w]
		data[which(data=="TRUE")] = "1"
		data[which(data=="FALSE")] = "0"
		clin[, w] = data
		action[w] = "keep"
	}
}

map = data.frame(variable=NULL, value=NULL, alias=NULL)
recode = which(action=="recode" & class=="character")
for(i in 1: length(recode)){
	w = recode[i]
	temp = names(table(clin[, w]))
	map = rbind(map, data.frame(variable = colnames(clin)[w], 
								value = names(table(clin[, w])), 
								alias = ""))
}
write.table(map, "map.txt", col.names=T, row.names=F, sep="\t", quote=F)

# manual editing
map = fread("map.txt")
for(i in 1: nrow(map)) {
	data = clin[, map$variable[i]]
	data[which(data==map$value[i])] = map$alias[i]
	clin[, map$variable[i]] = data
}
write.table(clin, "tcga.clin.patient.coded.all.txt", col.names=T, row.names=F, sep="\t", quote=F)
clin = clin[, which(action != "ignore")]
write.table(clin, "tcga.clin.patient.coded.pruned.txt", col.names=T, row.names=F, sep="\t", quote=F)
saveRDS(clin, "tcga.clin.patient.coded.pruned.rds")

clin2 = fread("tcga.clin.patient.coded.pruned.txt")
temp = NULL
class = sapply(clin, class)
w = which(class %in% c("logical", "character"))
for(i in 1: length(w)) {
	data = unlist(clin[, w[i], with=F])
	if(table(data)[order(-table(data))][2] < 20) {
		temp = c(temp, w[i])
	}
}




for(i in 1:ncol(clin)){
	data = clin[[i]]
	temp = unique(data)
	temp = temp[!is.na(temp)]
	if(length(temp)==4 & sum(temp %in% c("0", "1", "FALSE", "TRUE")) ==4) {
		data[which(data=="TRUE")] = "1"
		data[which(data=="FALSE")] = "0"
		clin[[i]] = data
	}


for(i in 1:ncol(clin)){
	data = clin[[i]]
	temp = unique(data)
	temp = temp[!is.na(temp)]
	if(length(temp)==4 & sum(temp %in% c("0", "1", "FALSE", "TRUE")) ==4) {
	print(i)
	}
}

