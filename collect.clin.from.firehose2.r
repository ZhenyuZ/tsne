library(data.table)
clin = fread("tcga.clin.patient.full.txt")
clin[clin==""] = NA
