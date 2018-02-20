library(data.table)
source("../tsne.util.r")

setwd("~/SCRATCH/tsne/autosome")

##############################
# methylation
##############################
load("../meth/TCGA.methylation.beta.sd0.2.rda")
# filter for automasomal probes only
probe.definition = fread("~/SCRATCH/dev/met/met/hm450.manifest.hg38.gencode.v22.tsv")
automasomal.probes = probe.definition[V1 %in% paste0("chr", 1:22)]$V4
met3 = met2[rownames(met2) %in% automasomal.probes, ]

# pca
ptm <- proc.time()
pca = pca(met3)
proc.time() - ptm
save(pca, file="TCGA.met.autosome.pca.rda")

# t-SNE
library(Rtsne)
ptm <- proc.time()
t = tsne(pca$x, permute=1000, use=50)
proc.time() - ptm
# about 2 hours, seed = 318
save(t, file="TCGA.met.autosome.tsne.perm1000.rda")

clin = fread("../clin/tcga.clin.patient.full.txt")
t$aliquot = rownames(t)
t$patient = substr(t$aliquot, 1, 12)
patients = intersect(t$patient, clin$bcr_patient_barcode)
m = match(patients, t$patient)
t = data.table(t[m, ])
m = match(patients, clin$bcr_patient_barcode)
clin = clin[m, ]
t = cbind(t, clin[, c("gender", "project", "history_neoadjuvant_treatment", "vital_status", "tissue_source_site", "neoplasm_histologic_grade", 
		"calculated_death",  "tumor_tissue_site", "age_at_initial_pathologic_diagnosis", "calculated_death_or_last_contact_days_to", 
		"birth_days_to", "days_to_initial_pathologic_diagnosis", "initial_pathologic_dx_year", "histological_type", 
		"history_other_malignancy", "race", "tumor_status", "last_contact_days_to", "ethnicity", "ajcc_pathologic_tumor_stage", 
		"tobacco_smoking_history_indicator", "death_days_to", "clinical_stage", "menopause_status", "tobacco_smoking_pack_years_smoked",
		"targeted_molecular_therapy", "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc"), with=F])

write.table(t, "TCGA.met.autosome.tsne.pdata.perm1000.txt", col.names=T, row.names=F, sep="\t", quote=F)


##############################
# mirna
##############################
count = read.table("../mirna3/mirna3.count.txt", h=T, check.names=F)

# filter low expression and gender chr
gff = fread("../meta/mirbase.v21.transcript.txt")
automasomal.genes = gff[seqname %in% paste0("chr", 1:22)]$Name

count = count[, which(substr(colnames(count), 1, 4) == "TCGA")]
count = count[rownames(count) %in% automasomal.genes, ]
mean.expr = apply(count, 1, mean)
count = count[which(mean.expr >= 1), ]
save(count, file="TCGA.mirna.autosome.lowfilter.count.rda")
tpm = t(t(count)/apply(count, 2, sum))*1e6
save(tpm, file="TCGA.mirna.autosome.lowfilter.tpm.rda")
tpm = log2(tpm+1)
save(tpm, file="TCGA.mirna.autosome.lowfilter.log2tpm.rda")


# pca
ptm <- proc.time()
pc = pca(tpm)
proc.time() - ptm
save(pc, file="TCGA.mirna.autosome.pca.rda")

# t-SNE
library(Rtsne)
library(parallel)
ptm <- proc.time()
t = tsne(pc$x, permute=1000, use=50)
proc.time() - ptm
# about 5 hours, seed = 595
save(t, file="TCGA.mirna.autosome.tsne.perm1000.rda")

clin = fread("../clin/tcga.clin.patient.full.txt")
t$aliquot = rownames(t)
t$patient = substr(t$aliquot, 1, 12)
patients = intersect(t$patient, clin$bcr_patient_barcode)
m = match(patients, t$patient)
t = data.table(t[m, ])
m = match(patients, clin$bcr_patient_barcode)
clin = clin[m, ]
t = cbind(t, clin[, c("gender", "project", "history_neoadjuvant_treatment", "vital_status", "tissue_source_site", "neoplasm_histologic_grade", 
		"calculated_death",  "tumor_tissue_site", "age_at_initial_pathologic_diagnosis", "calculated_death_or_last_contact_days_to", 
		"birth_days_to", "days_to_initial_pathologic_diagnosis", "initial_pathologic_dx_year", "histological_type", 
		"history_other_malignancy", "race", "tumor_status", "last_contact_days_to", "ethnicity", "ajcc_pathologic_tumor_stage", 
		"tobacco_smoking_history_indicator", "death_days_to", "clinical_stage", "menopause_status", "tobacco_smoking_pack_years_smoked",
		"targeted_molecular_therapy", "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc"), with=F])

write.table(t, "TCGA.mirna.autosome.tsne.pdata.perm1000.txt", col.names=T, row.names=F, sep="\t", quote=F)


##############################
# rna
##############################
load("../rna/rna.count.rda")

# annotate sample barcode
meta = fread("../meta/rna.meta.txt")
m = match(colnames(count), meta$htseq_count_id)
colnames(count) = meta$barcode[m]
count = count[, which(meta$study[m]=="TCGA")]
save(count, file="TCGA.rna.count.rda")

# filter low expression and non-autosomal genes
mean.expr = apply(count, 1, mean)
gtf = fread("../meta/gencode.v22.genes.txt")
automasomal.genes = gtf[seqname %in% paste0("chr", 1:22)]$gene_id
count = count[which(mean.expr >= 10 & rownames(count) %in% automasomal.genes), ]
save(count, file="TCGA.rna.autosome.lowfilter.count.rda")
tpm = t(t(count)/apply(count, 2, sum))*1e6
save(tpm, file="TCGA.rna.autosome.lowfilter.tpm.rda")
tpm = log2(tpm+1)
save(tpm, file="TCGA.rna.autosome.lowfilter.log2tpm.rda")

# pca
ptm <- proc.time()
pc = pca(tpm)
proc.time() - ptm
saveRDS(pc, file="TCGA.rna.autosome.lowfilter.pca.rds")

# t-SNE
library(Rtsne)
library(parallel)
ptm <- proc.time()
t = tsne(pc$x, permute=1000, use=50, dims=2)
proc.time() - ptm
# about 2 hours, seed = 20
saveRDS(t, file="TCGA.rna.autosome.lowfilter.tsne.perm1000.rds")

clin = fread("../clin/tcga.clin.patient.full.txt")
t$aliquot = rownames(t)
t$patient = substr(t$aliquot, 1, 12)
patients = intersect(t$patient, clin$bcr_patient_barcode)
m = match(patients, t$patient)
t = data.table(t[m, ])
m = match(patients, clin$bcr_patient_barcode)
clin = clin[m, ]
t = cbind(t, clin[, c("gender", "project", "history_neoadjuvant_treatment", "vital_status", "tissue_source_site", "neoplasm_histologic_grade", 
		"calculated_death",  "tumor_tissue_site", "age_at_initial_pathologic_diagnosis", "calculated_death_or_last_contact_days_to", 
		"birth_days_to", "days_to_initial_pathologic_diagnosis", "initial_pathologic_dx_year", "histological_type", 
		"history_other_malignancy", "race", "tumor_status", "last_contact_days_to", "ethnicity", "ajcc_pathologic_tumor_stage", 
		"tobacco_smoking_history_indicator", "death_days_to", "clinical_stage", "menopause_status", "tobacco_smoking_pack_years_smoked",
		"targeted_molecular_therapy", "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc"), with=F])

write.table(t, "TCGA.met.autosome.tsne.pdata.perm1000.txt", col.names=T, row.names=F, sep="\t", quote=F)




##############################
# cnv
##############################
cnv = read.table("../cnv/TCGA.cnv.window_1M.txt", check.names=F)
segmean.sd = apply(cnv, 1, sd)
cnv = cnv[which(!is.na(segmean.sd)), ]

# pca
ptm <- proc.time()
pc = pca(cnv)
proc.time() - ptm
saveRDS(pc, file="TCGA.cnv.1m.pca.rds")

# t-SNE
library(Rtsne)
library(parallel)
ptm <- proc.time()
t = tsne(pc$x, permute=1000, use=50)
proc.time() - ptm
# about 5 hours, seed = 800
save(t, file="TCGA.cnv.window_1m.tsne.perm1000.rda")

clin = fread("../clin/tcga.clin.patient.full.txt")
t$aliquot = rownames(t)
t$patient = substr(t$aliquot, 1, 12)
patients = intersect(t$patient, clin$bcr_patient_barcode)
m = match(patients, t$patient)
t = data.table(t[m, ])
m = match(patients, clin$bcr_patient_barcode)
clin = clin[m, ]
t = cbind(t, clin[, c("gender", "project", "history_neoadjuvant_treatment", "vital_status", "tissue_source_site", "neoplasm_histologic_grade", 
		"calculated_death",  "tumor_tissue_site", "age_at_initial_pathologic_diagnosis", "calculated_death_or_last_contact_days_to", 
		"birth_days_to", "days_to_initial_pathologic_diagnosis", "initial_pathologic_dx_year", "histological_type", 
		"history_other_malignancy", "race", "tumor_status", "last_contact_days_to", "ethnicity", "ajcc_pathologic_tumor_stage", 
		"tobacco_smoking_history_indicator", "death_days_to", "clinical_stage", "menopause_status", "tobacco_smoking_pack_years_smoked",
		"targeted_molecular_therapy", "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc"), with=F])

write.table(t, "TCGA.cnv.window_1m.tsne.pdata.perm1000.txt", col.names=T, row.names=F, sep="\t", quote=F)


# merge data

load("TCGA.autosome.pc10.tsne2.perm1000.met_cnv_rna_mirna.rda")

cnv.pc = data.frame(cnv.pc) 
cnv.pc$cnv.tsne1 = cnv.tsne$feature_1
cnv.pc$cnv.tsne2 = cnv.tsne$feature_2
cnv.pc = cnv.pc %>% 
	mutate(cnv.aliquot = rownames(cnv.pc), 
			sample = substr(cnv.aliquot, 1, 15),
			bcr_patient_barcode = substr(cnv.aliquot, 1, 12)) %>%
	arrange(cnv.aliquot) %>% 
	filter(!duplicated(sample))
colnames(cnv.pc)[1:10] = paste0("cnv.pc", 1:10)


rna.pc = data.frame(rna.pc) 
rna.pc$rna.tsne1 = rna.tsne$feature_1
rna.pc$rna.tsne2 = rna.tsne$feature_2
rna.pc = rna.pc %>% 
	mutate(rna.aliquot = rownames(rna.pc), 
			sample = substr(rna.aliquot, 1, 15),
			bcr_patient_barcode = substr(rna.aliquot, 1, 12)) %>%
	arrange(rna.aliquot) %>% 
	filter(!duplicated(sample))
colnames(rna.pc)[1:10] = paste0("rna.pc", 1:10)

mirna.pc = data.frame(mirna.pc) 
mirna.pc$mirna.tsne1 = mirna.tsne$feature_1
mirna.pc$mirna.tsne2 = mirna.tsne$feature_2
mirna.pc = mirna.pc %>% 
	mutate(mirna.aliquot = rownames(mirna.pc), 
			sample = substr(mirna.aliquot, 1, 15),
			bcr_patient_barcode = substr(mirna.aliquot, 1, 12)) %>%
	arrange(mirna.aliquot) %>% 
	filter(!duplicated(sample))
colnames(mirna.pc)[1:10] = paste0("mirna.pc", 1:10)

met.pc = data.frame(met.pc) 
met.pc$met.tsne1 = met.tsne$feature_1
met.pc$met.tsne2 = met.tsne$feature_2
met.pc = met.pc %>% 
	mutate(met.aliquot = rownames(met.pc), 
			sample = substr(met.aliquot, 1, 15),
			bcr_patient_barcode = substr(met.aliquot, 1, 12)) %>%
	arrange(met.aliquot) %>% 
	filter(!duplicated(sample))
colnames(met.pc)[1:10] = paste0("met.pc", 1:10)

temp1 = merge(rna.pc, mirna.pc, by=c("sample", "bcr_patient_barcode"), all=T)
temp2 = merge(met.pc, cnv.pc, by=c("sample", "bcr_patient_barcode"), all=T)
data = merge(temp1, temp2, by=c("sample", "bcr_patient_barcode"), all=T)
data = merge(data, clin, by="bcr_patient_barcode", all=T)
data$sample_code = substr(data$sample, 14, 15)

data = data.table(data)[(project == "LAML" & sample_code == "03") | 
			(project == "SKCM" & sample_code == "06") |
			(!project %in% c("LAML", "SKCM") &  sample_code == "01")]
data = data[!is.na(sample)]

saveRDS(data, "tcga.rna_mirna_met_cnv.combined.pc_tsne.rds")
write.table(data, "tcga.rna_mirna_met_cnv.combined.pc_tsne.txt", col.names=T, row.names=F, sep="\t", quote=F)

# t-SNE of randome sampleing
# pca3 = pca(met3)
# tsne3 = tsne(pca3$x)
# tsne3.pdata = annotate.tsne(tsne3, aliquot)
# tsne2.pdata = annotate.tsne(tsne2, aliquot)
# p3 = ggplot(tsne3.pdata, aes(x=feature_1, y=feature_2, col=project)) + geom_point()

png("TCGA.met.autosome.tsne.pdata.perm1000.png", width=1200, height=800)
p = ggplot(t, aes(x=feature_1, y=feature_2, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()

