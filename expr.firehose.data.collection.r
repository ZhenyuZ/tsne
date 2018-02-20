# Download data from firehose
# while read -r line; do wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/$line/20160128/gdac.broadinstitute.org_$line.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz ; done <"../../project"
# while read -r line; do zcat gdac.broadinstitute.org_$line.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz | head -n+2 | tail -1 | cut -f1 > $line.line2; done <"../../project"
# separate two input files into f2 and f13
# while read -r line; do zcat gdac.broadinstitute.org_$line.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz | sed '2d' | head -n+20532 > $line.normalized.rsem.txt; done <"f2"
# while read -r line; do zcat gdac.broadinstitute.org_$line.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz | sed '1d;3d' > $line.normalized.rsem.txt; done <"f13"
# ls *normalized.rsem.txt > files


library(data.table)

############################################################
# RNA
############################################################

disease = fread("~/SCRATCH/tsne/project", h=F)$V1
data = NULL
aliquot = data.frame(matrix("", 0, 2))
names(aliquot) = c("disease", "aliquot")

for(d in disease) {
	print(d)
	f = paste0(d, ".normalized.rsem.txt")
	# need to remove some wired characters in the raw file
	tt <- tempfile(tmpdir = "/home/ubuntu/SCRATCH/temp") 
	system(paste0("tr < ", f, " -d '\\000' >", tt))
	my.data = fread(tt, h=T, sep="\t")
	genes = unlist(my.data[, 1, with=F])
	names(genes) = NULL
	my.data = as.matrix(my.data[, 2:ncol(my.data), with=F])
	rownames(my.data) = genes

	data = cbind(data, my.data)
	aliquot = rbind(aliquot, cbind(d, colnames(my.data)))
}

write.table(aliquot, "disease_aliquot.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(data, "tcga.gaf.rsem.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(data, file="tcga.gaf.rsem.rda")


############################################################
# miRNA (LAML only as GAII; others are HiSeq)
############################################################

# while read -r line; do wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/$line/20160128/gdac.broadinstitute.org_$line.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2016012800.0.0.tar.gz ; done <"../../../project"
# while read -r line; do wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/$line/20160128/gdac.broadinstitute.org_$line.Merge_mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2016012800.0.0.tar.gz ; done <"../../../project"


disease = fread("~/SCRATCH/tsne/project", h=F)$V1
data1 = NULL
data2 = NULL
aliquot1 = data.frame(matrix("", 0, 2))
names(aliquot1) = c("disease", "aliquot")
aliquot2 = aliquot1

for(d in disease) {
	print(d)
	f = paste0(d, ".mirna.gene.txt")
	my.data = fread(f, h=T, sep="\t")
	genes = unlist(my.data[, 1, with=F])
	names(genes) = NULL
	my.data = as.matrix(my.data[, seq(2, ncol(my.data), 3), with=F])
	rownames(my.data) = genes
	if (nrow(my.data) == 1046) {
		data1 = cbind(data1, my.data)
		aliquot1 = rbind(aliquot1, cbind(d, colnames(my.data)))
	} else if (nrow(my.data) == 705) {
		data2 = cbind(data2, my.data)
		aliquot2 = rbind(aliquot2, cbind(d, colnames(my.data)))
	} else {
		cat("not all data are of 1046 or 705 miRNAs")
	}
}

aliquot = rbind(aliquot1, aliquot2)
write.table(aliquot, "disease_aliquot.txt", col.names=T, row.names=F, sep="\t", quote=F)

write.table(data1, "tcga.mirna.count.1046.txt", col.names=T, row.names=T, sep="\t", quote=F)
write.table(data2, "tcga.mirna.count.705.txt", col.names=T, row.names=T, sep="\t", quote=F)

# some data have 1046 genes, while others only have 705
common.genes = intersect(rownames(data1), rownames(data2))
data1 = data1[match(common.genes, rownames(data1)), ]
data2 = data2[match(common.genes, rownames(data2)), ]
data = cbind(data1, data2)
write.table(data, "tcga.mirna.count.txt", col.names=T, row.names=T, sep="\t", quote=F)
save(data, file="tcga.mirna.count.rda")



