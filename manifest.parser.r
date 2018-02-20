library(data.table)


# rnaseq
cghub = fread("~/SCRATCH/tsne/meta/LATEST_MANIFEST_20160619.tsv")
rna.meta = fread("~/SCRATCH/tsne/meta/rnaseq_for_import.csv")
rna.meta = merge(cghub, rna.meta, by.x="analysis_id", by.y="cghub_analysis_id", suffix = c("", "_rna"))
write.table(rna.meta, "rna.meta.txt", col.names=T, row.names=F, sep="\t", quote=F)

# mirnaseq
mirna.meta = fread("~/SCRATCH/tsne/meta/mirnaseq_for_import.csv")
mirna.meta = merge(cghub, mirna.meta, by.x="analysis_id", by.y="cghub_analysis_id", suffix = c("", "_mirna"))
write.table(mirna.meta, "mirna.meta.txt", col.names=T, row.names=F, sep="\t", quote=F)