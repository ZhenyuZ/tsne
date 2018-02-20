# somehow refGenome library does not work. It's probably faster to parse myself instead of debug
library(data.table)

# Extract the index column of GTF attribute column
StripeValue = function(attribute, index, delim=" ") {
	keyvalue = sapply(attribute, function(x) strsplit(x, ";")[[1]][index])
	keyvalue = gsub("^ ", "", keyvalue)
	key = sapply(keyvalue, function(x) strsplit(x, delim)[[1]][1])
	print(table(key))
	value = sapply(keyvalue, function(x) strsplit(x, delim)[[1]][2])
	value = gsub("\"", "", value)
	return(value)
}

######################
# Gencode v22
######################
setwd("~/SCRATCH/tsne/meta/")

gtf.file = "gencode.v22.annotation.gtf"
length.file = "gencode.v22.gene.length.txt"

# read gtf
gtf = fread(gtf.file, h=F)
names(gtf) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# subset genes only, and extract attribute 
gene = gtf[feature == "gene"]
gene$gene_id = StripeValue(gene$attribute, 1)
gene$gene_type = StripeValue(gene$attribute, 2)
gene$gene_status = StripeValue(gene$attribute, 3)
gene$gene_name = StripeValue(gene$attribute, 4)

# add length
len = fread(length.file, h=T)
len$aggregate_length = as.integer(len$aggregate_length)
m = match(gene$gene_id, len$gene_id)
sum(is.na(m))
gene$exon_length = len$aggregate_length[m]

gene = gene[, -"attribute", with=F]
write.table(gene, "gencode.v22.genes.txt", col.names=T, row.names=F, sep="\t", quote=F)


######################
# miRBase v21 
######################
setwd("~/SCRATCH/tsne/meta/")

gtf.file = "hsa.gff3"

# read gtf
gtf = fread(gtf.file, h=F)
names(gtf) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# extract transcript and parsing
transcript = gtf[feature == "miRNA_primary_transcript"]
transcript$ID = StripeValue(transcript$attribute, 1, delim="=")
transcript$Alias = StripeValue(transcript$attribute, 2, delim="=")
transcript$Name = StripeValue(transcript$attribute, 3, delim="=")

transcript = transcript[, -"attribute", with=F]
write.table(transcript, "mirbase.v21.transcript.txt", col.names=T, row.names=F, sep="\t", quote=F)
