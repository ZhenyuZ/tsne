library(data.table)
library(ggplot2)
library(dplyr)
library(GGally)
library(survival)


d = readRDS("tcga.rna_mirna_met_cnv.combined.pc_tsne.rds")
signif = fread("tcga.rna_mirna_met_cnv_tsne.clin.significant.correlation.txt")

We examined correlation betweeen TSNE coordinates of TCGA molecular data (RNA expression, miRNA expression, Methylation and CNV) and clinical variables, and found both RNA and miRNA expression patten are highly correlated with Tumor_Status (TUMOR_FREE or WITH_TUMOR)

signif[project=="COAD", c("type", "tsne", "variable", "pval")]

Let's started with COAD miRNA expression data first
data = d[project == "COAD" & !is.na(tumor_status) & !is.na(mirna.tsne1) & !is.na(retrospective_collection)]
theme_set(theme_gray(base_size = 18))
ggplot(data, aes(mirna.tsne1, mirna.tsne2, col=tumor_status)) + geom_point()

Does it mean that miRNA data is a good predictor of tumor free status after surgery? 

We also observed correlation between the same molecular feature with retrospective_collection or prospective_collection. 
From Tara "These data were collected to see if there was a difference between the quality of samples that were collected with the TCGA protocols in mind (prospective cases) vs. the samples that were banked (sometimes for many years) and later consented for use in TCGA (retrospective cases). Let me know if that makes sense or if you have additional questions."

Ok, so is the correlation between TSNE features and tumor_status due to their correlation to retrospective_collection? 
xtabs(~ tumor_status + retrospective_collection, data)
            retrospective_collection
tumor_status   0   1
  TUMOR FREE 137  50
  WITH TUMOR  18 174
This is very obvious that retrospective_collection has a lot more 'with tumor' patient after surgery
  
> summary(lm(mirna.tsne1 ~ tumor_status + retrospective_collection, data=data))$coefficients
                         Estimate Std. Error   t value     Pr(>|t|)
(Intercept)              6.776069  0.3467193 19.543383 3.352022e-59
tumor_statusWITH TUMOR   2.739407  0.5725582  4.784504 2.466596e-06
retrospective_collection 3.351078  0.5822388  5.755504 1.794464e-08

How about subsetting data
summary(lm(mirna.tsne1 ~ tumor_status, data=data[retrospective_collection==1]))$coefficients
retrospective collection samples only 
> summary(lm(mirna.tsne1 ~ tumor_status, data=data[retrospective_collection==1]))$coefficients
                       Estimate Std. Error   t value     Pr(>|t|)
(Intercept)            9.055797  0.5555757 16.299843 8.322915e-40
tumor_statusWITH TUMOR 4.118616  0.6303658  6.533692 4.329094e-10

It's significant even with non-parametric test
np = wilcox.test(data[retrospective_collection==1 & tumor_status=="TUMOR FREE"]$mirna.tsne1, 
	 		data[retrospective_collection==1 & tumor_status=="WITH TUMOR"]$mirna.tsne1, 
	 		alternative = "two.sided")
np$p.val
2.398144e-09

prospective_collection samples only
> summary(lm(mirna.tsne1 ~ tumor_status, data=data[retrospective_collection==0]))$coefficients
                         Estimate Std. Error    t value     Pr(>|t|)
(Intercept)             7.1670725  0.3825806 18.7334985 1.918210e-41
tumor_statusWITH TUMOR -0.6275683  1.1226710 -0.5589958 5.769819e-01

So, correlation between TSNE features with tumor_status is gone in prospectively-collected samples. 


surv.model = with(data, Surv(calculated_death_or_last_contact_days_to/30.4167, vital_status=="Dead") ~ tumor_status)
ggsurv(survfit(surv.model), xlab = "Time (month)", ylab = "Survival Portion")
survdiff(surv.model)

coxph_obj = coxph(Surv(data$calculated_death_or_last_contact_days_to, data$vital_status=="Dead") ~ data$mirna.tsne1)




temp = data
temp$tumor_status = factor(temp$tumor_status)
temp$retrospective_collection = temp$retrospective_collection == 1
ggboxplot(temp, x = "tumor_status", y = "mirna.tsne1", color = "retrospective_collection", palette = c("#00AFBB", "#E7B800"))

temp$type = interaction(temp$retrospective_collection, temp$tumor_status)
ggplot(temp, aes(x=type, y=mirna.tsne1, col=type)) + geom_boxplot()
ggplot(temp, aes(x=mirna.tsne1, y=mirna.tsne2, col=type)) + geom_point()




