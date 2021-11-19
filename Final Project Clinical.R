library(TCGAbiolinks)
library(BiocManager)
clin_query <- GDCquery(project = "TCGA-KIRC", data.category = "Clinical", file.type = "xml")
#GDCdownload( clin_query ) #Only use this line ONCE! Comment out after you have downloaded the data. 
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
names( clinic )[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up" #

clinic_read_in <- read.csv("Downloads/kirc_clinical_data.csv", row.names = 1)

colnames(clinic)
unique(clinic$gender)


num_female <- sum( clinic$gender == "FEMALE" )
num_male <- sum( clinic$gender == "MALE" )


female_clinic <- clinic[ clinic$gender== "FEMALE", ]
male_clinic <- clinic[ clinic$gender == "MALE", ]

deseq_query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts")
#GDCdownload(deseq_query)

sum_exp <- GDCprepare(deseq_query)


if(!requireNamespace("HDF5Array")) BiocManager::install("HDF5Array")
library(HDF5Array)
BiocManager::install("SummarizedExperiment")
saveHDF5SummarizedExperiment(sum_exp, 
                             dir="htseq_h5_sumexp", 
                             prefix = "", 
                             replace = FALSE, 
                             chunkdim = NULL, 
                             level = NULL, 
                             as.sparse = NA, 
                             verbose = NA)
if (!require(DESeq2)) BiocManager::install("DESeq2")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
colnames(colData(sum_exp))

counts <- assays(sum_exp)$"HTSeq - Counts"

patient_data <- colData(sum_exp)

write.csv(data.frame(results_subset), "C:/Users/haleyboren/Downloads/results_subset.csv", )


patient_data$gender = factor (patient_data$gender, levels = c("female" , "male" ))

patient_data$gender

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patient_data,
                             design = ~gender)

dds_obj = DESeq(dds)
resultsNames(dds_obj)

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

plot(results$log2FoldChange, -log10(results$padj))
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

#make a volcano plot using ggplot()

volcano_plot = ggplot(data = data.frame(results), aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point() + 
  theme_minimal() + # make things pretty
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green")

volcano_plot

