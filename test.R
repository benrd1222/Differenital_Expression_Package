rm(list=ls())
source("DE_tools.R")

# new data
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")

# Download metadata
meta = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")

factor1<-colnames(meta)[10]

f<-sprintf("~%s",factor1)


# out<-DE_prep(counts,meta,f)
# 
# count_out<-out[[1]]
# meta_out<-out[[2]]
# formula_out<-out[[3]]
# genes<-out[[4]]
 
# working as intended, we get the counts and a metadata structure with a factor 
# of two levels which matches this experiment as a knockout experiment

out<-DE_run(counts,meta,f,g=T)

dds<-out[[1]]
genes<-out[[2]]

res = results(out[[1]], contrast=c("sample.characteristic.genotype.", "wildtypegenotype", "snai1knockout"), alpha=1e-5)
res

res_df<-DE_summary(res,genes,p=0.01)
head(res_df)


c_genes<-res_df[,1]
out<-DE_cluster(dds,counts,c_genes,heatmap = T)



# you can then use the output from DE_cluster to visualize over time if that was
# the contrast or just break your candidate genes into clusters for more visualization

# can also use the number of clusters you see in the dendrogram as starting guess
# for a k-means clustering approach... TOADD:

# x,m <- .csv: unknown

# this is a check against the data in the tutorial from OMGenomics
# https://github.com/omgenomics/youtube/blob/main/2024-12-17-deseq2-intro/tutorial.R
# res_df = as.data.frame(res)
# head(res_df) #row names exist
# head(genes)
# row.names(genes)<-genes[,1]
# genes<-genes[,-1,drop=FALSE]
# 
# res_df = merge(res_df, genes, by='row.names')
# head(res_df)
# 
# genes_to_check = c("thy1", "sfmbt2", "pasd1", "snai1")
# res_df[res_df[,8] %in% genes_to_check, ]

#they look exactly the same

# plotMA(res) #basic MA plot
# 
# BiocManager::install('EnhancedVolcano')
# library(EnhancedVolcano)
# EnhancedVolcano(res, lab=res_df$V1, x='log2FoldChange', y='pvalue') # more customizable plot


while(TRUE){
  usr_in <- readline(prompt = "Enter the number of clusters as visualized by the dendrogram: ")
  
  usr_in<-suppressWarnings(as.numeric(usr_in))
  if(!is.na(usr_in)){
    break
  }else{
    print("Please enter a numeric value")
  }
}

dds <- makeExampleDESeqDataSet(n=100,m=12)
dds$genotype <- factor(rep(rep(c("I","II"),each=3),2))

design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds) 
resultsNames(dds)

# the condition effect for genotype I (the main effect)
results(dds, contrast=c("condition","B","A"))

# the condition effect for genotype II
# this is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in genotype II compared to genotype I).
results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))

# the interaction term, answering: is the condition effect *different* across genotypes?
results(dds, name="genotypeII.conditionB")
