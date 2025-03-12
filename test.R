rm(list=ls()) 
source("DE_tools.R")

x<- "test_data/counts.tsv"
m<- "test_data/meta.tsv"

#bad example need to delete some columns

counts<-read.delim(x)
meta<-read.delim(m)

meta <- meta[-c(8,11:13),] #removing extra columns don't replicate

factor1<-colnames(meta)[6]
factor2<-colnames(meta)[8]

f<-""

# out<-DE_prep(counts,meta,f)
# working as intended
# 
# count_out<-out[[1]]
# meta_out<-out[[2]]
# genes<-out[[4]]


out<-DE_run(counts,meta)

 
# oh I literally chose the worst test case, this example data is randomly missing controls
# oh it might be read.delim just ignores empty values and they may have left it empty
# instead of righting control

#test:
# f=regular: works
# f=~1: works
# f="": works in that it throws an error because this is illegal, and it works if you simply don't assign f
# 
# x,m <- df: works
# x,m <- .csv: unknown
# x,m <- .tsv: works