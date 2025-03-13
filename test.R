rm(list=ls()) 
source("DE_tools.R")

x<- "test_data/counts.tsv"
m<- "test_data/meta.tsv"

# bad example need to delete some columns with this jank data should really consider
# a better test.

counts<-read.delim(x)
meta<-read.delim(m)

meta <- meta[-c(8,11:13),] #removing extra columns don't replicate

factor1<-colnames(meta)[6]
factor2<-colnames(meta)[12]

f<-sprintf("~%s+%s",factor1,factor2)

# Main functionality test -----

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

# Visualization: never done unit testing before but would be better to have small
# test cases baked into like a vignette of the package
#
# something like
# test_DE_prep<-function(exampleargs)){
#   if(expected_results == DE_prep(exampleargs)){print("Working as intended")}
#   }
#
# Should probably read other packages to see how they handle testing but it's kind
# of fun to wing it and see how others do it after

## Visualization Tests ----------

# assuming you have output
# out<-DE_run(x,m,f)

results(out,contrast("something useful with new dataset"))
# yeah it's time to find a new dataset that actually has some differences


