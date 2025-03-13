# TODO: update docstring of DE_prep and DE_run to include gene as something I am returning

# TODO: deal with the logic, probably save the lengths of f_fac and f_join to memory to imporve speed
# actually the only thing that was slowing this down was calling clean(counts), added an all=F option
# for future use (probably should turn this into the columns you would like cleaned)

# TODO: add some visualization functions

# TODO: review the S3 classes documentation, would be good for future projects

if (!require(DESeq2) | !require(dplyr) | !require(tidyr) | !require(stringr)) {
  writeLines(readLines("requirements.txt"))
  stop("You don't have all of the necessary libraries loaded... please see the
         requirements.txt file for details")
}

strip<-function(n){
  # just a nice alias for a very repeatable action
  # names<-c("a BuNCh oF", "dumB nAmes", "wiTh too much", "VaRIENce")
  # strip(names)
  # [1] "abunchof"    "dumbnames"   "withtoomuch" "varience"
  
  n <- gsub(" ","",tolower(n)) #strip of whitespace and put into lowercase
  return(n)
}

clean<-function(df,all=T){
  #takes in a dataframe and uses strip to clean it
  row.names(df)<-strip(row.names(df))
  colnames(df)<-strip(colnames(df))
  if(all){
  df <- lapply(df, function(x) sapply(x, strip))
  }
  return(as.data.frame(df))
}

DE_prep <- function(x,m,f=FALSE){
  ## IN: x- is a counts matrix of gene expression
  # required to have the first 2 columns be Gene.ID and Gene.Name in that order,
  # or a filepath to the counts matrix
  # 
  # y- is a metadata df required to have a first column named run, or a filepath to the metadata. 
  # Other entries in metadata can be anything but if a column needs to be merged 
  # for the formula passed to DESeq the user will need to use that before hand.
  # 
  # f- is the structure of the formula desire to be run in with DESeq()
  # default is "~1" when passed from DE_run. Want to consider allowing DE_prep()
  # to be run without a formula though.
  #
  # DOES: Some repetitive pre-processing of common count expression data with
  # the ability to customize the formula used with DESeq
  #
  ## OUT: returns a list containing counts, metadata, and formula as cleaned version for use with 
  # DE_run(). DE_run() calls DE_prep() but this allows DE_prep() to be used 
  # independently if one should want
  
  
  #convert matrix inputs into dataframes because that's easy
  if(is.matrix(x)){x<-as.data.frame(x)}
  if(is.matrix(m)){m<-as.data.frame(m)}
  
  #load in the count matrix
  if(!is.data.frame(x)){
    #takes in a variety of files and read them, commonly csv or tsv add in more if necessary
    if(length(grep("csv$",x)) != 0){ #checks if the assumed filepath ends in csv
      counts<-read.csv(basename(x))
    }
    if(length(grep("tsv$",x) != 0)){
      counts<-read.delim(x)
    }
    else{
      stop("the input for the count matrix is not a dataframe and is not a supported file type")
    }
  }
  else{
    counts<-x
  }
  
  colnames(counts)<-strip(colnames(counts)) #just stripping in steps because I don't actually want to clean the numerical entries because it changes types
  
  # check requirements that the first 2 columns be Gene.ID and Gene.Name
  if(colnames(counts)[1]!= "gene.id" | colnames(counts)[2] != "gene.name"){stop("the names of the first two columns of the count dataframe do not match the required gene.id and gene.name format")}
  
  #now that counts is clean and order how we expect
  genes<-cbind(counts[,1],counts[,2]) 
  
  row.names(counts)<-counts[,1]
  counts<-counts[,-c(1,2)]
  row.names(counts)<-strip(row.names(counts))
  
  #load in the metadata
  if(!is.data.frame(m)){
    #takes in a variety of files and read them, commonly csv or tsv add in more if necessary
    if(length(grep("csv$",m)) != 0){ #checks if the assumed filepath ends in csv
      meta<-read.csv(basename(m))
    }
    if(length(grep("tsv$",m)) != 0){
      meta<-read.delim(m)
    }
    else{
      stop("the input for the count matrix is not a dataframe and is not a supported file type")
    }
  }
  else{
    meta<-m
  }
  
  meta<-clean(meta)
  
  # check requirements of first column being called run
  if(colnames(meta)[1]!="run"){stop("the first column of your metadata isn't named run, make sure it is a column named run filled with the sample IDs")}
  
  row.names(meta)<-meta[,1]
  
  ##check you aren't going to get a different sized matrix error down the road
  if(!all(row.names(meta)%in% colnames(counts))){
    print(colnames(counts[row.names(meta)%in% colnames(counts),]))
    stop("The sample names in the counts dataframe do not match the sample names listed in the metadata dataframe. The sample IDs listed above are the culprits")
    }
  
  
  # TODO: there has to be a better way, consider writing the logic out on paper
  # right now it seems as though it will work
  
  ## formula logic
  if(!is.null(f) && f != FALSE && !is.na(f) && f != ""){
    
    if(!is.character(f)){stop("Your formula design must be a string")}
    
    f<-strip(f) # strip whitespace to match dataframe cleaning protocol
    f<-sub("^.*~", "", f) #get rid of everything before and including the ~ if it exists
    
    # grab the factors to use to construct the formula later
    f_fac<-unlist(str_split(f, "[+*:]"))
    
    f_join<-unlist(str_extract_all(f,"[+,*,:]")) #will have character(0) if no results
    
    # stop if given an invalid formula
    if(f_fac[1]!=1){
    if(!all(f_fac %in% colnames(meta))==1){stop("One or more of the factors in your expression does not match a valid column name in the provided metadata")}
    }
    
    #stop if they are trying to break things on purpose
    if(length(f_join)>length(f_fac)){stop("check your formula, you're trying to break things!")}
    
    # we have the pieces of the design separated as f_fac and f_join
    # where f_join should be length(f_fac)-1
    
    if(length(f_fac)>1){
      formula<-"~"
      for(i in 1:length(f_fac)){
        
        formula <- paste0(formula,f_fac[i])
        
        if(i!=length(f_fac)){
          
          formula <- paste0(formula,f_join[i])
          
        }
      }
    }else{
      formula<-paste0("~",f_fac[1]) #no extra symbols and only one factor (this also encompasses default pass f="~1")
    }
    
    formula<-as.formula(formula)
  }
  
  #use f_fac to easy filter the columns we would like to keep if we have a formula
  if(f != FALSE && f_fac[1] != "1"){
    meta <- meta[,colnames(meta) %in% f_fac, drop = FALSE]
  }
  
  meta <- as.data.frame(lapply(meta, factor)) #turn all the columns into factors
  
# the above specifically doesn't work when only one factor is retained. I believe
# because it is being coereced into a vector and losing the ability to have
  
  if(f != FALSE){return(list(counts,meta,formula,tolower(genes)))}
  else{return(list(counts,meta,tolower(genes)))}
}

DE_run <- function(x,m,f="~1",g=F,cut=10){
  ## IN: x- is a counts matrix of gene expression
  # required to have the first 2 columns be Gene.ID and Gene.Name in that order,
  # or a filepath to the counts matrix
  # 
  # y- is a metadata df required to have a first column named run, or a filepath to the metadata. 
  # Other entries in metadata can be anything but if a column needs to be merged 
  # for the formula passed to DESeq the user will need to use that before hand.
  # 
  # f- is the structure of the formula desire to be run in with DESeq(). assumed to run
  # without any structure to the experimental factors if not passed a formula
  #
  # note: if you want just the 
  #
  # g- is an optional indicator, triggering the attachment of the gene.name to gene.id
  # key matrix. defaults off
  #
  # cut- is the criteria for low gene counts, defaults to 10
  #
  # DOES: runs DE_prep and runs DESeq with the desired formula
  #
  ## OUT: returns a DESeq datastructure, or a list with a DESeq datastructure and genes ID matrix if g=True
  
  prepped <- DE_prep(x,m,f) #check for requirements, and clean the data set to return a list(counts, metadata, formula)
  
  counts<-prepped[[1]]
  meta<-prepped[[2]]
  formula<-prepped[[3]]
  if(g){genes<-prepped[[4]]}
  
  
  # we now have counts, metadata, and formula available to us

  dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=formula)
    
  dds <- dds[rowSums(counts(dds)) > cut, ] 
  
  dds <- DESeq(dds)
  
  
  
  if(g){return(list(dds,genes))}
  else{return(dds)}
  
}

DE_cluster <- function(counts,res,p=0.01,genes=FALSE,normalize=TRUE,ab=0.5,heatmap=FALSE){
  ## IN: 
  # counts- the count matrix as output by DE_prep, only numbers.
  # res- the results object as determined by the user with the dds output from DE_run
  # p- p-value to determine cutoff of candidate gene
  # genes- a list of candidate genes if already determined or you want more specific control
  #        of who is part of the clustering process
  # normalize- default is to normalize the count matrix as it is likely coming in
  #            straight from DE_run
  # abline- customize cutoff for heirarchical clustering denrogram
  # heatmap- optional way of presenting data
  #
  ## DOES: organizes your count matrix by normalizing if asked, creating candidate genes
  #        from a results object unless provided, and creates a cluster tree and returns
  #        relevant data to making said tree
  #
  ## OUT: returns a the output of hclust() after some data fanangiling
  
  
  # if not passed a set of candidate genes, determine the candidate genes using
  # the given results. I am not bothering to check whether the results passed
  # match the count count matrix passed... just don't do that... please
  if(genes==FALSE){
    candidate_genes <- res %>% 
      filter(padj < p) %>%
      pull(gene) %>%
      unique()  
  }
  
  #normalization
  if(normalize==TRUE){
    # I forget this right now but something over each column
    # norm_cts <- lapply(counts, use_function(x),(x,x-mean/sd)) probably looks vaguely like that
  }
  
  
  hclust_matrix <- norm_cts %>% 
    select(-gene) %>% 
    as.matrix()
  
  # ENSURE THAT THE ROWNAMES HAVE MADE IT THIS FAR OR ADD THEM BACK IN
  
  # get the entries as z-scores
  hclust_matrix <- hclust_matrix %>% 
    t() %>% 
    scale() %>% 
    t()
  
  # now we need distance, generally genomics data is used with correlation
  require(amap) # no need to reinvent the wheel but TODO: make sure this does what you think it does
  gene_dist <- Dist(hclust_matrix,method="pearson")
  
  gene_hclust <- hclust(gene_dist, method = "complete")
  
  plot(gene_hclust, labels = FALSE)
  abline(h = ab, col = "#E34234", lwd = 2) # add horizontal line to illustrate cutting dendrogram
  
  # have to determine your own clustering grouping likely doesn't need to be in this function
  # gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  #  enframe() %>% 
  #  rename(gene = name, cluster = value)
  
  return(gene_hclust)
}






