# TODO: update docstring of DE_prep and DE_run to include gene as something I am returning

# TODO: deal with the logic, probably save the lengths of f_fac and f_join to memory to imporve speed
# actually the only thing that was slowing this down was calling clean(counts), added an all=F option
# for future use (probably should turn this into the columns you would like cleaned)

# TODO: review how other packages handle dependencies

# TODO: review the S3 classes documentation, would be good for future projects

# TODO: when a factor has more than 2 levels the comparison requires a reference
# level to be specified to make results easy to use after passing to DE_run
# Therefore, have a check within DE_run that checks if the pass from DE_prep
# has >2 factors then display options to set reference and give numerical choice
# to choose which one is set as the reference
## WAIT: THERE'S A POTENTIAL ALL OF THIS WAS BASED ON A TYPO, yeah I just had a typo
# although being able to designate a reference level could be useful

# TODO: make DE_prep() deal with continuous metadata variables or note in the docstring
# it assumes and only works with categorical variables

# TODO: update the difference set to display to the user when they have mismatched
# colnames and rownames. My previous logic doesn't displat true difference

if (!require(DESeq2) | !require(dplyr) | !require(tidyr) | !require(stringr)) {
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
  
  meta <- meta %>% mutate_all(factor)
  

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
  
  # realizing that important step here is to also allow the user to enter
 # if(length(dds$condition)>2){ ask the user for the reference condition and
  # set the input to the reference level}
 
  dds <- DESeq(dds)
  
  
  
  if(g){return(list(dds,genes))}
  else{return(dds)}
}

DE_summary <- function(res,genes,p=FALSE){
  ## IN: 
  #
  # res- the results object as determined by the user with the dds output from DE_run
  # p- p-value to determine cutoff of candidate gene, this also indicates you
  #    handed DE_viz to handle the full list of genes and geneIDs to filter the results object
  #
  # genes- a list of candidate genes if already determined or the full genes matrix
  #         both need to be a dataframe with column 1 being ids and column 2 being names
  #     note:
  #     genes is not a required argument and DE_summary(res, p) will just filter
  #     on the p-value but it will not append gene names to the data frame
  #
  ## DOES: creates candidate genes or takes in candidate genes and organizes the
  #        results into a useful dataframe. Essentially just a quick filtering function
  #        based on a p-value or a vector of gene names
  #        
  #
  ## OUT: filters by p-value of results if given a p-value or filters by listed genes
  #       otherwise
  
  if(class(res)!="DESeqResults"){stop("Please ensure that res is of class DESeqResults")}
  
  res_df <- as.data.frame(res)
  
  if(p!=FALSE){
    res_df <- res_df %>%
      filter(padj < 0.01)
  }
  
  if(!missing(genes)){
    #quick checks
    tryCatch({colnames(genes)<-c("gene.id","gene.name")},error=function(e){
      message("Make sure genes=df[ids,names]")
    })
    
    
    
    # making genes a vector for merging by rowname
    rownames(genes)<-genes[,1]# consider adding this to the DE_prep flow
    genes<-genes[,-1, drop=FALSE]
    
    
    
    #otherwise we assume that the merging will filter with all = FALSE dropping unmatched ids
    res_df = merge(res_df, genes, by='row.names', all=FALSE)
  }
  
  return(res_df) # the fact I'm going to return so much means this may not make sense
}

DE_cluster <- function(dds,counts,c_genes,normalize=TRUE,cut=TRUE,heatmap=FALSE){
  ## IN:
  # dds: the DESeqDataSet for normalization
  # counts- the raw count matrix or the normalized count matrix
  # c_genes- the candidate genes. rownames() of res_df as output by DE_summary()
  # normalize- default is to normalize the count matrix as it is likely coming in
  #            straight from DE_run, but you can give a normalized count matrix
  #            just set normalize to FALSE (also means you don't need to pass a dds)
  # cut- set to let user input after seeing the dendrogram, input a number before
  #      to disable the user input from terminal
  # heatmap- optional way of presenting data
  
  #nope this is going to be less user friendly
  # we will be assuming the raw counts are being submitted (raw data)
  # the dds that the raw counts were used on (DE_run)
  # the candidate_gene id's (DE_summary)
  
  require(amap)
  
  # it was going to be annoying to stack this in to one visualization call and also
  # probably just not helpful, so I'm not
  stopifnot(is.logical(heatmap))
  stopifnot(class(dds)=="DESeqDataSet")
  if(is.null(dds) && normalize==TRUE){stop("please provide a DESeqDataSet to estimate size factor for normalization")}
  
  
  #DESEQ has a normalization function rlog() that we could use instead
  if(normalize==TRUE){
    #dealing with raw counts if we need to normalize
    rownames(counts)<-strip(rownames(counts))
    colnames(counts)<-strip(colnames(counts))
    stopifnot(colnames(counts)[1]=="gene.id")
    rownames(counts)<-counts[,1]
    counts[,-c(1,2)]
    
    sf<-estimateSizeFactors(dds)
    hclust_matrix <- as.data.frame(counts(sf, normalized=TRUE))
  }
  else{
    #assume if they turned normalize false then counts is as they want
    stopifnot(all(sapply(counts, is.numeric))) # at least check they passed a matrix
    hclust_matrix <- counts
    }

  hclust_matrix<-hclust_matrix[rownames(hclust_matrix) %in% c_genes,]
  
  hclust_matrix <- hclust_matrix %>% 
    t() %>% 
    scale() %>% 
    t()
  
  gene_dist <- amap::Dist(hclust_matrix, method="pearson")
  
  gene_hclust <- hclust(gene_dist, method = "complete")
  
  plot(gene_hclust, labels = FALSE)
  abline(h = 0.5, col = "brown", lwd = 2) #make the height definable maybe???
  
  # user will determine the amount of clusters as determined by visual analysis of
  # the dendrogram
  
  while(TRUE && cut==TRUE){
    usr_in <- readline(prompt = "Enter the number of clusters as visualized by the dendrogram: ")
    usr_in<-suppressWarnings(as.numeric(usr_in))
    
    if(!is.na(usr_in)){
      break
    }else{
      print("Please enter a numeric value")
    }
  }
  
  if(cut==TRUE){ #if defaulted take the input
    cut <- as.numeric(usr_in) 
    
  }else{ #check input: have to check here because this is if input is passed as argument rather than live input
    stopifnot(is.numeric(cut))
    }

  
  #heatmap is working, would like to have the gene names instead of IDs
  if(heatmap==TRUE){
    suppressPackageStartupMessages(library(ComplexHeatmap))
    require(ComplexHeatmap)
    plot(Heatmap(hclust_matrix, show_row_names = FALSE))
  }
  
  return(as.data.frame(cutree(gene_hclust, k = cut)))
  
}
