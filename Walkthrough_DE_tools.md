How to Use DE_tools
================

- [Setup](#setup)
- [Full Example](#full-example)
- [Extra Considerations](#extra-considerations)

This is a document that will walk you through using my collection of
functions DE_tools that handles the data organization and some repetitve
tasks. The goal is to make differential expression analysis approachable
for someone totally new to R and to make the repeated analysis of RNASeq
data easier.

All you need is the count matrix, the metadata matrix, and a design
formula for the GLM in the back end.

My goal is to provide a brief walkthrough of differential expression
analysis using these tools, the different ways in which the functions
can be used to break off from the example, and to provide a list of the
functions and a brief overview of their intended use and ways to
personalize further to your project.

## Setup

In this section we will discuss the basics of installing R and starting
a project, if you already know how to do this navigate to the Full
Example.

### Installing R and RStudio

Navigate to <https://cran.r-project.org/bin/windows/base/> for the
windows download or to <https://cran.r-project.org/bin/macosx/> for the
mac download. This is the the download for the actual R language and
software. You could do the analysis through the R shell but that would
be painful.

Rstudio is the preferred integrated development environment for R and
can be downloaded here <https://posit.co/download/rstudio-desktop/>.
This is what most people mean when they refer to writing scripts or
running analyses in R, it provides an nice graphical user interface for
coding in R.

### Setting up your project

Before doing anything it will be helpful if you create a new project. In
your R-studio window navigate to file -\> New Project. This will bring
up a window prompting you with various ways to create a project. You can
create a new directory (folder) to start your project in or choose an
existing directory, for example a project folder on your Desktop with
the data for your analysis stored in it.

<figure>
<img src="walkthrough_files\new_proj.png" alt="Should look like this" />
<figcaption aria-hidden="true">Should look like this</figcaption>
</figure>

### A note on the working directory

Below is a bit of a lengthy tangent to try and give you context into how
your machine handles files which is useful but not critical to
understand.

Here’s the main take away:

It is important to know where R is working when performing an analysis,
both for reading and writing data.

If you see this,

<figure>
<img src="walkthrough_files\no_project.png"
alt="No project highlighted in red" />
<figcaption aria-hidden="true">No project highlighted in
red</figcaption>
</figure>

or check your working directory like this

![](walkthrough_files\getwd.png) and it doesn’t end in the name of your
project folder.

Set your working directory to the project directory like this,

![](walkthrough_files\set_directory.png)

#### A Tangent About File Paths

When you consider files on your computer you may think that they are
relatively easy to find. You probably have some system of organization
that includes grouping by related projects or file types. In reality you
have dense network of directories or folders each with files and
directories of their own all anchored to some root directory (/). In a
windows systems this root is usually the C drive (C:\\, which you can
see for yourself if you hit Win+R, type cmd, and hit enter to open your
command terminal. The full file path that pops up probably looks
something like this C:\Users\usr_name. If you then type tree and hit
enter you will see a lot folders and files connected with lines. These
are all of your files and representations of the paths from the current
working directory (C:\Users\usr_name or whatever is printed before the
\> in your terminal).

If you are on mac and want to follow along hit command+space, search
terminal, and type in ls and hit enter in the terminal. This only
displays the directories directly in the current working directory so
won’t be as dramatic.

Every path presented here has a relative path based on the current
working directory which is .\path\to\file where the . represents the
current working directory. However every directory and file on your
computer also has an absolute path listing every directory you need to
walk through from the root.

![](walkthrough_files\simple_tree.png)

For example if we added some example files to our newly created R
project a the tree from the project directory, indicated by the . after
C:, our tree might look like this. To access these files we now only
need to take one step along the path, which in this case just means type
out the name of the file of interest. However the absolute file path for
any of these file will be much longer walking all the way from C:\\
through Users\\ into \usr_name likely into \Onedrive where your \Desktop
is stored until finally we reach \example_project. For DE_tools.R the
relative path would be .\DE_tools.R while the absolute file path would
be C:\Users\usr_name\Onedrive\Desktop\example_project.

Your working directory in R is just the same as above, a starting point
for walking to local files without listing all of the directories from
the root. This also means if you need to access something like data
elsewhere than your R project you can copy the absolute file path.

Similarly important, is that R uses the working directory as the
location to store environmental variables to your .Rhistory when you
close a session and choose to save this information. So if you were to
process a bunch of projects in the default R working directory you could
have conflicting variables loaded in and overwrite or use incorrect data
if you swapped between projects without realizing where you were
working.

Quick Note: Windows does file paths with the backslash \\ and nearly
everything else uses \\ as an escape character to signal a special
symbol (even writing this right now I have to write \\ to display one
backslash). This means when using an absolute windows file path you need
to write \\ for every \\ in the file path or better yet just ignore
window’s wishes and use forward slashes (/). Otherwise, R will throw a
fit as it almost immediately assumes you are trying and failing to write
a hexadecimal digit (\u{HHH}) when it hits the the \U in \Users and just
assumes you forgot to give it numbers.

## Full Example

### Getting Access to the tools

The next step is to copy DE_tools.R into the project directory. This
makes obtaining the function stored in the tools file easier. In the end
your project directory should look something like the following:

![](walkthrough_files\example_directory_startup.png)

Now we can open up a new file to begin a differential expression
analysis using DE_tools.

However to have access to any of the tools within the file you need to
source it. I suggest doing this at the top of the file to ensure you
don’t try to use a function before it is read in. If we place DE_tools.R
into our project directory and have set our working directory then we
can just run the following.

``` r
source("DE_tools.R")
```

This is similar to using a package except packages are stored on The
Comprehensive R Archive Network (CRAN) and downloaded locally when you
use install() and the package functions are added to package
environments with library(). The main point is you now have access to
the functions written in DE_tools as they have been added to your global
environment.

``` r
# source("C:/absolute/filepath/to/DE_tools.R")
```

You can source from an absolute file path if you want to store
DE_tools.R elsewhere, if you don’t want to copy and paste it into each
project.

Sourcing this you may have hit an error telling you that you don’t have
the required packages installed. Read through the requirements.txt file
to see which packages you need and use the Tools -\> Install Packages in
R Studio to download dependencies.

### Loading in data

In this example I show how to read in data from the European
Bioinformatics Institute.

For differential expression you should have two matrices: one of counts
and one of metadata information.

``` r
#Download Count Data
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")

# Download metadata
meta = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
```

Let’s take a look at the structure of each of these dataframes.

``` r
head(counts[,-8]) #I'm cutting off the last column to fit the output within the margins
```

    ##           Gene.ID Gene.Name ERR1736465 ERR1736466 ERR1736467 ERR1736468
    ## 1 ENSG00000000003    TSPAN6        434        829        334        779
    ## 2 ENSG00000000005      TNMD          0          0          0          0
    ## 3 ENSG00000000419      DPM1       2202       3461       1627       4241
    ## 4 ENSG00000000457     SCYL3         45         72         31        119
    ## 5 ENSG00000000460  C1orf112        320        558        235        516
    ## 6 ENSG00000000938       FGR          0         10          0          0
    ##   ERR1736469
    ## 1        578
    ## 2          0
    ## 3       3676
    ## 4         69
    ## 5        468
    ## 6          0

Counts is a datafrane of the genes and the the frequency of occurence in
each sample. The important thing to note is that counts has a column of
gene ids and a column of gene names as the first two columns, yours
should be in this order and similarlry named gene.id and gene.name.

``` r
head(meta[,10:11])
```

    ##   Sample.Characteristic.genotype. Sample.Characteristic.Ontology.Term.genotype.
    ## 1              wild type genotype                                            NA
    ## 2              wild type genotype                                            NA
    ## 3              wild type genotype                                            NA
    ## 4                  Snai1 knockout                                            NA
    ## 5                  Snai1 knockout                                            NA
    ## 6                  Snai1 knockout                                            NA

The metadata is the relevant descriptors of cell line, species,
ontology, and treatment groups for each sample in the run column. The
only important ordering here is that the first column contains the
sample names and the column is labeled run.

### Determining Your Formula

That is very nearly the only things you need to do and hopefully your
data was already pre-delivered in said format. The last thing we need to
do is to determine a formula for the generalized liner model running in
the background. This is very dependent on the question you are asking
and your experimental set up.

For this specific data set the experiment was to compare wild type to
Snai1 knockouts in breast cancer cells.

``` r
meta$Sample.Characteristic.genotype.
```

    ## [1] "wild type genotype" "wild type genotype" "wild type genotype"
    ## [4] "Snai1 knockout"     "Snai1 knockout"     "Snai1 knockout"

This means we care about testing for the effect of genotype. The formula
would be written something like this ~genotype. You may have more
treatments or groups in your own data and should consider carefully how
your design formula is impacting your tests. For example if you had an
experiment with a group variable and a set of treatments applied to all
groups then if you care about the full set of interactions on gene
expression your formula would be ~group\*treatment.

For the use case of DE_tools you simply need to provide the names of the
columns of importance in your metadata combined with +, \*, or : which
are used to define formulas in R.

``` r
# you can hand write your equation or you can store the column names you care
# about and then create the final formula with the fstrings implementation in R
# using sprintf()

factor1<-colnames(meta)[10]

f<-sprintf("~%s",factor1) 
f
```

    ## [1] "~Sample.Characteristic.genotype."

### Running the Analysis

We now have everything we need to call DE_run.

``` r
out<-DE_run(counts,meta,f,g=T)

dds<-out[[1]]
genes<-out[[2]]
```

Making the results out of the output of DESeq takes one more step
related to your experimental design. You can use the results function
from DESeq to determine what contrast you would like to extract from the
DESeq object.

The easiest way to visualize your options is to get resultsNames(dds)
and read the documentaion on the results() function in DESeq2.

``` r
resultsNames(dds)
```

    ## [1] "Intercept"                                                        
    ## [2] "sample.characteristic.genotype._wildtypegenotype_vs_snai1knockout"

``` r
# ?DESeq2::results() call the man page for DESeq2's results() function
```

In this example we only have a main effect for genotype, which we should
know from our setup of our formula. However, as you can see above we can
also see by listing the different effects modeled by the GLM with
resultsNames(dds).

``` r
# Therefore, we would like to use the contrast statement with results
res = results(dds, contrast=c("sample.characteristic.genotype.", "wildtypegenotype", "snai1knockout"), 
alpha=1e-5)
```

This is the raw results containing all of the significant and
non-significant comparisons of the main effect. I also provide and
DE_summary function that quickly filters the results dataframe by the
p-value cutoff you desire. It also optionally takes in a gene matrix of
the gene names and gene ids and appends the relevant gene names.

``` r
res_df<-DE_summary(res,genes,p=0.01)
head(res_df)
```

    ##         Row.names  baseMean log2FoldChange     lfcSE      stat       pvalue
    ## 1 ensg00000001084 2025.6506      0.6079598 0.1215264  5.002700 5.653297e-07
    ## 2 ensg00000001460  118.5259     -9.9919973 1.2042907 -8.296998 1.067754e-16
    ## 3 ensg00000001461  518.9701     -1.3789857 0.1440883 -9.570419 1.064730e-21
    ## 4 ensg00000001617   61.0675     -2.9498782 0.4592325 -6.423497 1.331790e-10
    ## 5 ensg00000001629 1286.0658      0.4983303 0.1535889  3.244572 1.176273e-03
    ## 6 ensg00000002330 1089.7154     -1.1370367 0.1178679 -9.646702 5.076181e-22
    ##           padj gene.name
    ## 1 5.639037e-06      gclc
    ## 2 4.168661e-15     stpg1
    ## 3 6.559060e-20    nipal3
    ## 4 2.496393e-09    sema3f
    ## 5 5.218940e-03    ankib1
    ## 6 3.258058e-20       bad

The difference in expression can subsequently be visualized with the
using the res and res_df dataframes.

### Visualization Options

Visualization of the differential expression depends greatly on your
experimental design and personal preference.

A general visualization can be made with MAplot for the whole results
dataframe.

``` r
# plotMA() is a more basic version of the enhanced volcano plot below

library(EnhancedVolcano)

EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue') # another useful function for visualization is the EnhancedVolcano plot
```

![](Walkthrough_DE_tools_files/figure-gfm/volcano-1.png)<!-- -->

``` r
# note: would be better to get the row names as the gene names instead of the gene ids for the
# res data frame

# PCAs are another option
```

The following is a tool I made to use hclust() to allow you to look at a
denrogram and then deside how many clusters you think there are and then
group the ids into clusters for further use.

``` r
c_genes<-res_df[,1]
cluster_df<-DE_cluster(dds,counts,c_genes,cut=8)
```

![](Walkthrough_DE_tools_files/figure-gfm/cluster-1.png)<!-- -->

``` r
# cut is an alternative parameter that circumvents the interactive input of
# this function cut=8 means I am setting the cluster size before seeing the dendrogram
```

    ##                 cutree(gene_hclust, k = cut)
    ## ensg00000001084                            1
    ## ensg00000001460                            2
    ## ensg00000001461                            2
    ## ensg00000001617                            3
    ## ensg00000001629                            4
    ## ensg00000002330                            2

To Add:

PCA and k-means clustering

## Extra Considerations

### Alternate Function Use Cases

A few of the functions have alternate use cases described more clearly
in the doc-strings (descriptions under the function declaration) within
DE_tools

- DE_prep is the function used within DE_run to data wrangle, it can be
  used to just help structure the data and then you can customize your
  own script from there.

- DE_summary can just take a gene input if you have certain genes you
  ineterested in filtering out of your results. This isn’t hard to do
  outside of this function its just an extra feature because I coded it
  flexibly

- DE_cluster is meant to be used without a cut parameter to allow one to
  look at the dendrogram and then make a decision

### A Note on Maintenance

I will not maintain this project (ensure it retains functionality into
the future), so if you run into an error that you just cannot solve it
may be because the dependancies are working differently in the future.
The requirements.txt file associated with this project will have the
version of R this was built for and the version of packages utilized so
you can roll your installs back and it should work.
