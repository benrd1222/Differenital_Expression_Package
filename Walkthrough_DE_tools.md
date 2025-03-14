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
use install() and loaded into projects with library(). The main point is
you now have access to the functions written in DE_tools.

If you don’t want to copy DE_tools into the project directory you can
store it somewhere else and source it using its absolute file path.

### Loading in data

In this example I show how to read in data from the European
Bioinformatics Institute, but you can read in your own data as .tsv or
.csv or even pass the filepath to the data to DE_run() or DE_prep. I
will cover filepaths later, for now I will work under the assumption
that you can read your data in.

For differential expression you should have two matrices: one of counts
and one of metadata information.

``` r
#Download Count Data
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")

# Download metadata
meta = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
```

``` r
# factor1<-colnames(meta)[10]
# 
# f<-sprintf("~%s",factor1)
# 
# out<-DE_run(counts,meta,f,g=T)
# 
# dds<-out[[1]]
# genes<-out[[2]]
# 
# res = results(out[[1]], contrast=c("sample.characteristic.genotype.", "wildtypegenotype", "snai1knockout"), alpha=1e-5)
# res
# 
# res_df<-DE_summary(res,genes,p=0.01)
# head(res_df)
# 
# c_genes<-res_df[,1]
# out<-DE_cluster(dds,counts,c_gene)
```

## Extra Considerations

### Alternate Function Use Cases

### More Function Detail

### A Note on Maintenance

I will not maintain this project (ensure it retains functionality into
the future), so if you run into an error that you just cannot solve it
may be because the dependancies are working differently in the future.
The requirements.txt file associated with this project will have the
version of R this was built for and the version of packages utilized so
you can roll your installs back and it should work.
