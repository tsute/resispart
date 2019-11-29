<h2>Table of Content</h2>

* [Part I Introduction](part1)
* [Part II Workshop Tutorial](part2)
   * [1. Introduction](part2/#1)
   * [2. Installation of the R software and packages for analyzing microbiome sequence data](part2/#2)
   * [3. Basic R operations](part2/#3)
   * [4. Use of the R DADA2 package to process the 16S rRNA gene sequence reads](part2/#4)
      [Step 1 Create a working directory and download a set of demo Illumina 16S rRNA sequences](part2/#4-1)
      [Step 2 Find the sequence data](part2/#4-2)
   * [3. Basic R operations](part2/#3)
   * [3. Basic R operations](part2/#3)
   * [3. Basic R operations](part2/#3)  


<h2 style="font-weight:bold"> Part I. Introduction - how to use this workshop tutorial</h2>

<h3 style="font-weight:bold">1. R Code Box</h3>

The text in the light-blue box, contains the R codes that you can highlight (by dragging your mouse), copy (use the control-C key combination) and paste (control-V) into your R command-line interface on your computer. For example:

``` R
#Below is a working R code that you can copy and paste into your R to execuate some R commands
message("Hello! Welcome to the RESISTPART Bioinformatics Workshop")
2+2
3*3
```

<h3 style="font-weight:bold"> 2. Alternative method box </h3>

Text inside the alternative method box, shows another or more different ways of doing the same thing:

<table bgcolor="#C8FBD3"><tr><td style="font-size:0.8em" width="100%">
<b>Alternative:</b><br>
This is where other ways of doing the same thing will be described
</td></tr></table>

<h3 style="font-weight:bold"> 3. Notes and comments </h3>

<table  bgcolor="#DDDDDD"><tr><td style="font-size:0.8em;font-style:italic;">
  <b>Note:</b> A note is something that needs to be mentioned but is apart from the context.
</td></tr></table>



<br><br><br>

<h2  style="font-weight:bold;color:navy"> Part II. Workshop Tutorial</h2>

<br><br>

<h3 style="font-weight:bold;color:navy"> 1. Introduction</h3>

The goal of this workshop is to gain hands-on experience for analyzing the 16S rRNA gene sequences generated by the next generation sequencing (NGS) platform, such as the Illumina sequencers. The 16S rRNA gene is commonly used by the research community as a marker to decipher the diversity of the microbial community in a sample. The reasons for this are first, 16S rRNA gene is universally present in all prokaryotes (inlcuding Bacteria and Archaea); secondly, although the 16S rRNA gene sequences are very conserved in both Bacteria and Archaea, there are enough variability within the genes (i.e,. 9 hypervariable regions) for different species. Hence by sequencing the 16S rRNA genes in a microbial community, and matchinng the sequences to a set of 16S rRNA gene sequences with know taxonomy information, one can determine how many and what are the species in the samples. Since the NGS has a limitation of how long a DNA sequence can be decoded, only a portion of the 16S rRNA gene can be targeted. Usually the first step is to PCR-amplify a region of the 16S rRNA gene from the DNA extracted from the samples. The amplicons are then sequenced by the NGS technology to produce millions of DNA reads. These sequence reads are then filtered and denoised to reduce the errors and artifacts. The filtered reads are then taxonomically aasigned to determine their genus or species. The outcome of these processes is a read count table with columns representing samples and rows representing organisms, either genus, species or OTUs - operational taxonomic units, if the taxnomy has yet to be determined). This table is usually called an OTU table.

In this workshop we will learn how to use several commonly used software tools to analysze a set of demo sequences. After this workshop, we will have pratical hand-on experience for the bioinformatics task that is required to process this type of data in order to understand the mirobial diversity of the samples. 

<br><br>

<h3 style="font-weight:bold;color:navy"> 2. Installation of the R software and packages for analyzing microbiome sequence data </h3>

R is short for "The R Project for Statistical Computing" and is a free software environment for statistical analysis and data visualization. R has been a very popular software development platform for many bioinformatics toosl, including those that are useful for analyziging NGS sequenes reads, as well as the microbiome diversity analysis. 

In this workshop we will use several sotftware packges developed on the R platform. Henct the first thing is to download and install R on your computer.

For this workshop we will pratice the analysis on Widnows 10 Pcs so we can download the latest R for Windows computers from this link:

<https://vps.fmvz.usp.br/CRAN/bin/windows/base/R-3.6.1-win.exe>

click above link and save the file to the default download location on your computer. Find the downloaded file, double-click to install R on your computer. During the installation just respond with default options.

<table bgcolor="#C8FBD3"><tr><td style="font-size:0.8em">
<b>Alternative:</b><br>
These are some other mirrored web sites that provide this download:<br><br>

Brazil<br>
<a href="https://nbcgib.uesc.br/mirrors/cran/">	Computational Biology Center at Universidade Estadual de Santa Cruz</a><br>
<a href="https://cran-r.c3sl.ufpr.br/">	Universidade Federal do Parana</a><br>
<a href="https://cran.fiocruz.br/">	Oswaldo Cruz Foundation, Rio de Janeiro</a><br>
<a href="https://vps.fmvz.usp.br/CRAN/"> University of Sao Paulo, Sao Paulo</a><br>
<a href="https://brieger.esalq.usp.br/CRAN/">	University of Sao Paulo, Piracicaba</a><br>
<br><br>
For a complete list of all the mirrored sites in different country, visit this link:<br>
<a href="https://cran.r-project.org/mirrors.html">https://cran.r-project.org/mirrors.html</a>
  
</td></tr></table>
<br>
<table bgcolor="#C8FBD3"><tr><td style="font-size:0.8em">
You can also downlaod another version of R, Microsoft R Open here:<br>
<a href="https://mran.microsoft.com/open">https://mran.microsoft.com/ope</a><br>
It's an enhanced version of the original R. The most important feature of R Open is that it can use the multiple computer CPU cores hence it can speed up many process during the analysis. The R version for Windows can only use a single core hence it can be very slow, especially with a larger (sequence) data set.
</td></tr></table> 

After R is installed on your computer, find the shortcut and launch the program. You should choose the "x64" version which is fastser. You will be presented with an R graphical interface:
<img src="https://i.gyazo.com/07ae9074dbbd826b3e741ca3407dbd97.png">

There are many online R tutorials for beginners, including many Youtube videos, such as:

<https://www.youtube.com/watch?v=iijWlXX2LRk>

Next we'll get our hands wet by testing out the R in Windows.

<br><br>

<h3 style="font-weight:bold;color:navy"> 3. Basic R operations </h3>

<h4 style="font-weight:bold;color:navy"> The R command-line interface </h4>

There are two types of input in the R command-line interface 
1) comments - any line that start with a # sign, is treated as comment and the R program will ingore the entire line.
2) commands - if a line does not beging iwth a # sign, it is treated as a command and whatever you type in must obey the R language and grammer. If you type something that R doesn't recognize, it will give you an error output.
3) multipe lines command input - if a line ends with a / (backward slash) sign, it tells the R program that you have done giving out your commands, and will continue your command input in the next line. R will then wait until the last line without the / sign and collect all the lines together and execuate the commands. This is because often time a command has many parameters, or a long parameter (such as a long flie name) and this multi-line command featuer will come in handy.

Now copy and paste the example R code below into your R command line interface:
``` R
#Below is a working R code that you can copy and paste into your R to execuate some R commands
message("Hello! Welcome to the RESISTPART Bioinformatics Workshop")
2+2
3*3
```

You should get an output like this:
<img src="https://gyazo.com/bcbdbca157c0f9e720a2a63a44844a6d.png">


<h4 style="font-weight:bold;color:navy"> Using R as a calculator</h4>

```R
# BASIC ARITHMETIC OPERATORS
2-5                # subtraction
6/3                # division
3+2*5              # note order of operations exists
(3+2)*5            # if you need, force operations using
                   # redundant parentheses
4^3                # raise to a power is ^
exp(4)             # e^4 = 54.598 (give or take)
log(2.742)         # natural log of 2.74
log10(1000)        # common log of 1000
pi                 # 3.14159... 
```

<h4 style="font-weight:bold;color:navy"> Store result (or anything) in a variable (object)</h4>

```R
# Assign a number in a variable
x <- 5
# This works exactly the same
x = 5
# Now just type the variable name to print out its content
5
#
# Let get more complicated
longvariablename = "This is the content of a variable with a longe name"
longvariablename
# Store the result of a calculation in the variable
y = 1+2+3+4+5+6+7+8
y
```

<h4 style="font-weight:bold;color:navy"> R data types </h4>

R has 5 different data types:
1) character
2) numeric (real or decimal)
3) integer
4) logical (true or false)
5) complex


<h4 style="font-weight:bold;color:navy"> R data structures </h4>

The data (of different types) can be combined to form data structures
1) vectors - list of multiple data of the same type
2) list - list of multiple data of various types
3) matrix 
4) data frame
5) factors

<h4 style="font-weight:bold;color:navy"> Beyound basic R - install R packages </h4>

Copy and paste the below R package installation codes to install the packages that we will be using in this workshop:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("ggplot2")  
BiocManager::install("ggpubr")
BiocManager::install("phyloseq")
BiocManager::install("Biostrings")
BiocManager::install("xlsx")

```

After installing these packages, we need to load them before we can use it, let's load them all in a single shot:

```R
library("dada2")
library("ggpubr")
library("ggplot2")
library("phyloseq")
library("Biostrings")
library("xlsx")
```

<br><br>

<h3 style="font-weight:bold;color:#008C23"> 4. Use of the R DADA2 package to process the 16S rRNA gene sequence reads</h3>

<h4 style="font-weight:bold;color:#008C23"> Step 1. First we need to create a working directory and  download a set of demo Illumina 16S rRNA sequences</h4>

```R
#create a directory to work in, so everything is contained in this folder
dir.create("resispart")
#set the working directory to resispart
setwd("resispart")
#double check the currrent working directory is resispart
getwd()
```
You should see something like this:

<img src="https://gyazo.com/b7474148496ee1058742822a21508595.png">

Next, do this to download and extract a set of demo files:

```R
#download from an FTP site
download.file("ftp://www.homd.org/pub/resistpart/RESISPART_Bioinformatics_Workshop_Demo.zip","temp.zip")
## then unzip it, exdir tells R to just extract to the current working directory "."
unzip("temp.zip",exdir=".")
#list the content of the unzipped files
list.files()
```
The file size is 1.1Gb so it may take some time to download, depending on the Internet connection speed.
<img src="https://i.gyazo.com/e6e4ebd1a04be63bd6cf1cdd3b222297.png">

<table bgcolor="#C8FBD3"><tr><td style="font-size:0.8em" width="100%">
If the above link doesn't work use this alternative link:<br>
http://www.homd.org/ftp/pub/resistpart/RESISPART_Bioinformatics_Workshop_Demo.zip<br>
  </td></tr></table>


<h4 style="font-weight:bold;color:#008C23"> Step 2. Find the sequence data</h4>

```R
dir.create("fastq", showWarnings = F)
path <- "fastq"
#list.files(path)
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.names 

```

<h4 style="font-weight:bold;color:#008C23"> Step 3. Quality plots</h4>

```R
plotQualityProfile(c(fnFs[1],fnRs[1]))
```

```R
dir.create("quality_plots", showWarnings = F)

quality_plots=list()
for (i in 1:length(fnRs)) {
quality_plots[[i]]=plotQualityProfile(c(fnFs[i],fnRs[i]))
}
quality_plots_all=ggarrange(plotlist=quality_plots,ncol=1,nrow=length(fnRs))
pdf("quality_plots/quality_plots_all.pdf", width=10, height=3*length(fnRs))
quality_plots_all
dev.off()

```

<h4 style="font-weight:bold;color:#008C23"> Step 4. Filter and trimm sequences</h4>

```R
dir.create("filtered", showWarnings = FALSE)
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#reads.filtered.start=Sys.time()
reads.filtered <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE)
#reads.filtered.finished=Sys.time()

```

<h4 style="font-weight:bold;color:#008C23"> Step 5. Learn sequence error rates to build error models for denoise purpose</h4>

```R
dir.create("err_plots", showWarnings = FALSE)

#learnErrors.start=Sys.time()
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#learnErrors.finished=Sys.time()

pdf("err_plots/errF.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf("err_plots/errR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

```

<h4 style="font-weight:bold;color:#008C23"> Step 6. De-replicate reads</h4>

```R
#derepFastq.start=Sys.time()
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#derepFastq.finished=Sys.time()

```

<h4 style="font-weight:bold;color:#008C23"> Step 7. Denoise sequences based on the error model to product amplicon sequence variants (ASVs)</h4>

```R
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#dada2.start=Sys.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dada2.finished=Sys.time()

```

<h4 style="font-weight:bold;color:#008C23"> Step 8. Merge the pair-end reads to single reads</h4>

```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
write.table(t(seqtab),file="seqtab.txt",sep="\t",quote=FALSE)
```

<h4 style="font-weight:bold;color:#008C23"> Step 9. Identify chimera and remove them</h4>

```R
#bimera.start=Sys.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#bimera.finished=Sys.time()
write.table(t(seqtab.nochim),file="seqtab.nochim",sep="\t",quote=FALSE)

```

<h4 style="font-weight:bold;color:#008C23"> Step 10. Make a summary table for the above processes</h4>

```R
library(xlsx)
getN <- function(x) sum(getUniques(x))
track <- cbind(reads.filtered, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track=rbind(track,colSums(track))
rownames(track)[nrow(track)]="Row Sum"
track=rbind(track,(track[nrow(track),]/track[nrow(track)])*100)
rownames(track)[nrow(track)]="Percentage"
write.xlsx(t(track),file="dada2_summary.xlsx",sheetName="dada2 summary")

```

<h4 style="font-weight:bold;color:#008C23"> Step 11. Assign ASVs with taxonomy</h4>

```R
#taxa.start=Sys.time()
taxa=assignTaxonomy(seqtab.nochim, "taxonomy/MOMDHOMDEXTGGNCBI_7S.fa.gz",tryRC=TRUE)
#taxa.finished=Sys.time()

taxa.print <- taxa
rownames(taxa.print) <- NULL
taxa.print
head(taxa.print)

```




<h3 style="font-weight:bold;color:brown"> 5. Use of the R Phyloseq package to study microbial diversity</h3>

<h4 style="font-weight:bold;color:brown"> Step 1. Load some necessary libraries</h4>

```R
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())    
```


<h4 style="font-weight:bold;color:brown"> Step 2. Convert the ASV from DADA2 into a Phyloseq object</h4>

```R
rownames(seqtab.nochim)
samples.out=rownames(seqtab.nochim)
meta=read.table(file="meta.txt",sep="\t",header=TRUE,colClasses="factor",row.names=1)
#meta
ps=phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),sample_data(meta),tax_table(taxa))
#ps
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

```


<h4 style="font-weight:bold;color:brown"> Step 3. Make some alpha-diversity plots</h4>

```R
alpha_plot=plot_richness(ps, x="Concentration", measures=c("Observed","Shannon", "Simpson"), title="Antibiotics Treatment") + geom_boxplot() + theme(plot.title = element_text(hjust = 0.5))
group_order=c("0 ug/ml","25 ug/ml","50 ug/ml","100 ug/ml")
alpha_plot$data$Concentration<-as.character(alpha_plot$data$Concentration)
alpha_plot$data$Concentration<-factor(alpha_plot$data$Concentration,levels=group_order)
alpha_plot
pdf(file="alpha_diversity_box_plot.pdf",width=8,height=6)
alpha_plot
dev.off()

```

<h4 style="font-weight:bold;color:brown"> Step 4. Make some beta-diversity plots</h4>

```R
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ps.prop.nmds.bray=ordinate(ps.prop, method="NMDS", distance="bray")
nmds1=plot_ordination(ps.prop,ps.prop.nmds.bray,type="samples",color="Concentration") + geom_point(size=2) + ggtitle("RESISPART DEMO NMDS Plot") + theme(plot.title = element_text(hjust = 0.5))
nmds1$data$Concentration<-as.character(nmds1$data$Concentration)
nmds1$data$Concentration<-factor(nmds1$data$Concentration,levels=group_order)
nmds1

pdf(file="beta_diversity_nmds_plot1.pdf",width=6,height=6)
nmds1
dev.off()

nmds2=plot_ordination(ps.prop,ps.prop.nmds.bray,type="samples",color="Concentration") + geom_point(size=2) + ggtitle("RESISPART DEMO NMDS") + theme(plot.title = element_text(hjust = 0.5)) + geom_polygon(aes(fill=Concentration))
nmds2$data$Concentration<-as.character(nmds2$data$Concentration)
nmds2$data$Concentration<-factor(nmds2$data$Concentration,levels=group_order)
nmds1

pdf(file="beta_diversity_nmds_plot2.pdf",width=6,height=6)
nmds1
dev.off()

```


<h4 style="font-weight:bold;color:brown"> Step 5. Plot genus and species level bar charts</h4>

```R
ps.genus=tax_glom(ps.prop,taxrank="Genus")
top10.genus <- names(sort(taxa_sums(ps.genus), decreasing=TRUE))[1:10]
ps.top10.genus <- prune_taxa(top10.genus,ps.genus)
plot_bar(ps.top10.genus,x="Concentration",fill="Genus")

ps.species=tax_glom(ps.prop,taxrank="Species")
top10.species <- names(sort(taxa_sums(ps.species), decreasing=TRUE))[1:10]
ps.top10.species <- prune_taxa(top10.species,ps.species)
plot_bar(ps.top10.species,x="Concentration",fill="Species")

```

<h4 style="font-weight:bold;color:brown"> Step 6. Identify differentially represented species</h4>

```R
#### Under construction

```

<h4 style="font-weight:bold;color:brown"> Step 7. Export tables for more downstream analysis using MicrobiomAnalyst</h4>

```R
library(tibble)
write.table(rownames_to_column(as.data.frame(t(otu_table(ps))),'#NAME'),file="otu_table.txt",sep="\t",quote=FALSE,row.names=F)
write.table(rownames_to_column(as.data.frame(meta),'#NAME'),file="meta_table.txt",sep="\t",quote=FALSE,row.names=F)
write.table(rownames_to_column(as.data.frame(tax_table(ps)),'#TAXONOMY'),file="tax_table.txt",sep="\t",quote=FALSE,row.names=F)

```


<h3 style="font-weight:bold;color:brown"> 6. Use of the MicrobiomeAnalyst online tools for statistical, visual and meta-analysis of microbiome data</h3>

Achal Dhariwal will lead this section to perform downstream analyses using the online tool MicrobiomeAnalyst that he and his colleages at McGill University developed. 

MicrobiomeAnalyst web site:
<https://www.microbiomeanalyst.ca/>

Online Tutorial for MicrobiomeAnalysis:
<https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/tutorials/MDP.pdf>

In case the site is not avaialble, Achal provided the R codes for the MicrobiomeAnalyst:
<ftp://www.homd.org/pub/resistpart/Downstream%20_analysis_resourcesMA.zip>


<h3 style="font-weight:bold;color:purple"> 7. Beyond microbiome sequence data – Meta-genomic and Meta-transcriptomic data</h3>

Under contruction ...
