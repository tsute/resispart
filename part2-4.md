---
permalink: /part2-4/
---
<h3 style="font-weight:bold;color:#008C23"> 4. Use of the R DADA2 package to process the 16S rRNA gene sequence reads</h3>
<a name="A1">
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

<a name="A2">
<h4 style="font-weight:bold;color:#008C23"> Step 2. Find the sequence data</h4>

```R
path <- "fastq"
#list.files(path)
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.names 

```

<a name="A3">
<h4 style="font-weight:bold;color:#008C23"> Step 3. Quality plots</h4>

Make a plot for the read pairs (R1 and R2) of the first sample:
```R
plotQualityProfile(c(fnFs[1],fnRs[1]))
```
Quality plot output:

<img src="https://i.gyazo.com/f78c3b831e38438c7f7b6b1ffb4ebdd0.png">

Now we make a loop to plot all the read pairs in a single shot:
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

You will see some warnings but we can safely ignore them:

<img src="https://i.gyazo.com/4a8a27806fec46380ec3cca81f6014f6.png">

This will generate a long PDF file with quality plots for all 12 samples. Go to windows desktop folder "resispart" and the PDF plot is called "quality_plots_all.pdf" in the "quality_plots" folder.

Open the pdf file you will see 12 pairs of quality plots in one pdf file:

<img src="https://gyazo.com/8a563c82c605253f3b131aa95fd616ab.png" border="2">

The purpose of examining these quality plots is to help us pick a trimming length for the reads. Usually the sequence quality degrades more rapidly toward the end of the read, this is especially so for R2.

<a name="A4">
<h4 style="font-weight:bold;color:#008C23"> Step 4. Filter and trimm sequences</h4>

Based on the above sequence plots, we observe that quality score starts to drop rapidly after base 280 for R1 and 220 for R2 so next we will trim the bases off the reads after these two positions:

```R
dir.create("filtered", showWarnings = FALSE)
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#reads.filtered.start=Sys.time()
reads.filtered <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE)
#reads.filtered.finished=Sys.time()

```

<table class="altbg"><tr><td style="font-size:0.8em" width="100%">
This may take a while depending on your computer's speed. If it is taking too long, we can cancel the "filterAndTrim" command by pressing the Escape key and then the Enter to to cancel the filtering process. We will used a pre-filtered result in the "filtered_backup" folder. Go to Windows explorer and rename this folder to just "filtered" so we can continue to the next step. 
</td></tr></table>


<a name="A5">
<h4 style="font-weight:bold;color:#008C23"> Step 5. Learn sequence error rates to build error models for denoise purpose</h4>

The next step is the first phase of the DADA2 amplicon denoising procedure. DADA2 is the second version of the DADA (The Divisive Amplicon Denoising Algorithm)([Rosen et al, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3563472/)). DADA was developed to model the 454-sequenced amplicon data, whereas DADA2 was designed specifically for the Illumina sequence reads ([Callahan et al, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)). In this step we use the quality-filtred reads for building the error models required by the DADA2 denoising algorithm.

```R
dir.create("err_plots", showWarnings = FALSE)

#learnErrors.start=Sys.time()
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#learnErrors.finished=Sys.time()

plotErrors(errF, nominalQ=TRUE)

pdf("err_plots/errF.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

plotErrors(errR, nominalQ=TRUE)

pdf("err_plots/errR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

```

The plots look like this:

<img src="https://i.gyazo.com/10a7c328b4a1568f67aa58f52e165740.png">

There are a total of 12 possible nucleotide transition (e.g., A→C, A→G, …). The plot shows the pair-wise transition of A, T, C and G (including self). The dots indicate the observed error rates for each consensus quality score. The black lines are the estimated error rates after convergence of the machine-learning algorithm. The red lines are the error rates expected under the nominal definition of the Q-score. Ideally we want to see a good fit between the black dots and the black line. As you can see the higher the quality scores the lower the error frequencies. 

<table class="altbg"><tr><td style="font-size:0.8em" width="100%">
This step is computationally intensive, we used the default setting for the learnErrors function and used only first 100M bases to learn error model. If you don't see a roughly good fit between dots and the lines, you can try increasing the "nbases" parameter to see if the fit improves.
</td></tr></table>


<a name="A6">
<h4 style="font-weight:bold;color:#008C23"> Step 6. De-replicate reads</h4>

```R
#derepFastq.start=Sys.time()
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#derepFastq.finished=Sys.time()

```

<a name="A7">
<h4 style="font-weight:bold;color:#008C23"> Step 7. Denoise sequences based on the error model to product amplicon sequence variants (ASVs)</h4>

```R
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#dada2.start=Sys.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dada2.finished=Sys.time()

```

<a name="A8">
<h4 style="font-weight:bold;color:#008C23"> Step 8. Merge the pair-end reads to single reads</h4>

```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
write.table(t(seqtab),file="seqtab.txt",sep="\t",quote=FALSE)
```

<a name="A9">
<h4 style="font-weight:bold;color:#008C23"> Step 9. Identify chimera and remove them</h4>

```R
#bimera.start=Sys.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#bimera.finished=Sys.time()
write.table(t(seqtab.nochim),file="seqtab.nochim",sep="\t",quote=FALSE)

```

<a name="A10">
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

<a name="A11">
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

## [Next ▶](/resispart/part2-5)
