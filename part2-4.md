---
permalink=/part2-4/
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
dir.create("fastq", showWarnings = F)
path <- "fastq"
#list.files(path)
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.names 

```

<a name="A3">
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

<a name="A4">
<h4 style="font-weight:bold;color:#008C23"> Step 4. Filter and trimm sequences</h4>

```R
dir.create("filtered", showWarnings = FALSE)
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#reads.filtered.start=Sys.time()
reads.filtered <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE)
#reads.filtered.finished=Sys.time()

```

<a name="A5">
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


