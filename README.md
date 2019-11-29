<h2>Table of Content</h2>

[Instruction: How to use this workshop tutorial](part1)

[Workshop Tutorial](part2)

* [1. Introduction](part2/#1)
* [2. Installation of the R software and packages for analyzing microbiome sequence data](part2/#2)
* [3. Basic R operations](part2/#3)
* [4. Use of the R DADA2 package to process the 16S rRNA gene sequence reads](part2-4/#4)
	* [Step 1. Create a working directory and download a set of demo Illumina 16S rRNA sequences](part2-4/#4-1)
	* [Step 2. Find the sequence data](part2-4/#4-2)
	* [Step 3. Quality plots]([part2-4/$4-3)
	* [Step 4. Filter and trimm sequences](part2-4/#4-4)
	* [Step 5. Learn sequence error rates to build error models for denoise purpose](part2/#4-5)
	* [Step 6. De-replicate reads](part2-4/$4-6)
	* [Step 7. Denoise sequences based on the error model to product amplicon sequence variants (ASVs)](part2-4/#4-7)
	* [Step 8. Merge the pair-end reads to single reads](part2-4/#4-8)
	* [Step 9. Identify chimera and remove them](part2-4/#4-9)
	* [Step 10. Make a summary table for the above processes](part2-4/#4-10)
	* [Step 11. Assign ASVs with taxonomy](part2-4/#4-11)
* [5. Use of the R Phyloseq package to study microbial diversity](part2-5)
	* [Step 1. Load some necessary libraries](part2-5/#5-1)
	* [Step 2. Convert the ASV from DADA2 into a Phyloseq object](part2-5/#5-2)
	* [Step 3. Make some alpha-diversity plots]([part2-5/$5-3)
	* [Step 4. Make some beta-diversity plots](part2-5/#5-4)
	* [Step 5. Plot genus and species level bar charts](part2-5/#5-5)
	* [Step 6. Identify differentially represented species](part2-5/$5-6)
	* [Step 7. Export tables for more downstream analysis using MicrobiomAnalyst](part2-5/#4-7)
* [6. Use of the MicrobiomeAnalyst online tools for statistical, visual and meta-analysis of microbiome data](part2-6)
* [7. Beyond microbiome sequence data – Meta-genomic and Meta-transcriptomic data](part2-7)  





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
