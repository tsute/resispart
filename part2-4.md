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
download.file("https://github.com/tsute/resispart/releases/download/demo_data/RESISPART_Bioinformatics_Workshop_Demo.zip","demo.zip")

## then unzip it, exdir tells R to just extract to the current working directory "."
unzip("demo.zip",exdir=".")

#list the content of the unzipped files
list.files()

```
The file size is 1.1Gb so it may take some time to download, depending on the Internet connection speed.
<img src="https://i.gyazo.com/e6e4ebd1a04be63bd6cf1cdd3b222297.png">

<table bgcolor="#C8FBD3"><tr><td style="font-size:0.8em" width="100%">
If the above link doesn't work try one of these alternative download links:<br><br>
<a href="https://bioinformatics.forsyth.org/ftp/pub/resistpart/RESISPART_Bioinformatics_Workshop_Demo.zip">Forsyth Bioinformatics Download Link</a><br>
  <a href="https://drive.google.com/open?id=1A5fHVqJ2Nfloxvs-Ej8E-F2MhawvVqxe">Google Drive</a><br>
  <a href="https://1drv.ms/u/s!Amu_vVYXMX9XicpMCEljLG-u59QPDw?e=vING0V">One Drive</a><br>
  <a href="https://www.dropbox.com/s/qm0gull017umtot/RESISPART_Bioinformatics_Workshop_Demo.zip?dl=0">Dropbox</a><br>
  <br>
  Click on one of the above links and follow the direction to download the file <b>RESISPART_Bioinformatics_Workshop_Demo.zip</b><br>
  Find the file on your computer, double-click to unzip it and then move the entire content to the "resispart" working directory.

  </td></tr></table>

<a name="A2">
<h4 style="font-weight:bold;color:#008C23"> Step 2. Find the sequence data</h4>

```R
path <- "fastq"
#list.files(path)
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE)) # store forward (R1) sequence files in this object
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE)) # store reverse (R2) sequence files in this object
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# type the object name to see what sample names we have
sample.names 

```

<img src="https://i.gyazo.com/af5b6e496e6b3e79476a26b02cf2b64d.png">

In this workshop we have all together 3,383,251 x 2 sequences in 12 samples (each sample has two files R1 and R2). These are called pair-end reads. The Illumina sequencer can sequence both ends of a DNA amplicon at the same time, expanding the sequence length and span. If the two reads overlap enough we can also merge them together into a single sequence from each pair of the reads. We will process them separately first for quality filtering and later merge them together for downstream analysis.


<a name="A3">
<h4 style="font-weight:bold;color:#008C23"> Step 3. Quality plots</h4>

Make a plot for the read pairs (R1 and R2) of the first sample:

```R
plotQualityProfile(c(fnFs[1],fnRs[1]))


```
Quality plot output:

<img src="https://i.gyazo.com/f78c3b831e38438c7f7b6b1ffb4ebdd0.png">

<table  class="notebg"><tr><td style="font-size:0.8em;">
In gray-scale is a heat map of the frequency of each quality score at each base position. The median quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position
</td></tr></table>

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
reads.filtered.start=Sys.time()
reads.filtered <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE)
reads.filtered.finished=Sys.time()


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
  
Dereplication is to remove the redundant/repeated sequences from a set of sequences and keep only one copy of all the unique seuqences. Dereplication serves two purposes:
1. Reduce computaton time by emiminating repeated sequences.
2. The abudnace (copy number) infomration will be used in the partitioning step to form amplicon sequence variants.

Perform the following commands to dereplicate the reads:

```R
#derepFastq.start=Sys.time()
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#derepFastq.finished=Sys.time()


```

The result should be self-explanatory:

<img src="https://i.gyazo.com/b167b45e3120d746be22b5f750decb97.png">)



<a name="A7">
<h4 style="font-weight:bold;color:#008C23"> Step 7. Denoise sequences based on the error model to product amplicon sequence variants (ASVs)</h4>

This step is the central part of the DADA2 pipeline. In a very simplified version of explanation, DADA2 uses a consensus quality profile with the abundance information to partition the reads into amplicaon sequence variants (ASVs). In a sample, if a sequence is real and abundant, there should be many identical copies of the same sequences. The more copies of the same sequences, the more likely it is to be true sequences. DADA2 partition the reads based on most abundant unique reads, and then forms ASVs using the consensus quality information associated with the reads. 

Copy and paste the following R codes to perform the denoising step

```R
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#dada2.start=Sys.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dada2.finished=Sys.time()

```

<img src="https://i.gyazo.com/8302531cf335be4cdb4f5a303cc891b8.png">



<a name="A8">
<h4 style="font-weight:bold;color:#008C23"> Step 8. Merge the pair-end reads to single reads</h4>
  
Next we merge the forward (R1) and reverse (R2) reads together to generate the full denoised sequences. This is done by aligning the denoised R1 reads with the reverse-complement of the corresponding R2 reads. 

Perform the following command to merge the pair-end reads, and generate a liste of the merged ASV sequences.

```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
write.table(t(seqtab),file="seqtab.txt",sep="\t",quote=FALSE)
```

<img src="https://gyazo.com/fd3ccb0ec292173bac2d8b50eb104130.png">


<a name="A9">
<h4 style="font-weight:bold;color:#008C23"> Step 9. Identify chimera and remove them</h4>

The PCR step prior to the sequencing often generates chimera - a hybrid artifial sequences dervied from two different parent sequeces. The ASVs may still have these hybrid sequences and should be removed. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences. It is easier and more accurate to identify chimera from the denoised ASVs than the noisy raw sequences.

```R
#bimera.start=Sys.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#bimera.finished=Sys.time()
write.table(t(seqtab.nochim),file="seqtab.nochim",sep="\t",quote=FALSE)
write_xlsx(data.frame(seqtab.nochim),"seqtab.nochim.xlsx")

```

modified R


```R
#derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)
load("errF.RData")
load("errR.RData")
names(derepFs) <- sample.names
names(derepRs) <- sample.names
load("dadaFs.RData")
load("dadaRs.RData")
load("mergers.RData")
load("seqtab.nochim.RData")

write.table(t(seqtab.nochim),file="seqtab.nochim",sep="\t",quote=FALSE)

write_xlsx(data.frame(seqtab.nochim),"seqtab.nochim.xlsx")


```

<a name="A10">
<h4 style="font-weight:bold;color:#008C23"> Step 10. Make a summary table for the above processes</h4>
  
Use the following R codes to generate a summary of all the processes and keep track of number and percentage of reads after each step.

```R
load("reads.filtered.RData")
getN <- function(x) sum(getUniques(x))
track <- cbind(reads.filtered, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track=rbind(track,colSums(track))
rownames(track)[nrow(track)]="Row Sum"
track=rbind(track,(track[nrow(track),]/track[nrow(track)])*100)
rownames(track)[nrow(track)]="Percentage"
write_xlsx(data.frame(t(track)),"dada2_summary.xlsx")
 
```

<img src="https://i.gyazo.com/c5ad59bbc9aa68b739dc7ea64cd0e908.png">


<a name="A11">
<h4 style="font-weight:bold;color:#008C23"> Step 11. Assign ASVs with taxonomy</h4>

The ASVs generated by DADA2 are the sequence variants and do not provide the taxonomy information. To find out which bacteria these ASVs may have come from, we need to compare these ASVS to a set of reference 16S rRNA sequences with taxonomy information. If an ASV is most similar to a reference sequence, the ASV can then be assigned to the taxonomy that is associated with the matching reference. 

There are two general methods for finding the closest reference sequences (for the ASVs):
1. *Alignment:* The closest reference sequence is identified based on highest sequence percent similarity calculated between aligned ASV and reference sequences.
2. *Non-alignment:* the cloest reference sequence is identified based on the highest probability that the K-mers (e.g., all 8-mers in the ASVs) can be found in the reference. 

The alignment based method provides percent similarity information thus once the cloest reference sequence is identified, we also know how close the ASV is to the reference, in term of percent sequence identity. This information is useful to determine whether the ASV is  the same species (i.e., if 100% identical) of the reference. However the alignment is slow and can take a lot more computation time, especially if we are talking about millions of NGS reads. Luckily, the denoising process of DADA2 often reduces the number of reads from millions to hundreds, and thus drastically reduced the amount of time needed to calculate the percent similarity based on aligned sequences.

The K-mer based non-alignment method, on the other hand does not consider the position of the K-mer in the sequence and only calculate the (highest) possible probility for the vocabuary of K-mer words in an ASV to belong to a reference. Eventhough this method requires less computation time, the time diferennce between the alignment and non-alignment based methods is not significant when dealing with only hundreds of ASVs, rathaer than thousands of OTUs or millions of raw read data.

In this tutorial we will use the non-alignment taxonomy assignment method provided by the DADA2 package. This default k-mer based assignment is based on the naive Bayesian classifier method developed by the Ribosomal Database Project (RDP) ([Wang et al, 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950982/). The fucntion is called "assignTaxonomy". 

The assignment will need a set of 16S rRNA gene reference sequences that is well curated with taxonomy information. There are several choices of the reference sequences. The [DATA2 web site](https://benjjneb.github.io/dada2/training.html) provides several choices:

1. Silva version 132, Silva version 128, Silva version 123 (Silva dual-license)
2. RDP trainset 16, RDP trainset 14
3. GreenGenes version 13.8
4. UNITE (use the General Fasta releases)

For this tutorial we will be using a set of reference sequences maintained by the HOMD project, this set of reference sequences includes:

1. HMT RefSeq V15.1: 998sequences
2. HOMD RefSeq Extended V1.11: 151 sequences
3. GreenGeneGold V1: 2,623 sequences
4. NCBI 16S rRNA Reference: 18,044sequences
5. Mouse RefSeq V0.1: 82 sequences

All togethere there are 21,898 unique sequences representing a total of 14,651 prokarotic species. In addition the Genus names were added to the species level name so when the comparison is done at species level, the result will show the full scientific names for the species (i.e., Genus + Species).

This set of reference seuqences is included in the demo dataset in the taxonomy folder. 

Copy and paste below R codes to perform the taxonomy assignement task:

```R
#taxa.start=Sys.time()
taxa=assignTaxonomy(seqtab.nochim, "taxonomy/MOMDHOMDEXTGGNCBI_7S.fa.gz",tryRC=TRUE)
#taxa.finished=Sys.time()

taxa.print <- taxa
rownames(taxa.print) <- NULL
taxa.print
head(taxa.print)

```

Beginning of the taxa assignment output:

<img src="https://i.gyazo.com/d3b925b82d10fd020ef54e7cf1ab61db.png">

Once all the ASVs are all assgned with the potential taxonomy we are ready for many possible diversity analyses.

## [◀ Previous](/resispart/part2) ... [Next ▶](/resispart/part2-5)
