<h2>Table of Content</h2>

[Instruction: How to use this workshop tutorial](part1)

[Workshop Tutorial](part2)

* [1. Introduction](part2/#A1)
* [2. Installation of the R software and packages for analyzing microbiome sequence data](part2/#A2)
* [3. Basic R operations](part2/A3)
* [4. Use of the R DADA2 package to process the 16S rRNA gene sequence reads](part2-4)
	* [Step 1. Create a working directory and download a set of demo Illumina 16S rRNA sequences](part2-4/#A1)
	* [Step 2. Find the sequence data](part2-4/#A2)
	* [Step 3. Quality plots]([part2-4/#A3)
	* [Step 4. Filter and trimm sequences](part2-4/#A4)
	* [Step 5. Learn sequence error rates to build error models for denoise purpose](part2/#A5)
	* [Step 6. De-replicate reads](part2-4/#A6)
	* [Step 7. Denoise sequences based on the error model to product amplicon sequence variants (ASVs)](part2-4/#A7)
	* [Step 8. Merge the pair-end reads to single reads](part2-4/#A8)
	* [Step 9. Identify chimera and remove them](part2-4/#A9)
	* [Step 10. Make a summary table for the above processes](part2-4/#A10)
	* [Step 11. Assign ASVs with taxonomy](part2-4/#A11)
* [5. Use of the R Phyloseq package to study microbial diversity](part2-5)
	* [Step 1. Load some necessary libraries](part2-5/#A1)
	* [Step 2. Convert the ASV from DADA2 into a Phyloseq object](part2-5/#A2)
	* [Step 3. Make some alpha-diversity plots]([part2-5/#A3)
	* [Step 4. Make some beta-diversity plots](part2-5/#A4)
	* [Step 5. Plot genus and species level bar charts](part2-5/#A5)
	* [Step 6. Identify differentially represented species](part2-5/#A6)
	* [Step 7. Export tables for more downstream analysis using MicrobiomAnalyst](part2-5/#A7)
* [6. Use of the MicrobiomeAnalyst online tools for statistical, visual and meta-analysis of microbiome data](part2-6)
* [7. Beyond microbiome sequence data – Meta-genomic and Meta-transcriptomic data](part2-7)  








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
