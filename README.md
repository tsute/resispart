<h2>Table of Content</h2>

[Instruction: How to use this workshop tutorial](part1)

[Workshop Tutorial](part2)

* [1. Introduction](part2/#A1)
* [2. Installation of the R software and packages for analyzing microbiome sequence data](part2/#A2)
* [3. Basic R operations](part2/A3)
* [4. Use of the R DADA2 package to process the 16S rRNA gene sequence reads](part2-4)
	* [Step 1. Create a working directory and download a set of demo Illumina 16S rRNA sequences](part2-4/#A1)
	* [Step 2. Find the sequence data](part2-4/#A2)
	* [Step 3. Quality plots]([part2-4/$A3)
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
	* [Step 3. Make some alpha-diversity plots]([part2-5/$A3)
	* [Step 4. Make some beta-diversity plots](part2-5/#A4)
	* [Step 5. Plot genus and species level bar charts](part2-5/#A5)
	* [Step 6. Identify differentially represented species](part2-5/#A6)
	* [Step 7. Export tables for more downstream analysis using MicrobiomAnalyst](part2-5/#A7)
* [6. Use of the MicrobiomeAnalyst online tools for statistical, visual and meta-analysis of microbiome data](part2-6)
* [7. Beyond microbiome sequence data – Meta-genomic and Meta-transcriptomic data](part2-7)  







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
