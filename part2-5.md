---
permalink: /part2-5/
---
<h3 style="font-weight:bold;color:brown"> 5. Use of the R Phyloseq package to study microbial diversity</h3>
<a name="A1">
<h4 style="font-weight:bold;color:brown"> Step 1. Load some necessary libraries</h4>

```R
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())    
```

<a name="A2">
<h4 style="font-weight:bold;color:brown"> Step 2. Convert the ASV from DADA2 into a Phyloseq object</h4>
  
If you look at the Excel file we generated at the end of the DADA2 step, we will see that the row names are the sample names and the column names are the long ASV sequences.

<img src="https://gyazo.com/d55451dfe4aac77cf72ae5282c1857c7.png">

It is very clumsy to use the long sequences as the names of the ASVs. The following R codes will help to rename the ASV using the ASV as the prefix and a sequential number for the IDs. Also we do not delete the ASV sequences in case we need them later. Hence the sequences are stroed in an object called "dna" using a function "DNAStringSet" provided by the Biostrings.

```R
# take a look at the row names of the seqtab.nochim object
rownames(seqtab.nochim)
# assign sample names to samples.out
samples.out=rownames(seqtab.nochim)
# Now we read in from a provided meta inforamtion file that tells us which sample is what treatment
meta=read.table(file="meta.txt",sep="\t",header=TRUE,colClasses="factor",row.names=1)
#Now we have 3 ingredient to build the phyloseq object, 1) ASV object, 2) sample information and 3) ASV taxonomy assignment
ps=phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),sample_data(meta),tax_table(taxa))
# Store the ASV DNA sequences to the object dna
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
# We can even store the dna seuqences in the phyloseq object ps, although they will not be used in this tutorial
ps <- merge_phyloseq(ps, dna)
# Here we rename the ASV in the format of ASV1, ASV2, ... ASV100 ... etc
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
 
```

Now the phyloseq object is ready for many downstream diversity analysis tools that are provided by the Phyloseq pacakge. We will practice a few of the basic ones.

<a name="A3">
<h4 style="font-weight:bold;color:brown"> Step 3. Make some alpha-diversity plots</h4>

Alpha diversity measures the number of taxonomy units (Genus, Species, OTUs, or ASVs) and their abundance within each sample. There are many different ways to measure/estimate the alpha diversity of the sample. The most common estimators for alpha diversity include 
1. Observed species - only count the total number of the species in the sample, does not care about the abundance)
2. Shannon estimator - considers the number of species that are more abundant, disregarding the very low abundant species)
3. Simpson estimator - considers the most abundant species, disregarding the low abudant species)

Depending on the question you ask, one may prefer one estimator over the other. Here we practice to calculate and plot the alpha diversity in our demo samples using these three most commonly used diversity measurements:


```R
# use the plot_richness function provided by Phyloseq package
alpha_plot=plot_richness(ps, x="Concentration", measures=c("Observed","Shannon", "Simpson"), title="Antibiotics Treatment") + geom_boxplot() + theme(plot.title = element_text(hjust = 0.5))
# make sure the plots are shown in the order we desire
group_order=c("0 ug/ml","25 ug/ml","50 ug/ml","100 ug/ml")
alpha_plot$data$Concentration<-as.character(alpha_plot$data$Concentration)
alpha_plot$data$Concentration<-factor(alpha_plot$data$Concentration,levels=group_order)
# now generate the plot
alpha_plot
# Also save a copy of the plot in the PDF format
pdf(file="alpha_diversity_box_plot.pdf",width=8,height=6)
alpha_plot
dev.off()

```

OUtput:

<img src="https://i.gyazo.com/12058110d020e3a507c92fb767148a7d.png">


<a name="A4">
<h4 style="font-weight:bold;color:brown"> Step 4. Make some beta-diversity plots</h4>

Beta diversity measures the similarity/dissimilarity between samples or groups of samples. Common ways to measure the similarity include:

1. Bray–Curtis dissimilarity
- based on abundance or read count data 
- differences in microbial abundances between two samples (e.g., at species level) 
    values are from 0 to 1
    0 means both samples share the same species at exact the same abundances
    1 means both samples have complete different species abundances

2. Jaccard distance
- based on presence or absence of species (does not include abundance information) 
- different in microbial composition between two samples
    0 means both samples share exact the same species
    1 means both samples have no species in common

3. UniFrac
- sequence distances (phylogenetic tree) 
- based on the fraction of branch length that is shared between two samples or unique to one or the other sample
unweighted UniFrac: purely based on sequence distances (does not include abundance information)
weighted UniFrac: branch lengths are weighted by relative abundances (includes both sequence and abundance information)


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

<a name="A5">
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

<a name="A6">
<h4 style="font-weight:bold;color:brown"> Step 6. Identify differentially represented species</h4>

```R
#### Under construction

```

<a name="A7">
<h4 style="font-weight:bold;color:brown"> Step 7. Export tables for more downstream analysis using MicrobiomAnalyst</h4>

```R
library(tibble)
write.table(rownames_to_column(as.data.frame(t(otu_table(ps))),'#NAME'),file="otu_table.txt",sep="\t",quote=FALSE,row.names=F)
write.table(rownames_to_column(as.data.frame(meta),'#NAME'),file="meta_table.txt",sep="\t",quote=FALSE,row.names=F)
write.table(rownames_to_column(as.data.frame(tax_table(ps)),'#TAXONOMY'),file="tax_table.txt",sep="\t",quote=FALSE,row.names=F)

```

## [Next ▶](/resispart/part2-6-7)
