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

<a name="A3">
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

<a name="A4">
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
