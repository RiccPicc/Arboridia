# ANALISI 16S CICALINE
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(DESeq2); packageVersion("DESeq2")
library(knitr)
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(DESeq2); packageVersion("DESeq2")
library("ape"); packageVersion("ape")
library(knitr)
library(stringr)
library("vegan")
library('ggpubr')
library('cowplot')
library('ampvis2')
library('ggvegan')
library('ggbiplot')
library("dichromat")
library("RColorBrewer")
library('UpSetR')
library(rstatix)
library("pheatmap") 
library(dplyr)
library("vsn")
library(writexl)
library('biomformat')
library("data.table")   # Also requires R.utils to read gz and bz2 files
library("phyloseq")
library("ALDEx2")
library("tidyverse")
library(readr)
library(dplyr)
source("/Users/riccardo/Library/CloudStorage/GoogleDrive-riccardo.piccinno@unipv.it/.shortcut-targets-by-id/1OLHPrObneyH2cOmEvB_ahqndO3QTORXC/H_halys/metaBUG/Sequenze_MiSeq/functions.r") # Mac path
source("G:/Il mio Drive/PhD.AES.Dsuz.Riccardo.Piccinno/H_halys/metaBUG/Sequenze_MiSeq/functions.r") # Windows path

setwd('/Users/riccardo/Library/CloudStorage/GoogleDrive-riccardo.piccinno@unipv.it/My Drive/Dottorato/Cicaline/16S') # Mac path
load('/Users/riccardo/Library/CloudStorage/GoogleDrive-riccardo.piccinno@unipv.it/My Drive/Dottorato/Cicaline/16S/dada2.Rdata') # Mac path

setwd('G:/Il mio Drive/PhD.AES.Dsuz.Riccardo.Piccinno/Cicaline/16S') # Windows path
load('G:/Il mio Drive/PhD.AES.Dsuz.Riccardo.Piccinno/Cicaline/16S/dada2.Rdata') # Windows path

###### dataset preparation #####
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "-"), `[`, 2)
population <- c(rep('Apulia',5),rep('Outgroup',5),rep('Crete',4),
                rep('Dalmatia',5),'Neg_K')
genus <- c(rep('Arboridia',5),rep('Empoasca',5),rep('Arboridia',9),'None')
samdf <- data.frame(Subject=subject, Population=population, Genus=genus, 
                    quant_reading=as.data.frame(track)$nonchim)
rownames(samdf) <- samples.out
# Generating a phyloseq object with all phylogenetic information on OTUs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
tree <- rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna, tree)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
df$is.neg <- df$Population == "Neg_K"
sample_data(ps)$is.neg <- sample_data(ps)$Population == "Neg_K"
ggplot(data=df, aes(x=Index, y=LibrarySize, color=is.neg)) + geom_point()
contamdf <- isContaminant(ps, method="either", conc="quant_reading", neg='is.neg', threshold=0.05)
table(contamdf$contaminant)
which(contamdf$contaminant)
# The default threshold for a contaminant is that it reaches a probability of 
# 0.1 in the statistical test being performed. In the prevalence test there is 
# a special value worth knowing, threshold=0.5, that will identify as contaminants 
# all sequences thare are more prevalent in negative controls than in positive samples.
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Population == "Neg_K", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Population != "Neg_K", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ps.noncontam <- prune_taxa(!contamdf$contaminant, ps)
ps.noncontam <- prune_samples(sample_data(ps.noncontam)$Population != "Neg_K", ps.noncontam)
ps <- ps.noncontam
ncolors <- length(levels(as.factor(sample_data(ps)$Population)))
primary_color_list <- brewer.pal(ncolors, 'Dark2') 
sample_data(ps)$Population <- factor(sample_data(ps)$Population, levels=c('Apulia','Crete','Dalmatia','Outgroup'))
ps.prop <- transform_sample_counts(ps, function(x) x/sum(x)*100)

###### save OTU table and sequences ######
OTU1 <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
OTUdf <- data.frame(t(OTU1)) # taxa rows, samples columns
write_xlsx(cbind(" "=rownames(OTUdf), OTUdf),'OTU_table.xlsx')
OTUbiom <- make_biom(OTUdf)
write_biom(OTUbiom, 'OTU_table.biom')
ps %>%
  refseq() %>%
  Biostrings::writeXStringSet("./asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

refseq.df <- as.data.frame(refseq(ps))

dev.off()

ps.genus <- tax_glom(ps, taxrank = "Genus")

wol_Apulia <- subset_taxa(ps, Genus=="Wolbachia")
wol_asvs <- colnames(otu_table(wol_Apulia))
wol_seqs <- as.list(refseq.df[wol_asvs,])
wol_names <- c()
for (i in 1:length(wol_seqs)){wol_names <- c(wol_names, paste("Wolbachia_Arboridia_",i,sep=""))}
names(wol_seqs) <- wol_names
library(seqinr)
write.fasta(sequences=as.list(wol_seqs), names=wol_names, file.out="wolbachia_arboridia_seqs.fasta", open = "w", nbchar = 60, as.string = FALSE)


ps.prop.renamed <- ps.prop
cabd <- sort(taxa_sums(ps.prop.renamed)/nsamples(ps.prop.renamed),
             decreasing = T)
# Lets remove some low abundance phyla
renamed.clade <- names(which(cabd < 0.1))
tax_table(ps.prop.renamed)[renamed.clade, 'Phylum'] <- '< 0.1% Abd'

png(paste('./pop_phylobars_full_','Phylum','.png',sep=''), width=9.5,units="in", height=5, res=1200)
print(my_plot_bar(ps.prop.renamed, fill='Phylum', x="Subject") + facet_wrap(~Population, scales="free_x", nrow=1))
dev.off()



# Agglomerate taxa at genus level
ps <- prune_samples(sample_data(ps)$Population != "Outgroup", ps)
ps.genus <- prune_samples(sample_data(ps.genus)$Population != "Outgroup", ps.genus)
sample_data(ps.genus)$Subject <- c(1,2,3,4,5,1,2,3,4,1,2,3,4,5)

# Top N taxa
N <- 20
top <- names(sort(taxa_sums(ps.genus), decreasing = TRUE))[1:N]

# Calculate relative abundance
ps.genus.prop <- transform_sample_counts(ps.genus, function(x) x / sum(x) )

# Subset object to top N taxa
ps.genus.prop.top <- prune_taxa(top, ps.genus.prop)

df_top <- ps.genus.prop.top %>% psmelt()
df <- as.data.frame(cbind(df_top$Population, df_top$Subject, df_top$OTU, df_top$Genus))
colnames(df) <- c("Population","Specimens","ASV name","Genus")
df <- df[order(df[,1], df[,2]),]
write_xlsx(df,"abundant_genus.xlsx")


ps.prop.renamed <- ps.prop
cabd <- sort(taxa_sums(ps.prop.renamed)/nsamples(ps.prop.renamed),
             decreasing = T)
# Lets remove some low abundance phyla
renamed.clade <- names(which(cabd < 1.5))
tax_table(ps.prop.renamed)[renamed.clade, 'Genus'] <- '< 1.5% Abd'

png(paste('./pop_phylobars_full_','Genus','.png',sep=''), width=9.5,units="in", height=5, res=1200)
print(my_plot_bar(ps.prop.renamed, fill='Genus', x="Subject") + facet_wrap(~Population, scales="free_x", nrow=1)+theme(text = element_text(size=18)))
dev.off()

all.div.plots_gen <- alpha_diversity_pop(ps, "Genus")
all.div.plots_phy <- alpha_diversity_pop(ps, 'Phylum')
beta_phy_pop_plot <- beta_diversity_pop_general(ps,clade="Phylum",loc="",prev_th=0.1,type_dist='bray')
beta_gen_pop_plot <- beta_diversity_pop_general(ps,clade="Genus",loc="",prev_th=0.1,type_dist='bray')


# Fig. S1 - alpha beta div plots

top_row <- plot_grid(all.div.plots_gen[[2]]+ theme(text = element_text(size = 16)),
                     all.div.plots_gen[[1]]+ theme(text = element_text(size = 16)),
                     beta_gen_pop_plot + geom_point(size=1) + scale_shape_manual(values = c(1,2,4), breaks=c('Apulia', 'Crete', 'Dalmatia')) + theme(text = element_text(size = 16)), 
                     labels = c('A', 'B', 'C'), label_size = 12, rel_widths = c(1.7,1.4,2), ncol=3, nrow=1)
bot_row <- plot_grid(all.div.plots_phy[[2]]+ theme(text = element_text(size = 16)),
                     all.div.plots_phy[[1]]+ theme(text = element_text(size = 16)),
                     beta_phy_pop_plot + geom_point(size=1) + scale_shape_manual(values = c(1,2,4), breaks=c('Apulia', 'Crete', 'Dalmatia')) + theme(text = element_text(size = 16)), 
                     labels = c('D', 'E', 'F'), label_size = 12, rel_widths = c(1.7,1.4,2), ncol=3, nrow=1)
png('FigureS1.png', width=15,units="in", height=13, res=300)
plot_grid(top_row,bot_row,nrow=2, ncol=1)
dev.off()

my_plot_bar <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                         facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack") + xlab('') + labs(x='')
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 1))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


ps.prop <- prune_samples(sample_data(ps.prop)$Population != "Outgroup", ps.prop)
sample_data(ps.prop)$Subject <- c(1,2,3,4,5,1,2,3,4,1,2,3,4,5)
ps.prop.renamed <- ps.prop
cabd <- sort(taxa_sums(ps.prop.renamed)/nsamples(ps.prop.renamed),
             decreasing = T)
# Lets remove some low abundance genera
renamed.clade <- names(which(cabd < 1.5))
tax_table(ps.prop.renamed)[renamed.clade, 'Genus'] <- '< 1.5% Abd'
legend_labels <- levels(factor(as.data.frame(tax_table(ps.prop.renamed))$Genus))
legend_labels[legend_labels != '< 1.5% Abd'] %>% 
  paste('italic(',.,')', sep = "") %>%
  parse(text = .) -> eNames
legend_labels <- c('< 1.5% Abd', eNames, "NA")

png(paste('./Fig_8.png',sep=''), width=9.5,units="in", height=5, res=1200)
print(my_plot_bar(ps.prop.renamed, fill='Genus', x="Subject") + facet_wrap(~Population, scales="free_x", nrow=1)) +
  theme(text = element_text(size=18)) +
  scale_fill_discrete(labels=legend_labels)
dev.off()
