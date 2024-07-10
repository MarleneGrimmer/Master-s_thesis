library(readr)
library(dplyr)
library(stringr)
library(ggnewscale)
library(ggplot2)

# choosing an outgroup based on ANI----
ANIb <- read_table("~/Uni/11.masters_thesis/Phylogeny/ANIb.csv", col_names = FALSE)
ANIb[1,] <- c(NA, ANIb[1,1:7]) #move names in first row one over

# replace genome IDs with strep names
ANIb[ANIb == "GCA_019048645.1_ASM1904864v1_genomic.fna"] <- "strep_mutans"
ANIb[ANIb == "GCA_000423725.1_ASM42372v1_genomic.fna"] <- "strep_devriesei"
ANIb[ANIb == "GCA_900475025.1_41906_G01_genomic.fna"] <- "strep_ferus"
ANIb[ANIb == "GCA_900459485.1_44738_A02_genomic.fna"] <- "strep_macacae"
ANIb[ANIb == "GCA_001431045.1_ASM143104v1_genomic.fna"] <- "strep_orisasini"
ANIb[ANIb == "GCA_000286075.1_ASM28607v1_genomic.fna"] <- "strep_ratti"
ANIb[ANIb == "GCA_002355215.1_ASM235521v1_genomic.fna"] <- "strep_troglodytae"

# make data frame
ANIb <- as.data.frame(ANIb)
rownames(ANIb)[2:8] <- ANIb$X1[2:8]
colnames(ANIb)[2:8] <- ANIb$X1[2:8]
ANIb$X1 <- NULL
ANIb <- ANIb[2:8,]

# phylogeny----
library(ape)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(panstripe)

# reading in gene presence/absence matrices
gene_pa_outgroup <- read_rtab("~/Uni/11.masters_thesis/Phylogeny/panaroo_outgroup/gene_presence_absence.Rtab")
gene_pa_nooutgroup <- read_rtab("~/Uni/11.masters_thesis/Phylogeny/panaroo_nooutgroup/gene_presence_absence.Rtab")

# tree with outgroup:
tree_outgroup <- read.newick("~/Uni/11.masters_thesis/iqtree/with_outgroup/iqtree.treefile")
rerooted_outgroup <- root(tree_outgroup, outgroup = "PROKKA_02152024", resolve.root = TRUE)
rerooted_outgroup_droptip <- drop.tip(rerooted_outgroup, "PROKKA_02152024")

# tree without outgroup
tree_nooutgroup <- read.newick("~/Uni/11.masters_thesis/iqtree/without_outgroup/nooutgroup.treefile")
midpointroot_nooutgroup <- midpoint_root(tree_nooutgroup)

# adding metadata to tree----
load("~/Uni/11.masters_thesis/metadata/metadata.Rds")
dataset <- read_delim("~/Uni/11.masters_thesis/ncbi_dataset.tsv", 
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dataset <- dataset[which(grepl("GCA_", dataset$`Assembly Accession`)),]
dataset <- dataset[, colSums(is.na(dataset)) <= 400]

colnames(dataset)[7] <- "BioSampleAccession"
colnames(metadata)[1] <- "BioSampleAccession"

metaphylogeny <- inner_join(dataset, metadata, by = "BioSampleAccession")

for (i in 1:length(metaphylogeny$`Assembly Accession`)) {
  metaphylogeny$`Assembly Accession`[i] <- unlist(strsplit(metaphylogeny$`Assembly Accession`[i], "[.]"))[1]
}

# number of BGCs per genome----
# unique genome IDs for matching
for (i in 1:length(rerooted_outgroup_droptip[["tip.label"]])) {
  rerooted_outgroup_droptip[["tip.label"]][i] <- unlist(strsplit(rerooted_outgroup_droptip[["tip.label"]][i], "[.]"))[1]
}

# BGCs (GECCO, AS) in each genome
htgcf <- read_delim("~/Uni/11.masters_thesis/Clustering/htgcf_clusters.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
htgcf <- select(htgcf, - c("filename")) #remove filename column

# select GECCO and AS BGCs:
clusters <- htgcf[which(grepl(paste(c("region", "cluster"), collapse = "|"), htgcf$cluster_id)),]
clusters$genome <- clusters$cluster_id
str1 <- sapply(strsplit(clusters$genome, "_"), `[`, 1)
str2 <- sapply(strsplit(clusters$genome, "_"), `[`, 2)
clusters$genome <- paste(str1, str2, sep = "_")

BGCS_per_genome <- data.frame(matrix(ncol = 3, nrow = 458))

colnames(BGCS_per_genome) <- c("genome", "GECCO", "antiSMASH")
BGCS_per_genome$genome <- rerooted_outgroup_droptip[["tip.label"]]

for (i in 1:length(BGCS_per_genome$genome)) {
  BGCS_per_genome$antiSMASH[i] <- sum(grepl("antiSMASH", unlist(as.vector(clusters[which(grepl(BGCS_per_genome$genome[i], clusters$cluster_id)),][,1]))))
  BGCS_per_genome$GECCO[i] <- sum(grepl("cluster", unlist(as.vector(clusters[which(grepl(BGCS_per_genome$genome[i], clusters$cluster_id)),][,1]))))
}

for (i in 1:length(BGCS_per_genome$genome)) {
  BGCS_per_genome$totalBGCs[i] <- BGCS_per_genome$GECCO[i] + BGCS_per_genome$antiSMASH[i]
}


ggtree(rerooted_outgroup_droptip) + 
  geom_treescale(x = 0.005, y = 0, fontsize = 3) +
  geom_facet(panel = "GECCO", data = BGCS_per_genome, 
             geom = geom_col, aes(x = GECCO, fill = "GECCO"), 
             color = NA, size = 0, width = 1, orientation = 'y') +
  geom_facet(panel = "antiSMASH", data = BGCS_per_genome, 
             geom = geom_col, aes(x = antiSMASH, fill = "antiSMASH"), 
             color = NA, size = 0, width = 1, orientation = 'y') +
  geom_facet(panel = "total BGCs", data = BGCS_per_genome, 
             geom = geom_col, aes(x = totalBGCs, fill = "totalBGCs"), 
             color = NA, size = 0, width = 1, orientation = 'y') +
  scale_fill_manual(values = c("GECCO" = "#5A7E33", "antiSMASH" = "#AA0000", "totalBGCs" = "black")) +
    theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        strip.text.x = element_text(size = 20)) +
  ggtitle("Number of BGCs detected per genome by each method")
  


# adding metadata----
  # empty tree
ggtree(rerooted_outgroup_droptip) + 
  geom_treescale(x = 0.005, y = 0, fontsize = 3)
  #tree + BGCs per genome
ggtree(rerooted_outgroup_droptip) + 
  geom_treescale(x = 0.005, y = 0, fontsize = 3) + 
  geom_fruit(data = BGCS_per_genome, geom = geom_col, 
             mapping = aes(x = GECCO, y = genome), fill = "#73A790", 
             color = "#73A790") +
  geom_fruit(data = BGCS_per_genome, geom = geom_col, 
             mapping = aes(x = antiSMASH, y = genome), fill = "#EABAB9", 
             color = "#EABAB9")

  #BGCs + continent + isolation source
metaphylogeny <- as.data.frame(metaphylogeny)
rownames(metaphylogeny) <- metaphylogeny[,1]
colnames(metaphylogeny)[1] <- "genome"

ggtree(rerooted_outgroup_droptip) + 
  geom_treescale(x = 0.005, y = 0, fontsize = 3) + 
  geom_fruit(data = BGCS_per_genome, geom = geom_col, 
             mapping = aes(x = GECCO, y = genome), fill = "#73A790", 
             color = "#73A790") +
  geom_fruit(data = BGCS_per_genome, geom = geom_col, 
             mapping = aes(x = antiSMASH, y = genome), fill = "#EABAB9", 
             color = "#EABAB9") +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = 1, y = genome, fill = continent), width = 0.0007,
             offset = 0, pwidth = 0.05) + 
  scale_fill_viridis_d(option = "D", name="Continent", na.value = "gray") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = 1, y = genome, fill = `isolation source`), 
             width = 0.0007) +
  scale_fill_viridis_d(option = "A", name="Isolation source", na.value = "gray") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = 1, y = genome, fill = `host disease`), 
             width = 0.0007) +
  scale_fill_viridis_d(option = "H", name="Isolation source", na.value = "gray")

# phylogeny + MIBiG mixed GCFs----
#7 GCFs have MIBiG and antiSMASH/GECCO BGCs -> absence/presence heatmap per genome
library(ggsci)
mixed_mibig_clusters <- c("GCF0000070", "GCF0000073", "GCF0000075", "GCF0000082",
                          "GCF0000092", "GCF0000100", "GCF0000923")
for(i in mixed_mibig_clusters) {
  for (x in 1:length(clusters$gcf_id)) {
    if(i == clusters$gcf_id[x]) {
      ind <- which(metaphylogeny$genome == clusters$genome[x])
      metaphylogeny[ind,i] <- i
    }
  }
}

ggtree(rerooted_outgroup_droptip) + 
  geom_treescale(x = 0.005, y = 0, fontsize = 3) + 
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000070), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "darkblue", na.value = "white") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000073), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "darkred", na.value = "white") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000075), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "darkgreen", na.value = "white") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000082), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "darkorange", na.value = "white") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000092), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "darkorchid", na.value = "white") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000100), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "gold3", na.value = "white") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000923), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "darkturquoise", na.value = "white") +
  theme(legend.position = "right")


for(i in 1:length(clusters$gcf_id)){
  if(clusters$gcf_id[i] == "GCF0000084"){
    ind <- which(metaphylogeny$genome == clusters$genome[i])
    metaphylogeny$GCF_phylogeny[ind] <- "GCF0000084"
  }
}

for(i in 1:length(clusters$gcf_id)){
  if(clusters$gcf_id[i] == "GCF0000084"){
    ind <- which(metaphylogeny$genome == clusters$genome[i])
    metaphylogeny$GCF_phylogeny[ind] <- "GCF0000084"
  }
}

ggtree(rerooted_outgroup_droptip) +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(y = genome, fill = GCF_phylogeny),
             width = 0.002, pwidth = 0.05)

# novel GCFs----
# GCFs that do not contain a MIBiG BGC
novel_GCFs <- c("GCF0000099", "GCF0000091", "GCF0000103", "GCF0000107", "GCF0000088")
for (i in novel_GCFs) {
  for (x in 1:length(clusters$gcf_id)) {
    if(i == clusters$gcf_id[x]) {
      ind <- which(metaphylogeny$genome == clusters$genome[x])
      metaphylogeny[ind,i] <- i
    }
  }
}

ggtree(rerooted_outgroup_droptip) + 
  geom_treescale(x = 0.005, y = 0, fontsize = 3) + 
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000099), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "#E2725B", na.value = "gray") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000091), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "#808000", na.value = "gray") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000103), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "#A0522D", na.value = "gray") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000107), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "#FFDB58", na.value = "gray") + new_scale_fill() +
  geom_fruit(data = metaphylogeny, geom = geom_tile,
             mapping = aes(x = , y = genome, fill = GCF0000088), width = 0.0002, pwidth = 0.05) +
  scale_fill_manual(values = "#708090", na.value = "gray") +
  theme(legend.position = "none")

gcf0000084 <- subset(clusters, gcf_id == "GCF0000084")

ggtree(rerooted_outgroup_droptip) +
  geom_treescale(x = 0.005, y = 0, fontsize = 3) +
  geom_fruit(data = gcf0000084, geom = geom_tile,
             mapping = aes(x = , y = genome), width = 0.0002, pwidth = 0.05,
             fill = "#A13330") +
  labs(title = "Presence of GCF0000084") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18))

# Legends----
library(cowplot)
library(gridExtra)
library(ggplotify)
  # decoy BGC legend
BGC_legend <- get_legend(ggplot(data = data.frame(c(NA, NA), c(NA,NA))) + 
  geom_tile(aes(x = NA, y = NA, fill = c("Number of antiSMASH BGCs", "Number of GECCO BGCs"))) +
  scale_fill_manual(values = c("#EABAB9","#73A790"),
                     guide = guide_legend(title = "BGC Detection Method")))
  # metadata legend
continent_legend <- get_legend(ggtree(rerooted_outgroup_droptip) + 
                            geom_treescale(x = 0.005, y = 0, fontsize = 3) + 
                            geom_fruit(data = metaphylogeny, geom = geom_tile,
                                       mapping = aes(x = 1, y = genome, fill = continent), width = 0.0007,
                                       offset = 0, pwidth = 0.05) + 
                            scale_fill_viridis_d(option = "A", name="Continent", na.value = "gray"))

source_legend <- get_legend(ggtree(rerooted_outgroup_droptip) + 
                              geom_treescale(x = 0.005, y = 0, fontsize = 3) +
                              geom_fruit(data = metaphylogeny, geom = geom_tile,
                                         mapping = aes(x = 1, y = genome, fill = `isolation source`), 
                                         width = 0.0007) +
                              scale_fill_viridis_d(option = "D", name="Isolation source", na.value = "gray"))

disease_legend <- get_legend(ggtree(rerooted_outgroup_droptip) + 
                              geom_treescale(x = 0.005, y = 0, fontsize = 3) +
                              geom_fruit(data = metaphylogeny, geom = geom_tile,
                                         mapping = aes(x = 1, y = genome, fill = `host disease`), 
                                         width = 0.0007) +
                              scale_fill_viridis_d(option = "C", name="Host disease", na.value = "gray"))

  # decoy GCF subset legends
mixed_MIBiG_legend <- get_legend(ggplot(data = data.frame(c(NA, NA, NA, NA, NA, NA, NA), c(NA, NA, NA, NA, NA, NA, NA))) + 
                                  geom_tile(aes(x = NA, y = NA, fill = c("GCF0000070", "GCF0000073", "GCF0000075", "GCF0000082",
                                                                         "GCF0000092", "GCF0000100", "GCF0000923"))) +
                                  scale_fill_manual(values = c("darkblue","darkred","darkgreen","darkorange","darkorchid","gold3","darkturquoise"),
                                                    guide = guide_legend(title = "MIBiG + antiSMASH/GECCO GCFs")))

  # novel GCFs (do not contain a MIBiG BGC)
novel_GCFs_legend <- get_legend(ggplot(data = data.frame(c(NA, NA, NA, NA, NA), c(NA, NA, NA, NA, NA))) + 
                                  geom_tile(aes(x = NA, y = NA, fill = c("GCF0000099", "GCF0000091", "GCF0000103", "GCF0000107", "GCF0000088"))) +
                                  scale_fill_manual(values = c("#E2725B", "#808000", "#A0522D", "#FFDB58", "#708090"),
                                                    guide = guide_legend(title = "novel GCFs (do not contain a MIBiG BGC)")))

# gheatmaps----

pcontinent <- gheatmap(ggtree(rerooted_outgroup_droptip) +
                         geom_treescale(x = 0.005, y = 0, fontsize = 3) +
                         theme(legend.position = "none",
                               plot.title = element_text(size = 20, hjust = 0.5)) +
                         labs(title = expression(paste(italic("Streptococcus mutans"), " phylogeny"))), subset(phymeta, select = continent), 
                   colnames = FALSE, width = 0.15, color = NULL) +
            scale_fill_viridis_d(option = "A", name="Continent") + new_scale_fill()


gheatmap((gheatmap(pcontinent, subset(phymeta, select = `host disease`), 
                   offset = 0.0015, colnames = FALSE, width = 0.15, color = NULL) +
            scale_fill_viridis_d(option = "C", name="Host Disease") + new_scale_fill()), 
         subset(phymeta, select = `isolation source`), colnames = FALSE, width = 0.15, color = NULL, offset=0.003) +
  scale_fill_viridis_d(option = "D", name="Isolation Source") + NoLegend()

----

# Panstripe without outgroup----
fit_nooutgroup <- panstripe(gene_pa_nooutgroup, midpointroot_nooutgroup)
fit_nooutgroup$summary

plot_pangenome_params(fit_nooutgroup)
plot_pangenome_cumulative(fit_nooutgroup)
plot_pangenome_branches(fit_nooutgroup)
plot_pangenome_curve(fit_nooutgroup)
plot_gain_loss(fit_nooutgroup)


