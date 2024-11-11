# make counts matrix----------
library(reticulate)

use_python("~/.virtualenvs/r-reticulate/Scripts/python.exe")
scipy_sparse <- import("scipy.sparse")

counts = scipy_sparse$load_npz("C:/Users/marle/Documents/Uni/11.masters_thesis/umap/compositions.npz")
dim(counts)
colnames(counts) <- read.table("~/Uni/11.masters_thesis/umap/domains.tsv")[,1]
rownames(counts) <- read.table("~/Uni/11.masters_thesis/umap/labels.tsv")[,1]
saveRDS(counts, "~/Uni/11.masters_thesis/umap/counts.RDS")


# making the UMAP-----
counts <- readRDS("~/Uni/11.masters_thesis/umap/counts.RDS")

library(Seurat)
  #QC, analysis, and exploration of single-cell RNA-seq data
  #needed for CreateChromatinAssay()
library(Signac)
  #analysis of single-cell chromatin data
  #needed for CreateSeuratObject() and subsequent functions

colnames(counts) <- paste0(colnames(counts),":1-10")

counts <- as.matrix(counts)

chrom_assay <- CreateChromatinAssay(
  counts = t(counts),
  sep = c(":", "-")
  #min.cells = 10,
  #min.features = 200
)

smutans <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

smutans <- RunTFIDF(smutans) #Term Frequency-Inverse Document Frequency
  #normalizes feature counts by frequency -> features are weighted based on importance
smutans <- FindTopFeatures(smutans, min.cutoff = 'q0') #selects top features that are most informative and variable
  #reduces noise by filtering out less informative features
  #mincutoff q0: features with non-zero expression will be retained
  #output: reduced set of top features
smutans <- RunSVD(smutans) #Singular Value Decomposition: dimensionality reduction
  #decompose matrix into orthogonal components
  #reduces dimensions while preserving most important patterns

DepthCor(smutans) #We should include dim=1 in this case!
  #normalization to adjust for sequencing depth

smutans <- RunUMAP(object = smutans, reduction = 'lsi', dims = 1:30, 
                   n.components = 2L, n.neighbors = 30L, min.dist = 0.3)
  #UMAP will be applied to the data after dimensionality reduction using latent semantic indexing (LSI) obtained from SVD
  #uses first 30 dimensions from SVD
    #playing around with dimensions (usually taken from PCA dimensions)
    #we're not doing that -> making sure that the clusters are stable (around 20 dims)
smutans_tree <- FindNeighbors(object = smutans, reduction = 'lsi', dims = 1:30)
  #finds nearest neighbors for each BGC based on low-dimensional representations obtained from LSI

library(clustree)
res = seq.int(0.1, 0.9, 0.05) 
for (i in res){
  smutans_tree = FindClusters(object = smutans_tree, resolution = i, algorithm = 3)
}
clustree(smutans_tree@meta.data, prefix = "peaks_snn_res.", highlight_core = TRUE)

smutans <- FindClusters(object = smutans_tree, verbose = FALSE, algorithm = 3, resolution = 0.5)
  #assigns BGCs into clusters based on similarity in low-dimensional space
  #algorithm used: SLM (smart local moving) algorithm: maximizes a so-called modularity function
  #why this one? -> says that Leiden (4, requires something) is better

DimPlot(object = smutans, label = FALSE, pt.size = 1) +
  labs(title = "Seurat Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 20),
        axis.title=element_text(size=18))

# Adding Metadata to smutans object----
library(readr)
gecco_bgcs <- read_delim("~/Uni/11.masters_thesis/Clustering/gecco_clusters_merged.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
# remove unwanted rows
gecco_bgcs <- subset(gecco_bgcs, sequence_id != "sequence_id")

addmeta <- htgcf[,c(1,2,3,4)]
colnames(addmeta) <- c("BGC", "cluster_length", "GCF", "GCF_rep")
addmeta[,c("class", "method", "genome", "GCF_method")] <- NA

for (i in 1:length(addmeta$BGC)) {
  if (grepl("antiSMASH", addmeta$BGC[i])){
    addmeta$method[i] <- "antiSMASH"
  } else if (grepl("cluster", addmeta$BGC[i])) {
    addmeta$method[i] <- "gecco"
  } else {
    addmeta$method[i] <- "MIBiG"
  }
}

#adding gecco data
for (i in 1:length(gecco_bgcs$cluster_id)) {
  id <- which(grepl(gecco_bgcs$cluster_id[i], addmeta$BGC))
  addmeta$class[id] <- gecco_bgcs$type[i]
}

#adding antiSMASH data
antismash_info <- read_csv("~/Uni/11.masters_thesis/umap/antismash_info.csv")
str_1 <- sapply(strsplit(antismash_info$`File Name`, "[.]"), `[`, 1)
str_2 <- sapply(strsplit(antismash_info$`File Name`, "[.]"), `[`, 2)
antismash_info$`File Name` <- paste(str_1, str_2, sep = ".")

antismash_info$`Protocluster Classes` <- gsub("'|\\[|\\]","", antismash_info$`Protocluster Classes`)
antismash_info$`Protocluster Classes` <- strsplit(antismash_info$`Protocluster Classes`, ", ")
for (i in 1:length(antismash_info$`Protocluster Classes`)) {
  if(antismash_info$`Protocluster Count`[i]>1){
    if(length(unique(unlist(antismash_info$`Protocluster Classes`[i])))>1){
      antismash_info$`Protocluster Classes`[i] <- "mixed"
    } else if(length(unique(unlist(antismash_info$`Protocluster Classes`[i])))==1){
      antismash_info$`Protocluster Classes`[i] <- unlist(antismash_info$`Protocluster Classes`[i])[1]
    }
  }
}

for (i in 1:length(antismash_info$`File Name`)) {
  id <- which(grepl(antismash_info$`File Name`[i], addmeta$BGC))
  addmeta$class[id] <- antismash_info$`Protocluster Classes`[i]
}

#adding MIBiG data
mibig_info <- read_csv("~/Uni/11.masters_thesis/umap/mibig_info.csv")
mibig_info$`File Name` <- sapply(strsplit(mibig_info$`File Name`, "[.]"), `[`, 1)
mibig_info$`Protocluster Classes` <- gsub("'|\\[|\\]","", mibig_info$`Protocluster Classes`)

mibig_info$`Protocluster Classes`[mibig_info$`Protocluster Classes`==""] <- "unknown"
mibig_info$`Protocluster Classes` <- strsplit(mibig_info$`Protocluster Classes`, ", ")
for (i in 1:length(mibig_info$`Protocluster Classes`)) {
  if(mibig_info$`Protocluster Count`[i]>1){
    if(length(unique(unlist(mibig_info$`Protocluster Classes`[i])))>1){
      mibig_info$`Protocluster Classes`[i] <- "mixed"
    } else if(length(unique(unlist(mibig_info$`Protocluster Classes`[i])))==1){
      mibig_info$`Protocluster Classes`[i] <- unlist(mibig_info$`Protocluster Classes`[i])[1]
    }
  }
}

for (i in 1:length(mibig_info$`File Name`)) {
  id <- which(grepl(mibig_info$`File Name`[i], addmeta$BGC))
  addmeta$class[id] <- mibig_info$`Protocluster Classes`[i]
}

# tweaking data frame
for(i in unique(addmeta$GCF_rep)){
  indexes <- which(addmeta$GCF_rep==i)
  methods <- unlist(addmeta$method[indexes])
  if(sum(grepl("MIBiG", methods))>0){
    addmeta$GCF_method[indexes] <- "MIBiG"
  } else if ("gecco" %in% methods & "antiSMASH" %in% methods){
    addmeta$GCF_method[indexes] <- "mixed"
  } else if ("gecco" %in% methods){
    addmeta$GCF_method[indexes] <- "Gecco"
  } else if ("antiSMASH" %in% methods){
    addmeta$GCF_method[indexes] <- "antiSMASH"
  } else {
    addmeta$GCF_method[indexes] <- "wrong"
  }
}

addmeta <- as.data.frame(addmeta)
rownames(addmeta) <- addmeta$BGC
addmeta[addmeta$class=="Unknown","class"] <- "unknown"
addmeta[addmeta$class=="NRP;Polyketide", "class"] <- "mixed"
addmeta[addmeta$class=="Polyketide","class"] <- "PKS"
addmeta[addmeta$class=="other","class"] <- "unknown"
addmeta[addmeta$class=="NRP","class"] <- "NRPS"

# adding GCF count information
GCF_info <- as.data.frame(matrix(ncol=3, nrow=1640))
colnames(GCF_info) <- c("GCF", "#ofBGCs", "genomes with BGC (%)")
rownames(GCF_info) <- rownames(smutans@meta.data)
GCF_info$GCF <- smutans$GCF

for (i in 1:length(GCF_info$GCF)) {
  GCF_info$`#ofBGCs`[i] <- sum(addmeta$GCF == GCF_info$GCF[i])
}

for (i in 1:length(GCF_info$GCF)) {
  GCF_info$`genomes with BGC (%)`[i] <- length(unique(subset(
    clusters, gcf_id == GCF_info$GCF[i])$genome))/458 * 100
}

# AddMetaData
smutans <- AddMetaData(smutans, addmeta)
smutans <- AddMetaData(smutans, GCF_info)

# Colour by data----
DimPlot(object = smutans, label = FALSE)
DimPlot(object = smutans, label = FALSE, group.by = "class")
DimPlot(object = smutans, label = FALSE, group.by = "GCF_method")

DimPlot(object = smutans, label = FALSE, group.by = "cluster_length") + NoLegend()

table(smutans$class, smutans$seurat_clusters)

DimPlot(object = smutans, label = FALSE, group.by = "class", pt.size = 2.5) + 
  theme(plot.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_discrete(name = "Biosynthetic Class", labels = c("mixed", "NRPs", "PKS", 
                                                              "RiPPs", "Saccharides", "Terpenes",
                                                              "other")) +
  labs(title = "UMAP clustering of GCFs") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 28),
        legend.text = element_text(size = 20),
        legend.position.inside = c(0.8,0.9))

DimPlot(object = smutans, label = FALSE, group.by = "class", pt.size = 2.5) + 
  scale_color_discrete(name = "Biosynthetic Class", 
                       labels = c("mixed", "NRPs", "PKS", "RiPPs", "Saccharides", "Terpenes", "other")) +
  labs(title = "UMAP clustering of GCFs") + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 40),
    legend.text = element_text(size = 32),
    legend.position = "right",
    legend.title = element_text(size = 36),
    plot.margin = unit(c(1, 5, 1, 1), "cm") # Add space on the right
  )

# GCFs that contain both MIBiG and GECOO/AS BGCs----
mixed_GCFs <- data.frame(matrix(ncol = 2, nrow = 7))
colnames(mixed_GCFs) <- c("GCF", "Number of BGCs")

for (i in unique(addmeta$GCF)) { #adding GCF and number of BGCs per GCF
  if(sum(addmeta[addmeta$GCF==i,]$method=="MIBiG")>0){
    if(sum(addmeta[addmeta$GCF==i,]$method=="gecco")>0 | sum(addmeta[addmeta$GCF==i,]$method=="antiSMASH")>0){
      mixed_GCFs$GCF[x] <- i
      mixed_GCFs$`Number of BGCs`[x] <- sum(addmeta$GCF==i)
      x <- x + 1 
    }
  }
}

GCF0000073 <- addmeta$BGC[which(addmeta$GCF=="GCF0000073")]

#which MIBiG BGCs are in each GCF
for (i in 1:length(mixed_GCFs$GCF)) {
  mixed_GCFs$MIBiG_clusters[i] <- list(subset(addmeta, GCF == mixed_GCFs$GCF[i] & method == "MIBiG")$BGC)
}

#number of S. mutans genomes that contain at least one of the BGCs
for (i in 1:length(mixed_GCFs$GCF)) {
  mixed_GCFs$`genomes that contain one BGC (%)`[i] <- 
    length(unique(subset(clusters, cluster_id %in% subset(addmeta, GCF==mixed_GCFs$GCF[i])$BGC)$genome))
}

for (i in 1:length(mixed_GCFs$`genomes that contain one BGC (%)`)) {
  mixed_GCFs$`genomes that contain one BGC (%)`[i] <- mixed_GCFs$`genomes that contain one BGC (%)`[i]/458 * 100
}

# novel GCFs----

novel_GCFs <- addmeta[which(addmeta$GCF_method != "MIBiG"),]
colnames(novel_GCFs) <- c("cluster_id", "cluster_length", 
                          "gcf_id", "gcf_representative", "class", "method", "GCF_method")
novel_GCFs <- merge(novel_GCFs, clusters, all.x = TRUE, 
                    by.x = c("cluster_id", "cluster_length", "gcf_id", "gcf_representative"),
                    by.y = c("cluster_id", "cluster_length", "gcf_id", "gcf_representative"))


GCA_019048645_BGCs_novel <- novel_GCFs[which(novel_GCFs$genome == "GCA_019048645"),]$gcf_representative
GCA_019048645_GCFs <- unique(subset(novel_GCFs, genome == "GCA_019048645")$gcf_id)
  
for(i in GCA_019048645_GCFs){
  print(paste(i,": ", sum(clusters$gcf_id == i)))
}


ggplot(subset(smutans@meta.data, GCF %in% novelGCFs), 
       aes(x = factor(class, level = c("NRPS", "PKS", "RiPP", "mixed", "unknown")),
           fill = GCF_method)) +
  geom_bar() +
  scale_fill_manual(values = c("#AA0000", "#5A7E33", "#2A4765")) +
  labs(title = "Novel GCF predicted biosynthetic classes", y = "Number of GCFs", x = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0, size = 20)) +
  scale_x_discrete(label = c("NRP", "PKS", "RiPP", "mixed", "other"))


#clinker samples


#3D UMAP----
TS_plot_3D <- function (seurat_object, dims = 1:30, dot_size) {
  require(plotly)
  require(scales)
  temp1 <- RunUMAP(object = seurat_object, dims = dims, n.components = 3,
                   verbose = F, reduction='lsi')
  df <- data.frame(umap1 = temp1@reductions$umap@cell.embeddings[,
                                                                 1], umap2 = temp1@reductions$umap@cell.embeddings[, 2],
                   umap3 = temp1@reductions$umap@cell.embeddings[, 3], cell = Idents(temp1))
  plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~cell,
          colors = hue_pal()(length(levels(df$cell))), marker = list(size = dot_size)) %>% 
    add_markers() %>%
    layout(scene = list(xaxis = list(title = "umap 1"), yaxis = list(title = "umap 2"),
                        zaxis = list(title = "umap 3")))
}

TS_plot_3D(smutans, dims = 1:30, 5)


