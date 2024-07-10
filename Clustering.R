library(readr)
library(ggplot2)
library(stringr)
library(dplyr)
library(pheatmap)

# BGCs and classes per method----
 # antiSMASH----
antismash_info <- read_csv("~/Uni/11.masters_thesis/umap/antismash_info.csv")
#clean up dataframe
str_1 <- sapply(strsplit(antismash_info$`File Name`, "[.]"), `[`, 1)
str_2 <- sapply(strsplit(antismash_info$`File Name`, "[.]"), `[`, 2)
antismash_info$`File Name` <- paste(str_1, str_2, sep = ".") #remove .gbk

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

antismash_info$`Protocluster Classes` <- unlist(antismash_info$`Protocluster Classes`)
antismash_info[antismash_info=="NRPS"] <- "NRP"

# plot: BGC classes predicted by antiSMASH
ggplot(antismash_info, aes(x = factor(`Protocluster Classes`, 
                                      level = c("NRP", "PKS", "RiPP", "mixed", "other")))) + 
  geom_bar(fill = "#AA0000") +
  labs(title = "BGC classes predicted by antiSMASH",
       y = "Number of BGCs",
       x = "") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2, size = 6) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 2600))

 # GECCO----
gecco_bgcs <- read_delim("~/Uni/11.masters_thesis/Clustering/gecco_clusters_merged.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
# remove unwanted rows/clean up dataframe
gecco_bgcs <- subset(gecco_bgcs, sequence_id != "sequence_id")
gecco_bgcs[gecco_bgcs=="Unknown"] <- "other"
gecco_bgcs[gecco_bgcs=="NRP;Polyketide"] <- "mixed"
gecco_bgcs[gecco_bgcs=="Polyketide"] <- "PKS"

ggplot(gecco_bgcs, aes(x = factor(`Protocluster Classes`, 
                                  level = c("NRP", "PKS", "RiPP", "mixed", "other")))) +
  geom_bar(fill = "#5A7E33") +
  labs(title = "BGC classes predicted by GECCO",
       y = "Number of BGCs",
       x = "") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2, size = 6) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 2600))

 # MIBiG----
mibig_info <- read_csv("~/Uni/11.masters_thesis/umap/mibig_info.csv")
#clean up dataframe
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

mibig_info$`Protocluster Classes` <- unlist(mibig_info$`Protocluster Classes`)
mibig_info[mibig_info=="unknown"] <- "other"
mibig_info[mibig_info=="NRPS"] <- "NRP"

#plot: BGC classes in MIBiG
ggplot(mibig_info, aes(x = `Protocluster Classes`)) + 
  geom_bar(fill = "#2A4765") +
  labs(title = "BGC classes of validated MIBiG BGCs",
       y = "Number of BGCs",
       x = "") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2) +
  theme(plot.title = element_text(hjust = 0.5))

 # Combined plots----
antismash_info$method <- "antiSMASH"
mibig_info$method <- "MIBiG"
gecco_bgcs$method <- "GECCO"
colnames(gecco_bgcs)[1] <- "File Name"
colnames(gecco_bgcs)[7] <- "Protocluster Classes"

BGC_classes <- rbind(antismash_info, mibig_info)
BGC_classes$`Protocluster Count` <- NULL
gecco_info <- gecco_bgcs[,c(1,7,16)]
BGC_classes <- rbind(BGC_classes, gecco_info)

# combined graph
ggplot(BGC_classes, aes(x = `Protocluster Classes`)) + 
  geom_bar(fill = "#000000") +
  labs(title = "BGC classes",
       y = "Number of BGCs", x = "") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2) +
  theme(plot.title = element_text(hjust = 0.5))

# colored by method
ggplot(BGC_classes, aes(x = `Protocluster Classes`, fill = method)) + 
  geom_bar() +
  scale_fill_manual(values=c("#EABAB9", "#73A790", "#2A4765")) +
  labs(title = "BGC classes",
       y = "Number of BGCs", x = "") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5))

# GCF clustering with HTGCF----
 # making the data frame
htgcf <- read_delim("~/Uni/11.masters_thesis/Clustering/htgcf_clusters.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
htgcf <- select(htgcf, - c("filename")) #remove filename column

MIBiG <- htgcf[which(grepl("BGC000", htgcf$cluster_id)),] #2515 MIBiG BGCs

sum(grepl("BGC000", htgcf$cluster_id)) #2502 MIBiG BGCs
sum(grepl("region", htgcf$cluster_id)) #3708 antiSMASH BGCs
sum(grepl("cluster", htgcf$cluster_id)) #2849 GECCO BGCs
length(unique(htgcf$gcf_representative)) #1640 GCFs

 # extract genome IDs:
for (i in 1:length(htgcf$cluster_id)) {
  if (grepl("GCA", htgcf$cluster_id[i])){
    htgcf$genome[i] <- paste(unlist(strsplit(htgcf$cluster_id[i], "_"))[1],
                             unlist(strsplit(htgcf$cluster_id[i], "_"))[2], sep = "_")
  }
  else if(grepl("BGC000", htgcf$cluster_id[i])){
    htgcf$genome[i] <- NA
  }
  else {
    htgcf$genome[i] <- htgcf$cluster_id[i]
  }
}

 # Plot: cluster length
ggplot(htgcf, aes(x = `cluster_length`)) + 
  geom_histogram(fill = "#2A4765", binwidth = 10000) +
  labs(title = "Cluster length",
       y = "Number of BGCs", x = "Cluster length") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5))

 # Plot: number of BGCs per GCF
count_gcf <- htgcf %>%
  group_by(gcf_id) %>%
  summarise(cluster_count = n())

ggplot(count_gcf, aes(x = gcf_id, y = cluster_count)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of BGCs per GCF",
       x = "",
       y = "Number of BGCss") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank())

#HTGCF GCFs that contain one or more MIBiG GCFs + our BGCs----------------------
library(gridExtra)
library(grid)

length(unique(htgcf$gcf_id)) #1640 GCFs

htgcf_gcf_grouped <- htgcf %>%
  count(gcf_id)

for (i in 1:length(htgcf_gcf_grouped$gcf_id)) {
  rows <- as.numeric(which(htgcf$gcf_id == htgcf_gcf_grouped$gcf_id[i]))
  htgcf_gcf_grouped$bgcs[i] <- list(htgcf$cluster_id[rows])
}

#clusters that contain both a MiBIG and a GECCO or antiSMASH BGC
mixed_clusters <- htgcf_gcf_grouped$gcf_id[which(grepl("BGC00", htgcf_gcf_grouped$bgcs) & 
                                 grepl("region|cluster", htgcf_gcf_grouped$bgcs))]

mixed_clusters_df <- subset(htgcf_gcf_grouped, grepl(paste(mixed_clusters, collapse = "|"), gcf_id))

grid.table(mixed_clusters_df[,1:2])




find_positions <- function(string, df_column) {
  positions <- which(df_column == string)
  if (length(positions) == 0) {
    return(NA)  # Return NA if the string is not found
  } else {
    return(positions)
  }
}

# Apply the function to each string in the vector
df <- sapply(smutans_GCFs, function(x) find_positions(x, count_gcf[order(-count_gcf$cluster_count), ]$gcf_id))




  