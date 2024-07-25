library(readr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(countrycode)

metadata <- read_csv("Uni/11.masters_thesis/metadata_phylogeny/biosample_table.csv")

# remove empty rows:
metadata <- metadata[rowSums(is.na(metadata)) != ncol(metadata), ]

# how many NAs in each column
sapply(metadata, function(x) sum(is.na(x)))

# removing all columns that have more than 300 NAs
metadata <- metadata[, colSums(is.na(metadata)) <= 300]

# geographic location----
# split location -> remove cities/states -> only country
metadata$`geographic location` <- sapply(strsplit(metadata$`geographic location`,":"), `[`, 1)
# change other missing locations to NA
metadata$`geographic location`[grep("not provided", metadata$`geographic location`)] <- NA
metadata$`geographic location`[grep("Not applicable", metadata$`geographic location`)] <- NA
metadata$`geographic location`[grep("Missing", metadata$`geographic location`)] <- NA
# group by continent
metadata$continent <- countrycode(sourcevar = metadata$`geographic location`,
                            origin = "country.name",
                            destination = "continent")

ggplot(metadata, aes(x = `geographic location`, fill = continent)) + 
  geom_bar(stat = "count") +
  theme_light() +
  labs(title = "Geographic distribution",
       x = "Location",
       y = "Number of genomes",
       fill = "Continent") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18))

# host disease-------------------------------------------------------------------
# split at ; to remove unnecessary information
metadata$`host disease` <- sapply(strsplit(metadata$`host disease`,";"), `[`, 1)
# group different spellings
metadata$`host disease` <- metadata$`host disease` %>%
  str_to_lower()
metadata$`host disease`[grep(paste(c("none", "not applicable", "unknown"), collapse="|"), 
                                      metadata$`host disease`)] <- NA
metadata$`host disease`[grep("caries", metadata$`host disease`, ignore.case = TRUE)] <- "dental caries"

ggplot(metadata, aes(x = `host disease`, fill = `host disease`)) + 
  geom_bar(show.legend = FALSE) +
  theme_light() +
  labs(title = "Host disease",
       y = "Number of genomes",
       x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18))

# isolation source---------------------------------------------------------------
metadata$`isolation source` <- metadata$`isolation source` %>%
  str_to_lower()
metadata$`isolation source`[grep(paste(c("not available", "not provided"), collapse="|"), 
                              metadata$`isolation source`)] <- NA

# group sources together
metadata$`isolation source`[grep(paste(c("caries", "carious"), collapse = "|")
                                  , metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- "caries"
metadata$`isolation source`[grep(paste(c("plaque", "teeth", "dental", "tooth", "guinea"), collapse = "|")
                                  , metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- "plaque"
metadata$`isolation source`[grep("gut", metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- "gut microbiome"
metadata$`isolation source`[grep("endocarditis", metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- "endocarditis"
metadata$`isolation source`[grep(paste(c("saliva", "tongue", "oral", "root"), collapse = "|")
                                  , metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- "oral cavity"
metadata$`isolation source`[grep(paste(c("feces", "stool"), collapse = "|")
                                  , metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- "feces"
metadata$`isolation source`[grep("atcc", metadata$`isolation source`, 
                                  ignore.case = TRUE)] <- NA

ggplot(metadata, aes(x = `isolation source`, fill = `isolation source`)) +
  geom_bar() +
  theme_light() +
  labs(title = "Isolation source", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18))

# collection date----------------------------------------------------------------
metadata$`collection date`[grep(paste(c("Unknown", "Missing", "Not collected", "Not applicable", "not provided"), collapse="|"), 
                                  metadata$`collection date`)] <- NA

# clean up data: two dates == NA
for (i in 1:length(metadata$`collection date`)) {
  if (grepl("/", metadata$`collection date`[i])){
    metadata$`collection date`[i] <- NA
    } else if (grepl("-", metadata$`collection date`[i])){
      metadata$`collection date`[i] <- unlist(strsplit(metadata$`collection date`[i],"-"))[1]
      } else {
        metadata$`collection date`[i] <- metadata$`collection date`[i]
      }
  }

metadata$`collection date`[which(metadata$`collection date` > 2023)] <- NA
colfunc <- colorRampPalette(c("#e8ffff", "#000137"))

ggplot(metadata, aes(x = `collection date`, fill = `collection date`)) +
  geom_bar(show.legend = FALSE) +
  theme_light() +
  labs(title = "Collection date",
       y = "Count",
       x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = colfunc(36)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18))

save(metadata, file = "Uni/11.masters_thesis/metadata/metadata.Rds")

