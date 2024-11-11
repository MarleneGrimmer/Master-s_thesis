library(readr)
library(ggplot2)
library(dplyr)

#CheckM--------------------------------------------------------------------------
checkm_out <- read_delim("~/Uni/11.masters_thesis/quality_assessment/checkm_out.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

#which marker lineages were used
unique(checkm_out$`Marker lineage`)
sum(checkm_out$`Marker lineage`=="o__Lactobacillales (UID544)") # -> 475
sum(checkm_out$`Marker lineage`=="k__Bacteria (UID203)") # -> 1
which(checkm_out$`Marker lineage`=="k__Bacteria (UID203)") # one genome couldn't be classified in more detail

ggplot(checkm_out, aes(x = Completeness, y = Contamination)) +
  geom_point() +
  scale_y_reverse() +  # Invert the y-axis
  labs(title = "Genome Completeness vs. Contamination",
       x = "Completeness (%)",
       y = "Contamination (%)") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +  # horizontal line at contamination = 5%
  geom_vline(xintercept = 90, linetype = "dashed", color = "red")   # vertical line at completeness = 90%

#bad genomes (completeness < 90, contamination > 5)
bad_checkm <- checkm_out[checkm_out$Completeness < 90 | checkm_out$Contamination > 5, ][,1]
which(checkm_out$Completeness < 90 | checkm_out$Contamination > 5) # which genomes have low completeness or high contamination
which(checkm_out$Contamination > 5)
checkm_trimmed <- checkm_out[checkm_out$Completeness >= 90 & checkm_out$Contamination <= 5, ]

#GTDB-Tk-------------------------------------------------------------------------
gtdbtk_bac120_summary <- read_delim("~/Uni/11.masters_thesis/quality_assessment/gtdbtk_out/gtdbtk.bac120.summary.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

gtdbtk_bac120_summary[gtdbtk_bac120_summary=="N/A"] <- NA #replace all "N/A" with NA
gtdbtk_trimmed <- gtdbtk_bac120_summary[ , colSums(is.na(gtdbtk_bac120_summary)) < nrow(gtdbtk_bac120_summary)] #remove all empty columns

unique(gtdbtk_trimmed$classification)
#all genomes have been classified as Strep mutans

unique(gtdbtk_trimmed$note) #all genomes were classified based on ANI alone

#QUAST---------------------------------------------------------------------------

quast_results <- read_delim("~/Uni/11.masters_thesis/quality_assessment/quast_results.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

#remove # from colnames
colnames(quast_results) <- lapply(colnames(quast_results), 
                                  gsub, pattern='# ', replacement='')

#Total number of contigs per genome
ggplot(quast_results, aes(x=contigs)) + geom_histogram(col = "black", fill = "black",
                                                       binwidth = 5) +
  labs(title = "Total number of contigs per genome",
       x = "",
       y = "Count") +
  theme_light() +  # Set white background
  theme(plot.title = element_text(hjust = 0.5)) +   #title in the middle
  scale_x_continuous(expand = c(0, 0), limits = c(0, 470)) +    #0,0 in the bottom left
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5))

#Genome size
ggplot(quast_results, aes(x=`Total length`)) + geom_histogram(col = "black", fill = "black",
                                                              binwidth = 5000) +
  labs(title = "Genome size",
       x = "",
       y = "Count") +
  theme_light() +  # Set white background
  theme(plot.title = element_text(hjust = 0.5)) +   #title in the middle
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2500000)) +    #0,0 in the bottom left
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 1500000, linetype = "dashed", color = "red")

which(quast_results$`Total length` < 1500000)
quast_results[473,1]

#N50
ggplot(quast_results, aes(x=N50)) + geom_histogram(col = "black", fill = "black",
                                                   binwidth = 10000) +
  labs(title = "N50",
       x = "",
       y = "Count") +
  theme_light() +  # Set white background
  theme(plot.title = element_text(hjust = 0.5)) +   #title in the middle
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2500000)) +    #0,0 in the bottom left
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) +
  geom_vline(xintercept = 20000, linetype = "dashed", color = "red") +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5))

sum(quast_results$N50 < 100000) #191
sum(quast_results$N50 < 20000)  #6
sum(quast_results$N50 < 5000)   #0

# combined plot
ggplot(quast_results, aes(x = `Total length`, y = N50)) + geom_point() +
  labs(title = "Genome size vs. N50",
       x = "Genome size (bp)",
       y = "N50") +
  theme_light() +
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=24),
        axis.title=element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  geom_hline(yintercept = 20000, linetype = "dashed", color = "red") +  # horizontal line at contamination = 5%
  geom_vline(xintercept = 1500000, linetype = "dashed", color = "red")   # vertical line at completeness = 90%

#bad genomes (total length < 1.5M, N50 < 20k)
bad_quast <- quast_results[quast_results$`Total length` < 1500000 | quast_results$N50 < 20000,]

#bad genomes---------------------------------------------------------------------

bad_quast_names <- pull(bad_quast, Assembly)
bad_checkm_names <- pull(bad_checkm, `Bin Id`)

bad_quast_names %in% bad_checkm_names # which genomes have been qualified as 'bad' with both methods

all_bad <- unique(c(bad_checkm_names, bad_quast_names))








