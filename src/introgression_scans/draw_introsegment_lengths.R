# draw distribution of length of introgressed segments in ab and ba directions

# the script will look for the files:
# - ab_introgressed_095.merged.bed
# - ba_introgressed_095.merged.bed
# inside of the folder specified as first positional arguments

# the second positional argument will define the output folder for the plot

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
ifolder <- args[1]
ofolder <- args[2]

ab_bed <- paste0(ifolder, "/ab_introgressed_095.merged.bed")
ba_bed <- paste0(ifolder, "/ba_introgressed_095.merged.bed")

ab <- read.table(ab_bed, col.names = c("chrom", "start", "end"))
ba <- read.table(ba_bed, col.names = c("chrom", "start", "end"))

ab$len <- ab$end - ab$start
ba$len <- ba$end - ba$start

hist <- ggplot() +
  geom_histogram(data = ab, aes(x = len / 1000, y = ..density..),
                 bins = 100, fill = "aquamarine", color = "black", alpha = .7) +
  geom_histogram(data = ba, aes(x = len / 1000, y = ..density..),
                 bins = 100, fill = "blue", color = "black", alpha = .7) +
  ggtitle(paste("Introgressed segment length")) +
  xlab("Thousands of base pairs") +
  ylab("Frequency") +
  theme_minimal()

ggsave(
  filename = paste0(ofolder, "/introsegment_lengths.pdf"),
  plot = hist, height = 4, width = 6
)
