library(ggplot2)
library(dplyr)

# read the missingness table
miss_table <- read.table("data/variant_filtering/missing/miss_table.txt",
                         header = TRUE)
miss_table$percent_out <- round(miss_table$n_filter/miss_table$n_tot*100, 2)

miss_plot <- ggplot(miss_table) +
  geom_line(data = miss_table %>% filter(population == "lpa"),
            aes(x = f_miss, y = percent_out), size = 0.5, color = "grey") +
  geom_line(data = miss_table %>% filter(population == "eel"),
            aes(x = f_miss, y = percent_out), size = 0.5, color = "grey") +
  geom_line(data = miss_table %>% filter(population == "wel"),
            aes(x = f_miss, y = percent_out), size = 0.5, color = "grey") +
  geom_line(data = miss_table %>% filter(population == "sel"),
            aes(x = f_miss, y = percent_out), size = 0.5, color = "grey") +
  geom_point(data = miss_table, 
             aes(x = f_miss, y = percent_out, fill = population),
             shape = 21, size = 2.5) +
  theme(panel.grid.major.y = element_line(color = "black", size = 0.25, linetype = 2)) +
  scale_y_continuous(n.breaks = 12) +
  xlab(paste0("Maximum proportion of missing genotypes")) +
  ylab(paste0("Percentage of SNPs filtered"))

ggsave(filename = "data/variant_filtering/missing/missing_plot.png",
       plot = miss_plot, width = 4, height = 3)
