library(tidyverse)
library(extrafont)

intros <- c("wel_and_sel_to_lpa", "lpa_to_wel", "lpa_to_sel")
intros_long <- c("Iberian lynx", "Western Eurasian lynx", "Southern Eurasian lynx")

get_obs_stat <- function(intro, stat){
  return(read.table(paste0("data/introgression_scans/bed_files/", intro, "_intro.merged.", stat, ".bed"),
             col.names = c("chrom", "start", "end", "pi"))[,4])
}

get_ran_stat <- function(intro, n, stat){
  return(read.table(paste0("data/introgression_scans/", intro, "_randomwins/random_windows.", n, ".", stat, ".bed"),
             col.names = c("chrom", "start", "end", "pi"))[,4])
}

for (stat in c("diversity", "t_dist", "gene_overlap", "cds_overlap", "n_genes")){
  if (stat == "diversity"){
    yaxname <- "Nucleotide diversity (π)"
  } else if (stat == "t_dist"){
    yaxname <- "Distance to telomere (bp)"
  } else if (stat == "gene_overlap"){
    yaxname <- "Overlap with genes (bp)"
  } else if (stat == "cds_overlap"){
    yaxname <- "Overlap with CDS (bp)"
  } else if (stat == "n_genes"){
    yaxname <- "Number of genes"
  }
  
  df <- data.frame()
    for (i in 1:3){
      obs <- data.frame(
        pop = intros_long[i],
        pop_n = i,
        treat = "introgressed",
        stat = get_obs_stat(intros[i], stat)
      )
      df <- rbind(df, obs)
      for (n in 0:99){
        ran <- data.frame(
          pop = intros_long[i],
          pop_n = i,
          treat = "random",
          stat = get_ran_stat(intros[i], n, stat)
        )
        df <- rbind(df, ran)
      }
    }
    
    # violins  '#12a0d7'
    violins <- ggplot() +
      geom_violin(data=df, aes(x=pop, y=stat, fill=treat),
                  position=position_dodge(0.7), alpha = 0.7, lwd = 0.4) +
      geom_boxplot(data=df, aes(x=pop, y=stat, fill=treat),
                   position=position_dodge(0.7), width=0.2, outlier.shape = NA, lwd = 0.3) +
      scale_fill_manual(values = c('#e3aa00', '#003c82')) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(y=yaxname, x="") +
      theme_minimal(base_size = 12, base_family = "sans") +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray80", size = 0.5),
        panel.grid.minor = element_line(color = "gray90", size = 0.25),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA)
      )
    ggsave(paste0("plots/introgression_scans/obs_v_rand/", stat, ".violins.pdf"),
           plot = violins, width = 8, height = 4)
}

#################

for (i in 1:3){
  obs <- data.frame(
    pop = intros_long[i],
    pop_n = i,
    treat = "introgressed",
    stat = get_obs_stat(intros[i], stat)
  )
  statt <- c()
  for (n in 0:99){
    ran <- data.frame(
      pop = intros_long[i],
      pop_n = i,
      treat = "random",
      stat = get_ran_stat(intros[i], n, stat)
    )
    statt <- c(statt, max(ran$stat))
  }
  hist(statt, breaks = 20, main = paste0("histogram of ", intros_long[i]))
  abline(v=max(obs$stat))
}
