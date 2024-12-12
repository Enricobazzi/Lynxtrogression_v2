library(tidyverse)
library(lsr)

get_intro_df <- function(intro){
  obs_pi <- read.table(paste0("data/introgression_scans/bed_files/", intro, "_intro.merged.diversity.bed"),
                       col.names = c("chrom", "start", "end", "pi"))
  obs_pi$t_dist <- read.table(paste0("data/introgression_scans/bed_files/", intro, "_intro.merged.t_dist.bed"),
                              col.names = c("chrom", "start", "end", "t_dist"))[,4]
  obs_pi$gol <- get_gol(read.table(paste0("data/introgression_scans/bed_files/", intro, "_intro.merged.gene_overlap.bed"),
                        col.names = c("chrom", "start", "end", "gol")))[,4]
  obs_pi$n <- -99
  obs_pi$treat <- "obs"
  return(obs_pi)
}

get_ran_df <- function(intro, n){
  ran_pi <- data.frame()
  t <- read.table(paste0("data/introgression_scans/", intro, "_randomwins/random_windows.", n, ".diversity.bed"),
                  col.names = c("chrom", "start", "end", "pi")) %>% 
    mutate(n = {{n}})
  t <- t %>% left_join(get_gol(read.table(paste0("data/introgression_scans/", intro, "_randomwins/random_windows.",
                                                 n, ".gene_overlap.bed"), col.names = c("chrom", "start", "end", "gol"))),
                       by = c("chrom", "start", "end"))
  t <- t %>% left_join(read.table(paste0("data/introgression_scans/", intro, "_randomwins/random_windows.",
                                         n, ".t_dist.bed"), col.names = c("chrom", "start", "end", "t_dist")),
                       by = c("chrom", "start", "end"))
  ran_pi <- rbind(ran_pi, t)
  ran_pi$treat <- "ran"
  return(ran_pi)
}

get_gol <- function(a){
  b <- data.frame()
  for (i in 1:nrow(a)) {
    if (i == 1){
      b <- rbind(b, a[i,])
    } else {
      if (a[i,]$chrom == b[nrow(b),]$chr & a[i,]$start == b[nrow(b),]$start & a[i,]$end == b[nrow(b),]$end) {
        b[nrow(b),]$gol <- b[nrow(b),]$gol + a[i,]$gol
      } else {
        b <- rbind(b, a[i,])
      }
    }
  }
  return(b)
}

get_mean_ci <- function(x){
  # intercept(mean), 2.5%, 97.5%
  lmx <- lm(x ~ 1)
  m <- as.numeric(lmx$coefficients)
  lci <- as.numeric(confint(lmx, level=0.95)[1])
  hci <- as.numeric(confint(lmx, level=0.95)[2])
  return(c(m, lci, hci))
}

print_stuff <- function(intro, string, mean_ci){
  print(paste0(intro,": ", string, " = ", mean_ci[1], " [", mean_ci[2], " - ", mean_ci[3], "]"))
}


intros <- c("wel_and_sel_to_lpa", "lpa_to_wel", "lpa_to_sel")

for (intro_file in intros){
  
  pi_median_diff <- vector(mode = "list", length = 100)
  pi_prop <- vector(mode = "list", length = 100)
  tdist_median_diff <- vector(mode = "list", length = 100)
  tdist_prop <- vector(mode = "list", length = 100)
  gol_median_diff <- vector(mode = "list", length = 100)
  gol_prop <- vector(mode = "list", length = 100)
  pi_cds <- vector(mode = "list", length = 100)
  tdist_cds <- vector(mode = "list", length = 100)
  gol_cds <- vector(mode = "list", length = 100)
  
  for (n in 0:99){
    cat("\r", n)
    intro <- get_intro_df(intro_file)
    ran <- get_ran_df(intro_file, n)
    
    pi_median_diff[[n+1]] <- median(sample(intro$pi, length(ran$pi)) - ran$pi)
    pi_prop[[n+1]] <- median(sample(intro$pi, length(ran$pi)) / ran$pi)
    tdist_median_diff[[n+1]] <- median(sample(intro$t_dist, length(ran$t_dist)) - ran$t_dist)
    tdist_prop[[n+1]] <- median(sample(intro$t_dist, length(ran$t_dist)) / ran$t_dist)
    gol_median_diff[[n+1]] <- mean(sample(intro$gol, length(ran$gol)) - ran$gol)
    pi_cds[[n+1]] <- cohensD(log(intro$pi), log(ran$pi))
    tdist_cds[[n+1]] <- cohensD(log(intro$t_dist), log(ran$t_dist))
    gol_cds[[n+1]] <- cohensD(log(intro$gol[intro$gol!=0]), log(ran$gol[ran$gol!=0]))
  }
  
  pi_median_diff <- unlist(pi_median_diff)
  pi_prop <- unlist(pi_prop)
  tdist_median_diff <- unlist(tdist_median_diff)
  tdist_prop <- unlist(tdist_prop)
  gol_median_diff <- unlist(gol_median_diff)
  gol_prop <- unlist(gol_prop)
  pi_cds <- unlist(pi_cds)
  tdist_cds <- unlist(tdist_cds)
  gol_cds <- unlist(gol_cds)
  
  print_stuff(intro_file, "pi_median_diff", get_mean_ci(pi_median_diff))
  print_stuff(intro_file, "pi_prop", get_mean_ci(pi_prop))
  print_stuff(intro_file, "pi_cds", get_mean_ci(pi_cds))
  print_stuff(intro_file, "tdist_median_diff", get_mean_ci(tdist_median_diff))
  print_stuff(intro_file, "tdist_prop", get_mean_ci(tdist_prop))
  print_stuff(intro_file, "tdist_cds", get_mean_ci(tdist_cds))
  print_stuff(intro_file, "gol_median_diff", get_mean_ci(gol_median_diff))
  # print_stuff(intro_file, "gol_prop", get_mean_ci(gol_prop))
  print_stuff(intro_file, "gol_cds", get_mean_ci(gol_cds))
}

