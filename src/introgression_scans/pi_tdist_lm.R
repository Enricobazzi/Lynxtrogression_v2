library(tidyverse)
library(lme4)
library(DHARMa)
library("visreg")

get_ran_df <- function(intro){
  ran_pi <- data.frame()
  for (n in 0:99){
    t <- read.table(paste0("data/introgression_scans/", intro, "_randomwins/random_windows.", n, ".diversity.bed"),
                    col.names = c("chrom", "start", "end", "pi")) %>% 
      mutate(n = {{n}})
    t <- t %>% left_join(read.table(paste0("data/introgression_scans/", intro, "_randomwins/random_windows.",
                                           n, ".t_dist.bed"),
                                    col.names = c("chrom", "start", "end", "t_dist")))
    ran_pi <- rbind(ran_pi, t)
  }
  ran_pi$treat <- "ran"
  return(ran_pi)
}

get_intro_df <- function(intro){
  obs_pi <- read.table(paste0("data/introgression_scans/bed_files/", intro, "_intro.diversity.bed"),
                       col.names = c("chrom", "start", "end", "pi"))
  obs_pi$t_dist <- read.table(paste0("data/introgression_scans/bed_files/", intro, "_intro.merged.t_dist.bed"),
                              col.names = c("chrom", "start", "end", "t_dist"))[,4]
  obs_pi$n <- -99
  obs_pi$treat <- "obs"
  return(obs_pi)
}


lpa_ran <- get_ran_df("wel_and_sel_to_lpa")
lpa_intro <- get_intro_df("wel_and_sel_to_lpa")
lpa_ran_pi <- exp(mean(log(lpa_ran$pi)))
lpa_intro_pi <- exp(mean(log(lpa_intro$pi)))
lpa_ran_t_dist <- exp(mean(log(lpa_ran$t_dist)))
lpa_intro_t_dist <- exp(mean(log(lpa_intro$t_dist)))

wel_ran <- get_ran_df("lpa_to_wel")
wel_intro <- get_intro_df("lpa_to_wel")
wel_ran_pi <- exp(mean(log(wel_ran$pi)))
wel_intro_pi <- exp(mean(log(wel_intro$pi)))
wel_ran_t_dist <- exp(mean(log(wel_ran$t_dist)))
wel_intro_t_dist <- exp(mean(log(wel_intro$t_dist)))

sel_ran <- get_ran_df("lpa_to_sel")
sel_intro <- get_intro_df("lpa_to_sel")
sel_ran_pi <- exp(mean(log(sel_ran$pi)))
sel_intro_pi <- exp(mean(log(sel_intro$pi)))
sel_ran_t_dist <- exp(mean(log(sel_ran$t_dist)))
sel_intro_t_dist <- exp(mean(log(sel_intro$t_dist)))

# before model 
print(paste("lpa mean pi:", lpa_ran_pi))
print(paste("lpa intro pi:", lpa_intro_pi))
print(paste("lpa pi in intro is", lpa_intro_pi/lpa_ran_pi, "times pi in random"))
print(paste("wel mean pi:", wel_ran_pi))
print(paste("wel intro pi:", wel_intro_pi))
print(paste("wel pi in intro is", wel_intro_pi/wel_ran_pi, "times pi in random"))
print(paste("sel mean pi:", sel_ran_pi))
print(paste("sel intro pi:", sel_intro_pi))
print(paste("sel pi in intro is", sel_intro_pi/sel_ran_pi, "times pi in random"))

print(paste("lpa mean t_dist:", lpa_ran_t_dist))
print(paste("lpa intro t_dist:", lpa_intro_t_dist))
print(paste("wel mean t_dist:", wel_ran_t_dist))
print(paste("wel intro t_dist:", wel_intro_t_dist))
print(paste("sel mean t_dist:", sel_ran_t_dist))
print(paste("sel intro t_dist:", sel_intro_t_dist))

