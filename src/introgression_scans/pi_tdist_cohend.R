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
  pval_pi <- vector(mode = "list", length = 100)
  pval_tdist <- vector(mode = "list", length = 100)
  pval_gol <- vector(mode = "list", length = 100)
  
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
    # two tailed t-tests
    pval_pi[[n+1]] <- t.test(log(intro$pi), log(ran$pi))$p.value
    pval_tdist[[n+1]] <- t.test(log(intro$t_dist), log(ran$t_dist))$p.value
    pval_gol[[n+1]] <- t.test(log(intro$gol[intro$gol!=0]), log(ran$gol[ran$gol!=0]))$p.value
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
  pval_pi <- unlist(pval_pi)
  pval_tdist <- unlist(pval_tdist)
  pval_gol <- unlist(pval_gol)
  
  print_stuff(intro_file, "pi_median_diff", get_mean_ci(pi_median_diff))
  print_stuff(intro_file, "pi_prop", get_mean_ci(pi_prop))
  print_stuff(intro_file, "pi_cds", get_mean_ci(pi_cds))
  print(paste0(intro_file, "pi_pval = ", mean(pval_pi)))
  print_stuff(intro_file, "tdist_median_diff", get_mean_ci(tdist_median_diff))
  print_stuff(intro_file, "tdist_prop", get_mean_ci(tdist_prop))
  print_stuff(intro_file, "tdist_cds", get_mean_ci(tdist_cds))
  print(paste0(intro_file, "tdist_pval = ", mean(pval_tdist)))
  print_stuff(intro_file, "gol_median_diff", get_mean_ci(gol_median_diff))
  # print_stuff(intro_file, "gol_prop", get_mean_ci(gol_prop))
  print_stuff(intro_file, "gol_cds", get_mean_ci(gol_cds))
  print(paste0(intro_file, "gol_pval = ", mean(pval_gol)))
}

# 99[1] "wel_and_sel_to_lpa: pi_median_diff = 0.0003035261147395 [0.000302612962517584 - 0.000304439266961416]"
# [1] "wel_and_sel_to_lpa: pi_prop = 3.5097155650513 [3.48968356171155 - 3.52974756839106]"
# [1] "wel_and_sel_to_lpa: pi_cds = 1.81473702751558 [1.80759428336509 - 1.82187977166607]"
# [1] "wel_and_sel_to_lpapi_pval = 0"
# [1] "wel_and_sel_to_lpa: tdist_median_diff = -16677914.3525 [-16877817.3462318 - -16478011.3587682]"
# [1] "wel_and_sel_to_lpa: tdist_prop = 0.408292862165403 [0.404925934493006 - 0.411659789837799]"
# [1] "wel_and_sel_to_lpa: tdist_cds = 0.764682881689943 [0.756543346343594 - 0.772822417036293]"
# [1] "wel_and_sel_to_lpatdist_pval = 4.35875450152468e-72"
# [1] "wel_and_sel_to_lpa: gol_median_diff = -355.795306087785 [-613.606439337615 - -97.9841728379548]"
# [1] "wel_and_sel_to_lpa: gol_cds = 0.0988038939736316 [0.09191222676762 - 0.105695561179643]"
# [1] "wel_and_sel_to_lpagol_pval = 0.0845535218242214"
# 99[1] "lpa_to_wel: pi_median_diff = 0.0003463469642526 [0.000345014152302161 - 0.000347679776203039]"
# [1] "lpa_to_wel: pi_prop = 2.64885400983011 [2.63530897882956 - 2.66239904083066]"
# [1] "lpa_to_wel: pi_cds = 1.57546786381645 [1.56894360713591 - 1.58199212049698]"
# [1] "lpa_to_welpi_pval = 2.2548282924949e-286"
# [1] "lpa_to_wel: tdist_median_diff = -13090784.9 [-13279268.2081667 - -12902301.5918333]"
# [1] "lpa_to_wel: tdist_prop = 0.529066171906059 [0.525410476861176 - 0.532721866950941]"
# [1] "lpa_to_wel: tdist_cds = 0.562513168310943 [0.555879691905685 - 0.569146644716201]"
# [1] "lpa_to_weltdist_pval = 7.7897481072921e-43"
# [1] "lpa_to_wel: gol_median_diff = -1222.68480497812 [-1426.14042250334 - -1019.2291874529]"
# [1] "lpa_to_wel: gol_cds = 0.12304867301763 [0.117250756371846 - 0.128846589663415]"
# [1] "lpa_to_welgol_pval = 0.0192976966630414"
# 99[1] "lpa_to_sel: pi_median_diff = 0.00043813473276735 [0.000435431669534161 - 0.000440837796000539]"
# [1] "lpa_to_sel: pi_prop = 3.05411432065177 [3.02859086040853 - 3.07963778089501]"
# [1] "lpa_to_sel: pi_cds = 1.82628262780725 [1.81433271864141 - 1.83823253697309]"
# [1] "lpa_to_selpi_pval = 2.32625838625886e-122"
# [1] "lpa_to_sel: tdist_median_diff = -11257979.455 [-11565599.88507 - -10950359.02493]"
# [1] "lpa_to_sel: tdist_prop = 0.587537602975692 [0.580897253326475 - 0.594177952624908]"
# [1] "lpa_to_sel: tdist_cds = 0.488576743672081 [0.478272839356233 - 0.498880647987929]"
# [1] "lpa_to_seltdist_pval = 1.47229896184924e-09"
# [1] "lpa_to_sel: gol_median_diff = 712.032636679728 [348.429480850042 - 1075.63579250941]"
# [1] "lpa_to_sel: gol_cds = 0.0678479147911244 [0.0590078991441254 - 0.0766879304381233]"
# [1] "lpa_to_selgol_pval = 0.437028822625718"

