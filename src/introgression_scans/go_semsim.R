library(GOSemSim)
library(org.Hs.eg.db)

hsGO <- godata('org.Hs.eg.db', ont = "BP")
intros <- c("lpa_to_sel", "lpa_to_wel", "wel_and_sel_to_lpa")

for (intro in intros){
    go_table <- read.table(paste0("data/introgression_scans/genes/", intro, ".overrepresented_GOs.tsv"),
                           sep = '\t', header = T, as.is = T)
    go_terms <- go_table$GO.ID
    go_labels <- paste(go_table$GO.ID, go_table$Term, sep = " - ")
    print(paste(intro, go_labels, sep = " - "))
    
    simMatrix <- mgoSim(go_terms, go_terms, semData = hsGO, measure = "Wang", combine = NULL)
    dist_matrix <- as.dist(1 - simMatrix)
    
    hc <- hclust(dist_matrix, method = "complete") # add labels for plot
    membs1 <- cutree(hc, k=8)
    hc$labels <- paste(as.numeric(membs1), go_labels, sep = " - ")

    pdf(paste0("plots/introgression_scans/", intro, ".go_dendrogram.over.pdf"))
    plot(hc, main = "Dendrogram from Similarity Matrix", xlab = "", sub = "", cex = 0.6)
    dev.off()
}

for (intro in intros){
  go_table <- read.table(paste0("data/introgression_scans/genes/", intro, ".underrepresented_GOs.tsv"),
                         sep = '\t', header = T, as.is = T)
  go_terms <- go_table$GO.ID
  go_labels <- paste(go_table$GO.ID, go_table$Term, sep = " - ")
  print(paste(intro, go_labels, sep = " - "))

  simMatrix <- mgoSim(go_terms, go_terms, semData = hsGO, measure = "Wang", combine = NULL)
  dist_matrix <- as.dist(1 - simMatrix)
  
  hc <- hclust(dist_matrix, method = "complete") # add labels for plot
  membs1 <- cutree(hc, k=8)
  hc$labels <- paste(as.numeric(membs1), go_labels, sep = " - ")
  
  pdf(paste0("plots/introgression_scans/", intro, ".go_dendrogram.under.pdf"))
  plot(hc, main = "Dendrogram from Similarity Matrix", xlab = "", sub = "", cex = 0.4)
  dev.off()
}
