# ENRICHMENT #
library(topGO)
library(tidyverse)

get_cat_long <- function(cat){
  category_dictionary <- c(
    "0" = "Unclassified", 
    "1" = "Neural and synaptic structure and signaling",
    "2" = "Immune and inflammatory response",
    "3" = "Metabolism and catabolism",
    "4" = "Gene regulation",
    "5" = "Tissue differentiation, development and structuring",
    "6" = "Sensory and stimulus perception"
  )
  return(as.character(category_dictionary[cat]))
}

get_cat_color <- function(cat){
  color_dictionary <- c(
    "0" = "#d9cab3",
    "1" = "#005f99",
    "2" = "#a67500",
    "3" = "#f4d35e",
    "4" = "#407f7f",
    "5" = "#d96459",
    "6" = "#9b5094"
  )
  return(as.character(color_dictionary[cat]))
}

# Custom method to perform Fisher's exact test for underrepresentation (less frequent GO terms)
if (!isGeneric("GOFisherUnder")) {
  setGeneric("GOFisherUnder", function(object) standardGeneric("GOFisherUnder"))
}

setMethod("GOFisherUnder", "classicCount", function(object) {
  contMat <- contTable(object)
  if (all(contMat == 0)) {
    p.value <- 1
  } else {
    p.value <- fisher.test(contMat, alternative = "less")$p.value
  }
  return(p.value)
})

intros <- c(
  "lpa_to_wel",
  "lpa_to_sel",
  "wel_and_sel_to_lpa"
)
# Read the gene-to-GO mappings from a .tsv file
gene2GO_df <- read.delim("data/introgression_scans/genes/LYRU2_2A.FA.genego_table.tsv", header = TRUE, sep = "\t")
# Convert the gene2GO data frame into a named list where each gene ID points to a vector of GO terms
gene2GO <- setNames(strsplit(gene2GO_df$go_terms, ";"), gene2GO_df$gene_id)
# Get the list of all unique genes from the gene2GO mapping
allGenesList <- unique(gene2GO_df$gene_id)

# go to category table
go_cat_table <- read.table("data/introgression_scans/genes/go-categories.tsv", col.names = c("GO.ID", "Category"))
go_cat_table$Category_long <- sapply(as.character(go_cat_table$Category), get_cat_long)
go_cat_table$Color <- sapply(as.character(go_cat_table$Category), get_cat_color)

for (intro in intros){
  # Read the list of candidate genes for the current gene set
  candidateGenes <- readLines(paste0("data/introgression_scans/genes/", intro, ".geneids.txt"))
  
  # Create a named vector (geneList) for all genes:
  # 1 for candidate genes, 0 for non-candidate genes
  geneList <- setNames(as.numeric(allGenesList %in% candidateGenes), allGenesList)
  
  # Create the topGOdata object for both overrepresentation and underrepresentation tests
  GOdata <- new(
    "topGOdata",
    ontology = "BP",  # Choose "BP", "MF", or "CC" as appropriate
    allGenes = geneList,  # The named vector of all genes with 1 for candidates, 0 for non-candidates
    geneSel = function(p) p == 1,  # Function to select candidate genes (where p == 1)
    annot = annFUN.gene2GO,  # Annotation function for gene2GO mappings
    gene2GO = gene2GO  # Your gene-to-GO mappings
  )
  
  # Perform Overrepresentation Test using the "classic" algorithm and default Fisher's exact test
  resultFisherOver <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # Set up a new test for underrepresentation using the custom GOFisherUnder function
  test.stat <- new("classicCount", testStatistic = GOFisherUnder, name = "Fisher's exact test for underrepresentation")
  # Perform Underrepresentation Test using the "classic" algorithm and custom GOFisherUnder function
  resultFisherUnder <- getSigGroups(GOdata, test.stat)
  
  #Create a table with significant results with a p-value adjusted by FDR.
  result_over <- GenTable(GOdata, Fisher=resultFisherOver, topNodes=resultFisherOver@geneData[2], numChar=1000) %>% 
    as_tibble() %>% filter(Fisher<0.01)
  result_over$Fold_enrichment <- result_over$Significant / result_over$Expected
  
  # print for export:
  print(paste0(result_over$GO.ID, " - ", result_over$Term))
  # Add category
  result_over <- left_join(result_over, go_cat_table)

  # write tables  
  write.table(result_over, paste0("data/introgression_scans/genes/", intro, ".overrepresented_GOs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result_over[,c(1,6)], paste0("data/introgression_scans/genes/", intro, ".go_pval.over.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  result_under <- GenTable(GOdata, Fisher=resultFisherUnder, topNodes=resultFisherUnder@geneData[2], numChar=1000) %>% 
    as_tibble() %>% filter(Fisher<0.01)
  result_under$Fold_enrichment <- result_under$Significant / result_under$Expected
  write.table(result_under, paste0("data/introgression_scans/genes/", intro, ".underrepresented_GOs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result_over[,c(1,6)], paste0("data/introgression_scans/genes/", intro, ".go_pval.under.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# plot pie chart of categories
for (intro in intros){
  result_over <- read.table(paste0("data/introgression_scans/genes/", intro, ".overrepresented_GOs.tsv"),
                            sep = '\t', header = T, comment.char = "")
  
  # build table and colors
  pie_table <- data.frame()
  for (n in 0:6){
    # table
    group=get_cat_long(as.character(n))
    value=nrow(result_over[result_over$Category_long == group,])
    dfx <- data.frame(N=n, value=value, group=group)
    pie_table <- rbind(pie_table, dfx)
    }

  pie_table$group <- factor(pie_table$group, levels = c(
    "Unclassified", 
    "Neural and synaptic structure and signaling",
    "Immune and inflammatory response",
    "Metabolism and catabolism",
    "Gene regulation",
    "Tissue differentiation, development and structuring",
    "Sensory and stimulus perception"
    ))
  pie_table$color <- get_cat_color(as.character(pie_table$N))
  pie_table <- pie_table[pie_table$value > 0,]
  
  # pie-chart of Category
  pie <- ggplot(pie_table,
                aes(x = "", y = value, fill = group)) +
    geom_col(color = "black") +
    geom_text(aes(label = paste0(round(value/sum(value)*100, 1), "%")),
              position = position_stack(vjust = 0.5), size = 5) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = pie_table$color) +
    theme_void() #+ theme(legend.position = "none")
ggsave(paste0("plots/introgression_scans/genes/legend.pdf"), pie)
ggsave(paste0("plots/introgression_scans/genes/", intro, ".cat_pie.pdf"), pie,
       height = 3, width = 3)
  
}


