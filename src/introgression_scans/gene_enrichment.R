# ENRICHMENT #
library(topGO)
library(tidyverse)

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
gene2GO <- setNames(strsplit(gene2GO_df$go_terms, ","), gene2GO_df$gene_id)
# Get the list of all unique genes from the gene2GO mapping
allGenesList <- unique(gene2GO_df$gene_id)

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
  resultFisherOver <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Set up a new test for underrepresentation using the custom GOFisherUnder function
  test.stat <- new("classicCount", testStatistic = GOFisherUnder, name = "Fisher's exact test for underrepresentation")
  
  # Perform Underrepresentation Test using the "classic" algorithm and custom GOFisherUnder function
  resultFisherUnder <- getSigGroups(GOdata, test.stat)
  
  #Create a table with significant results with a p-value adjusted by FDR.
  result_over <- GenTable(GOdata, Fisher=resultFisherOver, topNodes=resultFisherOver@geneData[2], numChar=1000) %>% 
    as_tibble() %>% filter(Fisher<0.01)
  write.table(result_over, paste0("data/introgression_scans/genes/", intro, ".overrepresented_GOs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  result_under <- GenTable(GOdata, Fisher=resultFisherUnder, topNodes=resultFisherUnder@geneData[2], numChar=1000) %>% 
    as_tibble() %>% filter(Fisher<0.01)
  write.table(result_under, paste0("data/introgression_scans/genes/", intro, ".underrepresented_GOs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
}
