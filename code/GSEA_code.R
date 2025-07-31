# ==============================
# Batch GSEA Analysis Script in R
# ==============================
options(digits = 10)  # Set global display precision to 10 digits

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(ggplot2)
library(ggplotify)

# ==============================
# Step 1: Set Up the Environment
# ==============================

# Define directories
gene_list_file <- "/data/users/tfang/GSEA_ePC/geneList_forGSEA/GeneList_ID.csv"
gene_set_directory <- "/data/users/tfang/GSEA_ePC/geneList_forGSEA/Gene_1000days/data/data/GeneSet"
output_directory <- "/data/users/tfang/GSEA_ePC/geneList_forGSEA/GSEA_Results"
#dir.create(output_directory, showWarnings = FALSE)

# ==============================
# Step 2: Read and Prepare Gene List
# ==============================

gene_data <- read.csv(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
colnames(gene_data) <- c("gene", "weight")

gene_data$gene <- as.character(gene_data$gene)

gene_data <- gene_data[!is.na(gene_data$gene), ]
gene_data <- gene_data[order(gene_data$gene, -gene_data$weight), ]
gene_data <- gene_data[!duplicated(gene_data$gene), ]
gene_data$weight <- round(-gene_data$weight, digits = 10)

geneList <- gene_data$weight
names(geneList) <- as.character(gene_data$gene)
geneList <- sort(geneList, decreasing = TRUE)

# ==============================
# Step 3: Load and Process Gene Sets
# ==============================

gene_set_files <- list.files(gene_set_directory, pattern = "\\.txt$", full.names = TRUE)
gene_sets <- list()

convert_to_entrez <- function(genes) {
  converted <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(converted) == 0) {
    return(character(0))  # Avoiding NA
  }
  return(converted$ENTREZID)
}
TERM2GENE <- data.frame(TERM = character(), GENE = character(), stringsAsFactors = FALSE)

for (file in gene_set_files) {
  gene_set_name <- tools::file_path_sans_ext(basename(file))  # get the name of the gene set 
  genes <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  entrez_genes <- convert_to_entrez(genes)
  
  if (length(entrez_genes) > 0) {  
    TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = gene_set_name, GENE = entrez_genes))
  }
}

for (file in gene_set_files) {
  gene_set_name <- tools::file_path_sans_ext(basename(file))
  genes <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  entrez_genes <- convert_to_entrez(genes)
  gene_sets[[gene_set_name]] <- entrez_genes
}



# ==============================
# Step 4: Perform GSEA Analysis
# ==============================

gsea_results <- GSEA(
  geneList = geneList,
  TERM2GENE = TERM2GENE,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.5
)

# ==============================
# Step 4: Alternatively perform GSEA using fgsea
# ==============================

fgseaRes <- fgsea(
  pathways = gene_sets,  
  stats = geneList,      
  minSize = 10,         
  maxSize = 500,         
  nperm = 10000          
)


# ==============================
# Step 5: Visualize and Save Results
# ==============================
library(ggplotify)

# Set the option to avoid scientific notation
options(scipen = 999)

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

for (i in 1:nrow(gsea_results@result)) {
  gene_set_name <- fgseaRes$pathway[i]
  NES_value <- fgseaRes$NES[i]
  # Extract p-value and round it to 4 decimal places
  p_value <- fgseaRes$padj[i]
  
  
  gene_set_index <- which(gsea_results@result$ID == gene_set_name)
  if (length(gene_set_index) == 0 || is.na(NES_value) || is.na(p_value)) {
    next  
  }
  
  # generate GSEA plot
  gsea_plot <- gseaplot2(gsea_results, geneSetID = gene_set_index)
  
  
  if (!inherits(gsea_plot, "ggplot")) {
    gsea_plot <- as.ggplot(gsea_plot)
  }
  
  # add title
  gsea_plot <- gsea_plot + 
    ggtitle(paste(gene_set_name, "\nNES:", round(NES_value, 3), " P-value:", round(p_value, 4))) +
    theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 18, color = "darkblue", hjust = 0.5)
    )
  
  # save images
  ggsave(
    filename = file.path(output_directory, paste0(gene_set_name, "_GSEA_plot.svg")),
    plot = gsea_plot,
    width = 10,
    height = 6,
    dpi = 150,
    device = "svg"
  )
}

# Save results to a CSV file
fgseaRes[] <- lapply(fgseaRes, function(x) if (is.list(x)) sapply(x, toString) else x)
write.csv(fgseaRes, file.path(output_directory, "GSEA_results_fgsea.csv"), row.names = FALSE)
