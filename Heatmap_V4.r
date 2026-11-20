R
# This is the heatmap one 

# Very notably a work in Progress


# Ensure all are downloaded already 
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(sesame)
library(data.table)
library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

sesameDataCache("idatSignature")

# Load the custom manifest as a dataframe 
# !!!!!!!!!!!!!!!Replace with proper manifest!!!!!!!!!!!!!!!!!
excel_path <- "CUSTOM_MANIFEST.xlsx"
anno_df <- as.data.frame(read_excel(excel_path, sheet = 1))

# ERROR CHECK
head(anno_df)

anno_clean <- anno_df[!is.na(anno_df$CHR) & !is.na(anno_df$MAPINFO),]

# Clean the probe ID (removes suffix if present)
# (Many probes have _BC11 or _TC11)
anno_clean$Core_Probe <- sub("_.*", "", anno_clean$Probe_ID)


# Create the GenomicRanges 
gr <- makeGRangesFromDataFrame(
  df = anno_clean,
  keep.extra.columns = TRUE,
  seqnames.field = "CHR",       
  start.field = "MAPINFO",      
  end.field = "MAPINFO"         
)

# 2. Set up IDAT paths
# Replace with proper .IDAT files directory
data_dir <- "file/path/idats"
basenames <- list.files(data_dir, pattern = "_Red\\.idat$", recursive = TRUE, full.names = TRUE)
basenames <- sub("_Red\\.idat$", "", basenames)

# Read IDAT pair using manifest you generated
beta_list <- lapply(basenames, function(bn) {
  sdf <- readIDATpair(bn, manifest = anno_df)
  openSesame(sdf)
})

# Generate the Betas matrix
names(beta_list) <- basename(basenames)
betas <- do.call(cbind, beta_list)

# ERROR CHECK
dim(betas)

# Get the gene features for overlap 
# We use h19 i believe because we are basing this off epic V1
gene_features <- build_annotations(genome = 'hg19', annotations = c('hg19_basicgenes'))
overlap_hits <- findOverlaps(gr, gene_features, ignore.strand = TRUE)
intragenic_indices <- unique(queryHits(overlap_hits))

# Also run annotatr to get the overlapping symbol names for the conditional step
ann_regions <- annotate_regions(regions = gr, annotations = gene_features, ignore.strand = TRUE, quiet = TRUE)
annotated_df <- data.frame(ann_regions)

# 5. Get all gene ranges for distance calculation 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_gr <- genes(txdb)
nearest_indices <- nearest(gr, genes_gr)

# Map Entrez IDs of nearest genes to symbols
entrez_ids <- names(genes_gr)[nearest_indices]
nearest_symbols <- mapIds(org.Hs.eg.db, 
                          keys = entrez_ids, 
                          keytype = "ENTREZID", 
                          column = "SYMBOL", 
                          multiVals = "first")

# Calculate distance to the nearest gene
distances <- distance(gr, genes_gr[nearest_indices])

# Apply the conditional logic across all probes
anno_clean$Region_Type <- ifelse(seq_along(gr) %in% intragenic_indices, "Intragenic", "Intergenic")

# Create the final annotated output column based on the condition
anno_clean$Annotation_Result <- mapply(
  function(is_intra, symbol_overlap, near_sym, dist) {
    if (is_intra) {
      # If intragenic, show the gene symbol it overlaps
      return(paste0("Intragenic: ", symbol_overlap))
    } else {
      # If intergenic, show the nearest gene and distance
      # I have no idea if this is the best way of doing things 
      return(sprintf("Intergenic -> Nearest: %s (%d bp away)", near_sym, dist))
    }
  },
  anno_clean$Region_Type == "Intragenic",
  annotated_df$annot.symbol[match(anno_clean$Probe_ID, annotated_df$Probe_ID)], 
  nearest_symbols,
  distances
)

# View the results for the first 20 probes
# Change with whatever it doesnt matter
View(head(anno_clean[, c("Probe_ID", "MAPINFO", "Region_Type", "Annotation_Result")], 20))

# does betas matrix has a clean probe ID column for matching
betas_df <- data.frame(
  Core_Probe = sub("_.*", "", rownames(betas)), 
  Full_Probe_ID = rownames(betas), 
  betas, 
  check.names = FALSE
)

# Merge annotation table (anno_clean) with the betas table by Core_Probe
annotated_betas_master <- merge(
  anno_clean[, c("Core_Probe", "Region_Type", "Annotation_Result", "MAPINFO")], 
  betas_df, 
  by = "Core_Probe", 
  all.y = TRUE
)

# Replace with actual genes 
target_genes <- c("GENE1", "GENE2")

metadata_cols <- c("Core_Probe", "Full_Probe_ID", "Region_Type", "Annotation_Result", "MAPINFO")
sample_cols <- setdiff(colnames(annotated_betas_master), metadata_cols)
sample_col_indices <- which(colnames(annotated_betas_master) %in% sample_cols)

betas_matrix <- as.matrix(annotated_betas_master[, sample_col_indices])
class(betas_matrix) <- "numeric"
rownames(betas_matrix) <- annotated_betas_master$Full_Probe_ID

target_row_indices <- c()
target_labels <- c()

for (gene in target_genes) {
  # Search for the gene keyword in annotation results
  matches <- which(grepl(gene, annotated_betas_master$Annotation_Result, ignore.case = TRUE))
  
  if (length(matches) > 0) {
    # Takes the first matching probe for this gene 
    idx <- matches[1] 
    target_row_indices <- c(target_row_indices, idx)
    
    # Format a clean row label showing the gene and its classification
    lbl <- paste0(toupper(gene), " (", annotated_betas_master$Core_Probe[idx], "): ", annotated_betas_master$Annotation_Result[idx]) 
    target_labels <- c(target_labels, lbl)
  } else {
    message("Note: No probe found for target gene: ", gene)
  }
}

# Subset your matrix to include the chosen target genes
# edit the betas finally to plot the matrix 
if (length(target_row_indices) > 0) {
  matrix_to_plot_custom <- betas_matrix[target_row_indices, , drop = FALSE]
  matrix_to_plot_custom[is.na(matrix_to_plot_custom)] <- 0
  
  # Define sample groups for column annotation
  # This helps with the heatmap
  sample_names <- colnames(matrix_to_plot_custom)
  group_labels <- character(length(sample_names))
  group_labels[grepl("cntrl", sample_names, ignore.case = TRUE)] <- "Control"
  group_labels[grepl("treatment", sample_names, ignore.case = TRUE)] <- "Expanded"

  annotation_col = data.frame(
    Group = factor(group_labels)
  )
  rownames(annotation_col) <- sample_names


  pheatmap(
    matrix_to_plot_custom,
    annotation_col = annotation_col,
    show_rownames = TRUE,           
    labels_row = target_labels, 
    show_colnames = TRUE,          
    labels_col = group_labels,     
    cluster_cols = TRUE,
    cluster_rows = if(length(target_row_indices) > 2) TRUE else FALSE, 
    scale = "none",                 
    color = colorRampPalette(c("blue”, "white", "red"))(100),
    main = "Gene Methylation Signatures",
    fontsize_row = 10               
  )
} else {
  stop("None of the specified target genes were found in the master table.")
}

