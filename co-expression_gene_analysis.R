library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(limma)
library(gridExtra)

allowWGCNAThreads()  

# 1. Fetch Data ------------------------------------------------
data <- read.delim('C:/Users/prade/Downloads/GSE86356_DMseq_sample_converted_nonUTR_psi.txt.gz', header = T)

# get metadata
geo_id <- "GSE86356"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,8,43)]

# prepare data
data[1:10,1:10]

# Expression data = data
# Phenotype data = phenodata
# Rename to match your intended variable name
phenodata <- phenoData

# Set sample IDs as rownames in phenodata
rownames(phenodata) <- phenodata$title

# Expression dataset columns (skip gene ID column)
exprSamples <- colnames(data)[-1]    
phenoSamples <- rownames(phenodata)  

# Find common samples
commonSamples <- intersect(exprSamples, phenoSamples)

# Subset expression data
data_inner <- data[, c("X.Event", commonSamples)]

# Subset phenotype data
phenodata_inner <- phenodata[commonSamples, ]

# Check alignment
all(colnames(data_inner)[-1] == rownames(phenodata_inner))
library(tidyverse)
library(tibble)

# If your data already has a column "X.Event", use a new name for rownames
data_long <- data %>%
  rownames_to_column(var = "gene_event") %>%  # changed column name
  pivot_longer(-gene_event, names_to = "samples", values_to = "counts")

# Inner join with phenodata based on actual sample ID
data_inner <- data_long %>%
  inner_join(phenoData, by = c("samples" = "title")) %>%
  select(gene_event, counts, geo_accession) %>%
  pivot_wider(names_from = geo_accession, values_from = counts) %>%
  column_to_rownames(var = "gene_event")


# Pivot expression data longer (no need for rownames_to_column since X.Event exists)
data_long <- data %>%
  pivot_longer(-X.Event, names_to = "samples", values_to = "counts")

# Inner join with phenodata on sample IDs and store in a new dataset
data_inner2 <- data_long %>%
  inner_join(phenoData, by = c("samples" = "title")) %>%
  select(X.Event, counts, geo_accession) %>%
  pivot_wider(names_from = geo_accession, values_from = counts) %>%
  column_to_rownames(var = "X.Event")

data_working <- data_inner

#removing n/a part
data_working[data_working == "n/a"] <- NA
data_working <- as.data.frame(lapply(data_working, as.numeric))
rownames(data_working) <- rownames(data_inner2)

#2 QC - outlier detection

gsg <- goodSamplesGenes(data_working, verbose = 3)
summary(gsg)

if (!gsg$allOK) {
  data_working <- data_working[gsg$goodSamples, gsg$goodGenes]
}
data_working[is.na(data_working)] <- 0


gsg <- goodSamplesGenes(t(data_working))
summary(gsg)
gsg$allOK
data_working1<- data_working
rownames(data_working1) <- seq_len(nrow(data_working1))
table(gsg$goodGenes)
table(gsg$goodSamples)
data_working1 <- data_working1[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data_working1)), method = "average")
plot(htree)

# pca - method 2

pca <- prcomp(t(data_working1))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# exclude outlier samples
samples.to.be.excluded <- c('GSM2309583', 'GSM2309536', 'GSM2309535' ,'GSM2309570' ,'GSM2309561')

# 3. Normalization and Data Preparation ----------------------------------------------------------------------
# create a deseq2 dataset

# First, exclude outlier samples from data_working1
data_working1 <- data_working1[, !colnames(data_working1) %in% samples.to.be.excluded]

# exclude outlier samples from phenoData
colData <- phenoData %>%
  filter(!row.names(.) %in% samples.to.be.excluded)

# fixing column names in colData
names(colData)
names(colData) <- gsub('_ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# Check dimensions before matching
cat("Dimensions before matching:\n")
cat("data_working1 columns:", ncol(data_working1), "\n")
cat("colData rows:", nrow(colData), "\n")

# Find common samples between data and metadata
commonSamples <- intersect(colnames(data_working1), rownames(colData))
cat("Common samples found:", length(commonSamples), "\n")

# Subset both datasets to common samples
data_working1 <- data_working1[, commonSamples]
colData <- colData[commonSamples, ]

# Reorder colData to match data_working1 column order
colData <- colData[colnames(data_working1), ]

# Check they match now
cat("\nDimensions after matching:\n")
cat("data_working1 columns:", ncol(data_working1), "\n")
cat("colData rows:", nrow(colData), "\n")
cat("Row names match:", all(colnames(data_working1) == rownames(colData)), "\n")

#making the dataset into integers
head(data_working1[, 1:5])

#dataset cleaning step for PSI data (instead of Dseq)
library(WGCNA)

# ----------------------------------------
# Replace NAs with 0 (already done, repeat safe)
data_working1[is.na(data_working1)] <- 0

# Filter low-variance events (in-place)
varEvents <- apply(data_working1, 1, var)
data_working1 <- data_working1[varEvents > median(varEvents), ]

# ============================================================================
# 5. WGCNA PREPARATION WITH LIMMA-FILTERED DATA
# ============================================================================

# Option 1: Use all filtered data for WGCNA
datExpr <- t(data_working1)

# Option 2: Use only significant events from limma (more stringent)
# significant_events <- rownames(results_limma)
# data_significant <- data_working1[significant_events, ]
# datExpr <- t(data_significant)

# Check for good samples/genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# CRITICAL: Reduce dataset size to avoid memory issues
# Select top variance features to make analysis computationally feasible
cat("\n=== Reducing dataset size for memory management ===\n")
cat("Original number of features:", ncol(datExpr), "\n")

# Calculate variance for each feature
feature_vars <- apply(datExpr, 2, var)

# Select top 5000 most variable features (adjust as needed)
# You can increase this if you have more RAM, or decrease if still getting errors
n_features <- min(5000, ncol(datExpr))
top_var_indices <- order(feature_vars, decreasing = TRUE)[1:n_features]
datExpr <- datExpr[, top_var_indices]

cat("Reduced to top", ncol(datExpr), "most variable features\n")
cat("This represents", round(ncol(datExpr)/length(feature_vars)*100, 1), "% of original features\n")

# Optional: scale/center features for WGCNA
datExpr <- scale(datExpr)

# Free up memory
rm(feature_vars, top_var_indices)
gc()

# ============================================================================
# 6. WGCNA ANALYSIS
# ============================================================================

# Sample clustering
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Choose soft-thresholding power
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale independence and mean connectivity
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")

# Select optimal power (usually where curve plateaus, R^2 > 0.8)
picked_power <- sft$powerEstimate
if(is.na(picked_power)) {
  picked_power <- 6
}
cat("\nSelected power:", picked_power, "\n")

# Construct network and detect modules
# Use tryCatch to handle potential errors
cat("\n=== Starting Network Construction ===\n")

net <- tryCatch({
  blockwiseModules(datExpr,
                   power = picked_power,
                   TOMType = "unsigned",
                   minModuleSize = 30,
                   reassignThreshold = 0,
                   mergeCutHeight = 0.25,
                   numericLabels = TRUE,
                   pamRespectsDendro = FALSE,
                   saveTOMs = TRUE,
                   saveTOMFileBase = "PSI_TOM",
                   verbose = 3,
                   corType = "pearson",
                   maxBlockSize = 5000,
                   networkType = "unsigned",
                   deepSplit = 2,
                   detectCutHeight = 0.995)
}, error = function(e) {
  cat("\n!!! Error in blockwiseModules detected !!!\n")
  cat("Error message:", e$message, "\n")
  cat("\n=== Switching to Manual One-Step Network Construction ===\n")
  
  # IMPORTANT: Explicitly use WGCNA's cor function
  # This is the fix for the masking issue
  cat("Using manual network construction to avoid cor() masking issue...\n")
  
  # Reduce dataset size if too large (optional)
  if(ncol(datExpr) > 10000) {
    cat("Dataset has", ncol(datExpr), "features. Consider using top variance features.\n")
    # Optionally subset to top variance features
    # datExpr <- datExpr[, order(apply(datExpr, 2, var), decreasing=TRUE)[1:10000]]
  }
  
  cat("\nStep 1: Calculating adjacency matrix...\n")
  adjacency <- adjacency(datExpr, power = picked_power, type = "unsigned")
  
  cat("Step 2: Calculating TOM similarity...\n")
  TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
  dissTOM <- 1 - TOM
  
  # Clean up large objects
  rm(adjacency)
  gc()
  
  cat("Step 3: Hierarchical clustering...\n")
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  cat("Step 4: Detecting modules using dynamic tree cut...\n")
  dynamicMods <- cutreeDynamic(dendro = geneTree,
                               distM = dissTOM,
                               deepSplit = 2,
                               pamRespectsDendro = FALSE,
                               minClusterSize = 30)
  
  cat("Found", length(unique(dynamicMods)), "initial modules\n")
  
  # Convert to colors
  dynamicColors <- labels2colors(dynamicMods)
  
  cat("Step 5: Calculating module eigengenes...\n")
  # Use WGCNA's moduleEigengenes which handles the correlation internally
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  
  cat("Step 6: Merging close modules...\n")
  # mergeCloseModules also uses internal correlation
  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
  
  cat("After merging:", length(unique(merge$colors)), "modules\n")
  
  # Create dendrogram for merged modules
  MEDiss <- 1 - WGCNA::cor(merge$newMEs, use = "pairwise.complete.obs")
  METree <- hclust(as.dist(MEDiss), method = "average")
  
  cat("\n=== Manual network construction completed successfully ===\n\n")
  
  # Create output similar to blockwiseModules
  list(colors = merge$colors,
       unmergedColors = dynamicColors,
       MEs = merge$newMEs,
       dendrograms = list(geneTree),
       blockGenes = list(1:ncol(datExpr)),
       TOMFiles = "Manual_TOM_calculation")
})

cat("\nNetwork construction finished!\n")
cat("Number of modules detected:", length(unique(net$colors)), "\n")

# Module colors
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# Plot dendrogram with module colors
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
# ============================================================================
# 7. RELATE MODULES TO TRAITS
# ============================================================================

# Define traits from colData
traits <- colData

# Check what columns you have
print("Available columns in colData:")
print(colnames(traits))
print(str(traits))

# Convert ALL columns to numeric (handles both numeric and categorical)
traits_numeric <- data.frame(lapply(traits, function(x) {
  if (is.numeric(x)) {
    return(x)
  } else {
    # Convert factors/characters to numeric
    return(as.numeric(factor(x)))
  }
}))

# Restore row names
rownames(traits_numeric) <- rownames(traits)

# Ensure sample order matches
traits_numeric <- traits_numeric[rownames(datExpr), , drop = FALSE]

# Check if we have any traits
print(paste("Number of traits to analyze:", ncol(traits_numeric)))
print(colnames(traits_numeric))

# Stop if no traits available
if (ncol(traits_numeric) == 0) {
  stop("No traits available! Check your colData columns.")
}

# Calculate module eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

# Correlate modules with traits
moduleTraitCor <- cor(MEs, traits_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Visualize module-trait relationships
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits_numeric),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# ============================================================================
# 8. EXPORT RESULTS
# ============================================================================

# Export module membership
module_df <- data.frame(
  event_id = colnames(datExpr),
  module_color = moduleColors
)
write.csv(module_df, "WGCNA_module_assignments.csv", row.names = FALSE)

# Calculate gene (event) significance and module membership
# For a specific trait of interest
trait_of_interest <- traits_numeric[, 1]  # Adjust column index

GS <- as.numeric(cor(datExpr, trait_of_interest, use = "p"))
GS_pvalue <- as.numeric(corPvalueStudent(GS, nrow(datExpr)))

MM <- cor(datExpr, MEs, use = "p")
MM_pvalue <- corPvalueStudent(MM, nrow(datExpr))

# Export gene significance and module membership
gs_mm_df <- data.frame(
  event_id = colnames(datExpr),
  module = moduleColors,
  GS = GS,
  GS_pvalue = GS_pvalue
)

# Add MM for each module
for(mod in unique(moduleColors)) {
  me_col <- paste0("ME", mod)
  if(me_col %in% colnames(MM)) {
    gs_mm_df[[paste0("MM_", mod)]] <- MM[, me_col]
    gs_mm_df[[paste0("MM_pvalue_", mod)]] <- MM_pvalue[, me_col]
  }
}

write.csv(gs_mm_df, "WGCNA_gene_significance_module_membership.csv",
          row.names = FALSE)

# ============================================================================
# 9. HUB EVENT IDENTIFICATION
# ============================================================================

# Identify hub events in each module
for(module in unique(moduleColors)) {
  if(module == "grey") next  # Skip unassigned events
  
  # Get events in this module
  module_events <- (moduleColors == module)
  
  # Get module membership values
  me_col <- paste0("ME", module)
  if(me_col %in% colnames(MM)) {
    mm_values <- MM[module_events, me_col]
    
    # Get top 10 hub events
    top_hubs <- names(sort(abs(mm_values), decreasing = TRUE)[1:10])
    
    cat("\nTop 10 hub events in", module, "module:\n")
    print(top_hubs)
    
    # Save to file
    write.csv(data.frame(event = top_hubs,
                         MM = mm_values[top_hubs]),
              paste0("hub_events_", module, "_module.csv"))
  }
}


write.csv(MEs, "MEs_export.csv", row.names = TRUE)
getwd()
list.files()


