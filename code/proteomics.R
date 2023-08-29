
# Proteomics data analysis, based on "msqrob2" ------------------------------------

.libPaths(new = c("~/sbin/R/R-4.2.1", .libPaths()))
# BiocManager::install(c("QFeatures", "msqrob2"))

pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
      "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "limma",
      "QFeatures", "msqrob2", "plotly", "gridExtra", "msdata")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "proteomics"
dataset <- "test"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}")
workdir %>% checkdir() %>% setwd()

peptidesFile <- msdata::quant(pattern = "cptac_a_b_peptides", full.names = TRUE)

ecols <- grep("Intensity\\.", names(read.delim(peptidesFile)))

pe <- readQFeatures(table = peptidesFile, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t")

cond <- which(strsplit(colnames(pe)[[1]][1], split = "")[[1]] == "A") # find where condition is stored
colData(pe)$condition <- substr(colnames(pe), cond, cond) %>% unlist() %>% as.factor()

# Calculate how many non zero intensities we have per peptide
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)

# Peptides with zero intensities are missing peptides and should be represent with a NA value rather than 0
pe <- QFeatures::zeroIsNA(pe, "peptideRaw") # convert 0 to NA

# Inspect the missingness in our data with the plotNA() function provided with MSnbase
pdf("na_plot.pdf", width = 5, height = 5)
  MSnbase::plotNA(assay(pe[["peptideRaw"]])) + xlab("Peptide index (ordered by data completeness)")
dev.off() # Why would it occur twice?


# Pre-processing

# Log-transform
pe <- QFeatures::logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
pdf("peptide_log.pdf", width = 5, height = 5)
  limma::plotDensities(assay(pe[["peptideLog"]]))
dev.off()

# Filtering, handling the overlapped proteins
# A peptide can map to multiple proteins, 
# as long as there is none of these proteins present in a smaller subgroup
Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% 
  msqrob2::smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]
dim(pe[[1]])
# 10740     6

# Remove reverse sequences (decoys) and contaminants
pe <- filterFeatures(pe, ~ Reverse != "+")
# 'Reverse' found in 2 out of 2 assay(s)
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")
# 'Potential.contaminant' found in 2 out of 2 assay(s)

# Remove peptides of proteins that were only identified with modified peptides
# Large protein groups file needed for this

# Drop peptides that were only identified in one sample
nrow(pe[["peptideLog"]])
# 10678
pe <- filterFeatures(pe, ~ nNonZero >= 2)

nrow(pe[["peptideLog"]])
# 7011

# Normalize the data by median centering
pe <- QFeatures::normalize(pe, i = "peptideLog", name = "peptideNorm", method = "center.median")
# ‘"center.mean"’ and ‘"center.median"’ center the respective
# sample (column) intensities by subtracting the respective
# column means or medians. ‘"div.mean"’ and ‘"div.median"’
# divide by the column means or medians. These are equivalent
# to ‘sweep’ing the column means (medians) along ‘MARGIN = 2’
# with ‘FUN = "-"’ (for ‘"center.*"’) or ‘FUN = "/"’ (for ‘"div.*"’).

# Explore normalized data
pdf("peptide_norm.pdf", width = 5, height = 5)
  limma::plotDensities(assay(pe[["peptideNorm"]]))
  boxplot(assay(pe[["peptideNorm"]]), col = palette()[-1],
          main = "Peptide distribtutions after normalisation", ylab = "intensity")
  limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))
dev.off()

# Aggregate peptides to proteins
pe <- QFeatures::aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", 
                                   na.rm = TRUE, name = "protein")

# MDS plot of proteins
pdf("protein_mds.pdf", width = 5, height = 5)
  limma::plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))
dev.off()

# Estimation, rlm as default
pe <- msqrob(object = pe, i = "protein", formula = ~condition)
# Inference
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
# (Intercept)  conditionB
# -2.672396    1.513682
L <- msqrob2::makeContrast("conditionB=0", parameterNames = c("conditionB"))
pe <- msqrob2::hypothesisTest(object = pe, i = "protein", contrast = L)


volcano <- ggplot(rowData(pe[["protein"]])$conditionB,
  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) + geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() + ggtitle("Default workflow")
ggsave("volcano_protein.pdf", volcano, width = 5, height = 5)

# Heatmap of significant features
sigNames <- rowData(pe[["protein"]])$conditionB %>%
  rownames_to_column("protein") %>%
  filter(adjPval < 0.05) %>% pull(protein)
pdf("sig_pro_heatmap.pdf", width = 7, height = 8)
  ComplexHeatmap::Heatmap(assay(pe[["protein"]])[sigNames, ])
  heatmap(assay(pe[["protein"]])[sigNames, ])
dev.off()

# Detail plot of top 5 proteins
topN <- 5
if (length(sigNames) > topN) {
  pdf("detailed_plot_top5.pdf", width = 8, height = 5)
    for(protName in sigNames[1:topN]) {
      pePlot <- pe[protName, , c("peptideNorm", "protein")]
      pePlotDf <- data.frame(longFormat(pePlot))
      pePlotDf$assay <- factor(pePlotDf$assay,
                               levels = c("peptideNorm", "protein")
      )
      pePlotDf$condition <- as.factor(colData(pePlot)[pePlotDf$colname, "condition"])
    
      # plotting
      p1 <- ggplot(data = pePlotDf, aes(x = colname, y = value, group = rowname)) +
        geom_line() + geom_point() + theme_minimal() + facet_grid(~assay) + ggtitle(protName)
      print(p1)
      
      # plotting 2
      p2 <- ggplot(pePlotDf, aes(x = colname, y = value, fill = condition)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = position_jitter(width = .1), aes(shape = rowname)) +
        scale_shape_manual(values = 1:nrow(pePlotDf)) +
        labs(title = protName, x = "sample", y = "peptide intensity (log2)") +
        theme_minimal() + facet_grid(~assay)
      print(p2)
    }
  dev.off()
}

# Comparison with other workflows
# Median summarisation
pe <- aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                        name = "proteinMedian", fun = matrixStats::colMedians)
pe <- msqrob(object = pe, i = "proteinMedian", formula = ~condition)
pe <- hypothesisTest(object = pe, i = "proteinMedian", contrast = L)

pdf("protein_median_mds_plot.pdf", width = 5, height = 5)
  limma::plotMDS(assay(pe[["proteinMedian"]]), col = as.numeric(colData(pe)$condition))
dev.off()

volcanoMed <- rowData(pe[["proteinMedian"]])[[colnames(L)]] %>%
  ggplot(aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5) + theme_minimal() + ggtitle("median summarisation") +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  geom_vline(xintercept = log2(0.74 / .25), col = "red")
ggsave("volcano_median.pdf", volcanoMed)

# Robust summarisation followed by robust ridge regression
try(pe <- msqrob(object = pe, i = "protein", formula = ~condition,
  modelColumnName = "ridge", ridge = TRUE)) # note: intentional error
## The mean model must have more than two parameters for ridge regression.
# if you really want to adopt ridge regression when your factor has only two levels
# rerun the function with a formula where you drop the intercept. e.g. ~-1+condition
pe <- msqrob(object = pe, i = "protein", formula = ~ -1 + condition,
  modelColumnName = "ridge", ridge = TRUE)
Lridge <- makeContrast("ridgeconditionB - ridgeconditionA = 0",
  c("ridgeconditionB", "ridgeconditionA"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = Lridge, modelColumn = "ridge")
volcanoRidge <- rowData(pe[["protein"]])[[colnames(Lridge)]] %>%
  ggplot(aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5) + theme_minimal() + ggtitle(paste("robust ridge")) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  geom_vline(xintercept = log2(0.74 / 0.25), col = "red")
ggsave("volcano_ridge.pdf", volcanoRidge, width = 5, height = 5)
pdf("volcano_grid.pdf", width = 4, height = 10)
  grid.arrange(grobs = list(volcano, volcanoMed, volcanoRidge), ncol = 1)
dev.off()

rowData(pe[["protein"]])$ups <- grepl("UPS", rownames(pe[["protein"]]))
rowData(pe[["proteinMedian"]])$ups <- grepl("UPS", rownames(pe[["proteinMedian"]]))
logFC <- data.frame(default = rowData(pe[["protein"]])[[colnames(L)]][, 1],
  median = rowData(pe[["proteinMedian"]])[[colnames(L)]][, 1],
  ridge = rowData(pe[["protein"]])[[colnames(Lridge)]][, 1],
  ups = rowData(pe[["protein"]])$ups)

logFC <- logFC %>% pivot_longer(names_to = "method", values_to = "log2FC", cols = 1:3)
logFC$ups <- as_factor(logFC$ups)
boxplot_logfc <- logFC %>% ggplot(aes(x = method, y = log2FC, fill = ups)) + 
  geom_boxplot() + geom_hline(yintercept = log2(0.74 / .25), color = "red")
ggsave("boxplot_logfc.pdf", boxplot_logfc, width = 5, height = 5)

# Sensitivity and accuracy estimation: TPR (True Positive Rate) and FDP (False Discovery Proportion)
tprFdp <- function(pval, tp, adjPval) {
  ord <- order(pval)
  return(data.frame(pval = pval[ord], adjPval = adjPval[ord],
    tpr = cumsum(tp[ord]) / sum(tp), fdp = cumsum(!tp[ord]) / 1:length(tp)))
}
tprFdpDefault <- tprFdp(rowData(pe[["protein"]])[[colnames(L)]]$pval,
  rowData(pe[["protein"]])$ups, rowData(pe[["protein"]])[[colnames(L)]]$adjPval)
tprFdpMedian <- tprFdp(rowData(pe[["proteinMedian"]])[[colnames(L)]]$pval,
  rowData(pe[["proteinMedian"]])$ups, rowData(pe[["proteinMedian"]])[[colnames(L)]]$adjPval)

tprFdpRidge <- tprFdp(rowData(pe[["protein"]])[[colnames(Lridge)]]$pval,
  rowData(pe[["protein"]])$ups, rowData(pe[["protein"]])[[colnames(Lridge)]]$adjPval)

hlp <- rbind(cbind(tprFdpDefault, method = "default"), cbind(tprFdpMedian, method = "median"),
  cbind(tprFdpRidge, method = "ridge"))
tprFdpPlot <- hlp %>% ggplot(aes(x = fdp, y = tpr, color = method)) + geom_path()


