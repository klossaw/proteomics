
# Label-free proteomics data analysis -------------------------------------
.libPaths(new = c("~/sbin/R/R-4.2.1", .libPaths()))
pkgs <- c("fs", "stringr", "ggthemes", "msqrob2",
      "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "QFeatures", "plotly")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "proteomics"
dataset <- "tutorial"
species <- "test"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}")
workdir %>% fs::dir_create() %>% setwd()

# peptidesFile <- "https://raw.githubusercontent.com/statOmics/PDA/data/quantification/fullCptacDatasSetNotForTutorial/peptides.txt"
pep_file <- glue("~/projects/{project}/data/peptides.txt")
# grep the ion intensity
ecols <- grep("Intensity\\.", names(read.delim(pep_file)))
# read files as "QFeatures"
pe <- readQFeatures(table = pep_file, fnames = 1, ecol = ecols, name = "peptideRaw", sep="\t") # comma-separated values by default 
names(pe)
# "peptideRaw"
class(pe[[1]])
# "SummarizedExperiment"
rowData(pe[["peptideRaw"]]) %>% dim()
# Sequence information, 11466   143
colData(pe[["peptideRaw"]])
# DataFrame with 45 rows and 0 columns
# No information is stored yet on the design
pe[["peptideRaw"]] %>% colnames() %>% head()
# See the columns' names of all the samples

# Update the colData with information on the design
# "colData" returns a DFrame object, which could not be edited by tidyverse
# colData(pe) <- colData(pe) %>% mutate(lab = as.factor(rep(rep(paste0("lab", 1:3), each=3), 5)),
#                                       condition = pe[["peptideRaw"]] %>% colnames() %>% substr(12,12) %>% as.factor(),
#                                       spikeConcentration = rep(c(A = 0.25, B = 0.74, C = 2.22, D = 6.67, E = 20), each = 9))
# Error in UseMethod("mutate") :
#   no applicable method for 'mutate' applied to an object of class
# "c('DFrame', 'DataFrame', 'SimpleList', 'RectangularData', 'List', 'DataFrame_OR_NULL', 'Vector', 
# 'list_OR_List', 'ListorHits', 'Annotated', 'vector_OR_Vector')"
# 
colData(pe)$lab <- rep(rep(paste0("lab", 1:3), each=3), 5) %>% as.factor()
colData(pe)$condition <- pe[["peptideRaw"]] %>% colnames() %>% substr(12,12) %>% as.factor()
colData(pe)$spikeConcentration <- rep(c(A = 0.25, B = 0.74, C = 2.22, D = 6.67, E = 20), each = 9)
colData(pe) %>% typeof()

# Preprocessing ----
# Log-transformation
subset <- pe["AALEELVK", colData(pe)$lab=="lab1"]
plot_df <- data.frame(concentration = colData(subset)$spikeConcentration,
                      y = assay(subset[["peptideRaw"]]) %>% c())
plot_before_log <- plot_df %>% ggplot(aes(concentration, y)) + geom_point() +
  labs(x = "concentration (fmol/l)", title = "peptide AALEELVK in lab1")
plot_after_log <- plot_df %>% ggplot(aes(concentration, y)) + geom_point() +
  scale_x_continuous(trans = 'log2') + scale_y_continuous(trans = 'log2') +
  labs(x = "concentration (fmol/l)", title = "peptide AALEELVK in lab1 with axes on log2 scale")
# plot_after_ln <- plot_df %>% ggplot(aes(concentration, y)) + geom_point() +
#   scale_x_continuous(trans = 'log') + scale_y_continuous(trans = 'log') +
#   labs(x = "concentration (fmol/l)", title = "peptide AALEELVK in lab1 with axes on ln scale")
# plot_after_lg <- plot_df %>% ggplot(aes(concentration, y)) + geom_point() +
#   scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10') +
#   labs(x = "concentration (fmol/l)", title = "peptide AALEELVK in lab1 with axes on lg scale")
# plot_comb <- (plot_before_log + plot_after_log) / (plot_after_ln + plot_after_lg)
plot_comb <- (plot_before_log + plot_after_log)
ggsave("why_log.pdf", plot_comb, width = 10, height = 4)
# Note: it is the axis scale changed, x and y for each point are all the same
# when the base is 10 or e, it is similar when the base is 2

# How many non zero intensities we have for each peptide and this can be useful for filtering
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
# Convert 0 to NA
# Peptides with zero intensities are missing peptides and should be represent with a NA value 
pe <- zeroIsNA(pe, "peptideRaw") 
# Log transform data with base 2
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")

# Filtering
# Handling overlapping protein groups
pe[[1]] %>% rowData() %>% colnames() %>% grep(pattern = "Proteins", value = TRUE)
pe <- filterFeatures(pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins))
pe
# An instance of class QFeatures containing 2 assays:
# [1] peptideRaw: SummarizedExperiment with 10740 rows and 45 columns
# [2] peptideLog: SummarizedExperiment with 10740 rows and 45 columns
# Remove reverse sequences (decoys) and contaminants
rowData(pe[[1]])$Reverse %>% table()
# .
# 
# 10740
rowData(pe[[1]])$Potential.contaminant %>% table()
# .
# +
#   10678    62
pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")
pe
# An instance of class QFeatures containing 2 assays:
# [1] peptideRaw: SummarizedExperiment with 10678 rows and 45 columns
# [2] peptideLog: SummarizedExperiment with 10678 rows and 45 columns

# Drop peptides that were only identified in one sample
summary(rowData(pe[[1]])$nNonZero)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00   14.00   25.00   24.35   34.00   45.00
table(rowData(pe[[1]])$nNonZero > 3)
# FALSE  TRUE
# 516 10162
pe <- filterFeatures(pe, ~ nNonZero >= 2)
nrow(pe[["peptideLog"]])
# 10478
assays(pe[[1]]) %>% is.na() %>% table()

# Normalization
# Before normalization
density_cond_d <- pe[["peptideLog"]][, colData(pe)$condition=="D"] %>%
  assay() %>% as.data.frame() %>% 
  pivot_longer(names_to = "sample", values_to = "intensity", cols = 1:ncol(.)) %>%
  mutate(lab = colData(pe)[sample,"lab"]) %>%
  ggplot(aes(x = intensity, group = sample, color = lab)) +
  geom_density() + ggtitle("condition D -- before normalization")
density_lab2 <- pe[["peptideLog"]][, colData(pe)$lab == "lab2"] %>%
  assay %>% as.data.frame() %>%
  pivot_longer(names_to = "sample", values_to = "intensity", cols = 1:ncol(.)) %>%
  mutate(condition = colData(pe)[sample,"condition"]) %>%
  ggplot(aes(x = intensity, group = sample, color = condition)) +
  geom_density() + ggtitle("lab2 -- before normalization")
ggsave("density_before_norm.pdf", (density_cond_d + density_lab2), width = 9, height = 4)
# Normalization of the data by median centering
pe <- normalize(pe, i = "peptideLog", name = "peptideNorm", method = "center.median")
density_cond_d_norm <- pe[["peptideNorm"]][, colData(pe)$condition == "D"] %>%
  assay() %>% as.data.frame() %>%
  pivot_longer(names_to = "sample", values_to = "intensity", cols = 1:ncol(.)) %>%
  mutate(lab = colData(pe)[sample, "lab"]) %>%
  ggplot(aes(x = intensity,group = sample,color = lab)) +
  geom_density() + ggtitle("condition D -- after normalization")
density_lab2_norm <- pe[["peptideNorm"]][, colData(pe)$lab == "lab2"] %>%
  assay() %>% as.data.frame() %>%
  pivot_longer(names_to = "sample", values_to = "intensity", cols = 1:ncol(.)) %>%  
  mutate(condition = colData(pe)[sample, "condition"]) %>%
  ggplot(aes(x = intensity, group = sample, color = condition)) +
  geom_density() + ggtitle("lab2 -- after normalization")
ggsave("density_after_norm.pdf", (density_cond_d_norm + density_lab2_norm), width = 9, height = 4)

# Explore normalized data
pdf("peptide_norm.pdf", width = 5, height = 5)
  limma::plotDensities(assay(pe[["peptideNorm"]]))
  boxplot(assay(pe[["peptideNorm"]]), col = palette()[-1],
          main = "Peptide distribtutions after normalisation", ylab = "intensity")
  limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))
dev.off()

# Summarization, codes below are copied from package 'msqrob2' tutorial
summary_plot <- pe[["peptideNorm"]][
  rowData(pe[["peptideNorm"]])$Proteins == "P12081ups|SYHC_HUMAN_UPS",
  colData(pe)$lab == "lab2" & colData(pe)$condition %in% c("A", "E")] %>%
  assay() %>% as.data.frame() %>% rownames_to_column(var = "peptide") %>%
  gather(sample, intensity, -peptide) %>%
  mutate(condition = colData(pe)[sample,"condition"]) %>%
  ggplot(aes(x = peptide, y = intensity, color = sample, group = sample, label = condition)) + 
  geom_line() + geom_text() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  xlab("Peptide") + ylab("Intensity (log2)")
ggsave("summary_plot.pdf", summary_plot, width = 9, height = 5)
# Aggregate peptides to proteins
pe <- aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")

# MDS plot of proteins
pdf("protein_sum.pdf", width = 5, height = 5)
  limma::plotDensities(assay(pe[["protein"]]))
  boxplot(assay(pe[["protein"]]), col = palette()[-1],
        main = "Protein distribtutions after normalisation", ylab = "intensity")
  limma::plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))
dev.off()

# Estimation, rlm as default
pe <- msqrob(object = pe, i = "protein", formula = ~condition)
# Inference
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
# # Inference
# getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
# There were 50 or more warnings (use warnings() to see the first 50)
# (Intercept)  conditionB  conditionC  conditionD  conditionE
# -3.883933    1.531672    2.801236    4.588426    6.071540

L <- msqrob2::makeContrast("conditionB=0", parameterNames = c("conditionB"))
pe <- msqrob2::hypothesisTest(object = pe, i = "protein", contrast = L)


volcano <- ggplot(rowData(pe[["protein"]])$conditionB,
                  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) + geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() + ggtitle("Default workflow")
ggsave("volcano_protein_cond_b.pdf", volcano, width = 5, height = 5)

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
  pdf("detailed_plot_top5.pdf", width = 10, height = 6)
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
# It seems that the default processing not suitable for this complicated design
