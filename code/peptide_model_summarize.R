
# Peptide model for summarisation -----------------------------------------
.libPaths(new = c("~/sbin/R/R-4.2.1", .libPaths()))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "msqrob2", "gridExtra",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "QFeatures", "plotly")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "proteomics"
dataset <- "tutorial"
species <- "test"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/pep_model")
workdir %>% checkdir() %>% setwd()

# Data source: "https://raw.githubusercontent.com/statOmics/PDA/data/quantification/cptacAvsB_lab3/proteinGroups.txt"
protein_file <-  glue("~/projects/{project}/data/protein_groups.txt")
ecols <- grep("LFQ\\.intensity\\.", names(read.delim(protein_file)))
pe_lfq <- readQFeatures(table = protein_file, fnames = 1, ecol = ecols, name = "protein_raw", sep = "\t")
cond <- which(strsplit(colnames(pe_lfq)[[1]][1], split = "")[[1]] == "A") 
# find where condition is stored

colData(pe_lfq)$condition <- substr(colnames(pe_lfq), cond, cond) %>%
  unlist() %>% as.factor()

rowData(pe_lfq[["protein_raw"]])$nNonZero <- rowSums(assay(pe_lfq[["protein_raw"]]) > 0)

pe_lfq <- zeroIsNA(pe_lfq, "protein_raw") %>% logTransform(., base = 2, i = "protein_raw", name = "protein_log")
# Note: when a new assay is created, it could not be passed to the next
pe_lfq <- filterFeatures(pe_lfq, ~ Reverse != "+") %>% filterFeatures(., ~ Potential.contaminant != "+") %>%
  normalize(object = pe_lfq, i = "protein_log", name = "protein", method = "center.median")
pe_lfq <- msqrob(object = pe_lfq, i = "protein", formula = ~condition)

L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
pe_lfq <- hypothesisTest(object = pe_lfq, i = "protein", contrast = L, overwrite = TRUE)

volcanoLFQ <- ggplot(rowData(pe_lfq[["protein"]])$conditionB, 
                     aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5) + scale_color_manual(values = alpha(c("black", "red"), 0.5)) + theme_minimal() +
  ggtitle(paste0("maxLFQ: TP = ",
                 sum(rowData(pe_lfq[["protein"]])$conditionB$adjPval<0.05 & 
                       grepl(rownames(rowData(pe_lfq[["protein"]])$conditionB),pattern ="UPS"), na.rm=TRUE), 
                 " FP = ", 
                 sum(rowData(pe_lfq[["protein"]])$conditionB$adjPval<0.05 &! 
                       grepl(rownames(rowData(pe_lfq[["protein"]])$conditionB),pattern ="UPS"), na.rm=TRUE)))
ggsave("volcano_lfq.pdf", volcanoLFQ, width = 6, height = 6)


# Median and robust normalization -----------------------------------------
# Data source: "https://raw.githubusercontent.com/statOmics/SGA2020/data/quantification/cptacAvsB_lab3/peptides.txt"
peptides_file <- glue("~/projects/{project}/data/peptides_file.txt")
ecols <- grep("Intensity\\.", names(read.delim(peptides_file)))
pe <- readQFeatures(table = peptides_file, fnames = 1, ecol = ecols, name = "peptideRaw", sep="\t")
cond <- which(strsplit(colnames(pe)[[1]][1], split = "")[[1]] == "A") # find where condition is stored
colData(pe)$condition <- substr(colnames(pe), cond, cond) %>% unlist() %>% as.factor()

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
nrow(pe[["peptideRaw"]]) # 11466
pe <- zeroIsNA(pe, "peptideRaw") %>% logTransform(., base = 2, i = "peptideRaw", name = "peptideLog") 
pe <- filterFeatures(object = pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)) %>% 
  filterFeatures(., ~Reverse != "+") %>% filterFeatures(., ~ Potential.contaminant != "+") %>% 
  filterFeatures(., ~ nNonZero >=2)
nrow(pe[["peptideLog"]]) # 7011
pe <- normalize(pe, i = "peptideLog", name = "peptideNorm", method = "center.median")
# Different summarising method
pe <- aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "proteinMedian",
                        fun = matrixStats::colMedians)
pe <- aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "proteinRobust")

# Modeling and inference
pe <- msqrob(object = pe, i = "proteinMedian", formula = ~condition)
L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
pe <- hypothesisTest(object = pe, i = "proteinMedian", contrast = L)

pe <- msqrob(object = pe, i = "proteinRobust", formula = ~condition)
pe <- hypothesisTest(object = pe, i = "proteinRobust", contrast = L)

volcanoMedian <- ggplot(rowData(pe[["proteinMedian"]])$conditionB,
                        aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5) + scale_color_manual(values = alpha(c("black", "red"), 0.5)) + theme_minimal() +
  ggtitle(paste0("Median: TP = ",
                 sum(rowData(pe[["proteinMedian"]])$conditionB$adjPval < 0.05 & 
                    grepl(rownames(rowData(pe[["proteinMedian"]])$conditionB), pattern = "UPS"), na.rm = TRUE), 
                 " FP = ", 
                 sum(rowData(pe[["proteinMedian"]])$conditionB$adjPval < 0.05 &! 
                    grepl(rownames(rowData(pe[["proteinMedian"]])$conditionB), pattern = "UPS"), na.rm = TRUE)))

volcanoRobust<- ggplot(rowData(pe[["proteinRobust"]])$conditionB,
                       aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5) + scale_color_manual(values = alpha(c("black", "red"), 0.5)) + theme_minimal() +
  ggtitle(paste0("Median: TP = ",
                 sum(rowData(pe[["proteinRobust"]])$conditionB$adjPval < 0.05 & 
                       grepl(rownames(rowData(pe[["proteinRobust"]])$conditionB), pattern = "UPS"), na.rm = TRUE), 
                 " FP = ", 
                 sum(rowData(pe[["proteinRobust"]])$conditionB$adjPval < 0.05 &! 
                       grepl(rownames(rowData(pe[["proteinRobust"]])$conditionB), pattern = "UPS"), na.rm = TRUE)))

ylims <- c(0, ceiling(max(c(-log10(rowData(pe_lfq[["protein"]])$conditionB$pval),
                         -log10(rowData(pe[["proteinMedian"]])$conditionB$pval),
                         -log10(rowData(pe[["proteinRobust"]])$conditionB$pval)),
                       na.rm=TRUE)))

xlims <- max(abs(c(rowData(pe_lfq[["protein"]])$conditionB$logFC,
                   rowData(pe[["proteinMedian"]])$conditionB$logFC,
                   rowData(pe[["proteinRobust"]])$conditionB$logFC)),
             na.rm=TRUE) * c(-1, 1)

# Comparison summarization methods ----
pdf("grid_volcanos.pdf", width = 5, height = 10)
  grid.arrange(volcanoLFQ + xlim(xlims) + ylim(ylims), 
               volcanoMedian + xlim(xlims) + ylim(ylims), 
               volcanoRobust + xlim(xlims) + ylim(ylims),
               ncol = 1)
dev.off()

compBoxPlot <- rbind(rowData(pe_lfq[["protein"]])$conditionB %>% mutate(method = "maxLFQ") %>% rownames_to_column(var = "protein"),
                     rowData(pe[["proteinMedian"]])$conditionB %>% mutate(method = "median")%>% rownames_to_column(var = "protein"),
                     rowData(pe[["proteinRobust"]])$conditionB %>% mutate(method = "robust")%>% rownames_to_column(var = "protein")) %>%
  mutate(ups= grepl(protein, pattern="UPS")) %>%
  ggplot(aes(x = method, y = logFC, fill = ups)) + geom_boxplot() +
  geom_hline(yintercept = log2(0.74 / .25), color = "#00BFC4") +
  geom_hline(yintercept = 0, color = "#F8766D")
ggsave("comp_boxplot.pdf", compBoxPlot)



