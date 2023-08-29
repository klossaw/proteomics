
# Differential abundance of proteins --------------------------------------

.libPaths(new = c("~/sbin/R/R-4.2.1", .libPaths()))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "msqrob2", "gridExtra",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "QFeatures", "plotly")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "proteomics"
dataset <- "tutorial"
species <- "test"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/dif_abund")
workdir %>% checkdir() %>% setwd()

# Data source: "https://raw.githubusercontent.com/statOmics/PDA/data/quantification/francisella/peptides.txt"
pep_file <- glue("~/projects/{project}/data/peptides_dif_abund.txt")
ecols <- grep("Intensity\\.", names(read.delim(pep_file)))
# read files as "QFeatures"
pe <- readQFeatures(table = pep_file, fnames = 1, ecol = ecols, name = "peptideRaw", sep="\t")
names(pe)
dim(rowData(pe[[1]]))
dim(assay(pe[[1]]))
assay(pe[[1]])[1:4, 1:4]

rowData(pe[[1]])$nNonZero <- rowSums(assay(pe[[1]]) > 0)
colData(pe)$genotype <- pe[[1]] %>% colnames() %>% 
  substr(12, 13) %>% as.factor() %>% relevel("WT")
pe %>% colData()

pe[["protein"]] %>% assay() %>% head()

rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- pe %>% zeroIsNA(., "peptideRaw") %>% logTransform(., base = 2, i = "peptideRaw", name = "peptideLog")

pe <- filterFeatures(pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)) %>% 
  filterFeatures(., ~ Reverse != "+") %>% filterFeatures(., ~ Contaminant != "+") %>% 
  filterFeatures(., ~ nNonZero >= 2)

nrow(pe[[2]])
pe <- normalize(pe, i = "peptideLog", name = "peptideNorm", method = "center.median")

pe <- aggregateFeatures(pe, i = "peptideNorm", fcol = "Proteins", na.rm = TRUE, name = "protein")
gt_peptide_p <- pe[["peptideNorm"]] %>% assay %>% as.data.frame() %>%
  gather(sample, intensity) %>% mutate(genotype = colData(pe)[sample, "genotype"]) %>%
  ggplot(aes(x = intensity, group = sample, color = genotype)) + geom_density() + ggtitle("Peptide-level")
gt_protein_p <- pe[["protein"]] %>% assay %>% as.data.frame() %>%
  gather(sample, intensity) %>% mutate(genotype = colData(pe)[sample,"genotype"]) %>%
  ggplot(aes(x = intensity,group = sample,color = genotype)) + geom_density() + ggtitle("Protein-level")
ggsave(filename = "genotype_pept_prot_density.pdf", (gt_peptide_p + gt_protein_p), width = 8, height = 4)

prot <- "WP_003023392"
prot_df <- data.frame(
  intensity = assay(pe[["protein"]][prot, ]) %>% c(), 
  genotype = colData(pe)[, 1]) 

prot_pt <- prot_df %>% ggplot(aes(x=genotype, y=intensity)) + geom_point() + ggtitle(glue("Protein {prot}"))
t_test <- t.test(intensity ~ genotype, data = prot_df, var.equal=TRUE)

ttestMx <- function(y, group) {
  test <- try(t.test(y[group], y[!group], var.equal = TRUE), silent=TRUE)
  if(is(test,"try-error")) {
    return(c(log2FC = NA, se = NA, tstat = NA, p = NA))
  } else {
    return(c(log2FC = (test$estimate %*% c(1, -1)), se = test$stderr, tstat = test$statistic, pval = test$p.value))
  }
}

res <- apply(assay(pe[["protein"]]), 1, ttestMx,
  group = colData(pe)$genotype=="D8") %>% t() 
colnames(res) <- c("logFC", "se", "tstat", "pval")
res <- res %>% as.data.frame %>% na.exclude %>% arrange(pval)
res$adjPval <- p.adjust(res$pval, "fdr")
alpha <- 0.05
res$adjAlphaForm <- paste0(1:nrow(res)," x ",alpha,"/",nrow(res))
res$adjAlpha <- alpha * (1:nrow(res))/nrow(res) 
res$"pval < adjAlpha" <- res$pval < res$adjAlpha 
res$"adjPval < alpha" <- res$adjPval < alpha 

volcanoT <- res %>% 
  ggplot(aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) + geom_point(cex = 2.5) + 
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) + theme_minimal() 

problemPlots <- list() 
problemPlots[[1]] <- res %>% 
  ggplot(aes(x = logFC, y = se, color = adjPval < 0.05)) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_minimal() 

for (i in 2:3) {
  problemPlots[[i]] <- colData(pe) %>% 
    as.data.frame() %>% 
    mutate(intensity = pe[["protein"]][rownames(res)[i],] %>% assay() %>% c()) %>% 
    ggplot(aes(x=genotype,y=intensity)) + geom_point() + ylim(-3,0) + ggtitle(rownames(res)[i])
}


