
# Full CPTAC study ------------

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

# Data source: https://raw.githubusercontent.com/statOmics/PDA/data/quantification/fullCptacDatasSetNotForTutorial/peptides.txt
pep_file <- glue("~/projects/{project}/data/peptides_full_cptac.txt")
# grep the ion intensity
ecols <- grep("Intensity\\.", names(read.delim(pep_file)))
# read files as "QFeatures"
pe <- readQFeatures(table = pep_file, fnames = 1, ecol = ecols, name = "peptideRaw", sep="\t")
names(pe)
# "peptideRaw"

pe %>% colnames
# CharacterList of length 1
# [["peptideRaw"]] Intensity.6A_1 Intensity.6A_2 ... Intensity.6E_9

colData(pe)$lab <- rep(rep(paste0("lab", 1:3), each=3), 5) %>% as.factor
colData(pe)$condition <- pe[["peptideRaw"]] %>% colnames %>% substr(12,12) %>% as.factor
colData(pe)$spikeConcentration <- rep(c(A = 0.25, B = 0.74, C = 2.22, D = 6.67, E = 20), each = 9)

# log-transformed
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") %>% logTransform(., base = 2, i = "peptideRaw", name = "peptideLog")

# Filtering
pe <- filterFeatures(pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)) %>% 
  filterFeatures(., ~Reverse != "+") %>% filterFeatures(., ~ Potential.contaminant != "+") %>% 
  filterFeatures(., ~ nNonZero >= 2)
pe
# An instance of class QFeatures containing 2 assays:
# [1] peptideRaw: SummarizedExperiment with 10478 rows and 45 columns
# [2] peptideLog: SummarizedExperiment with 10478 rows and 45 columns

pe <- normalize(pe, i = "peptideLog", name = "peptideNorm", method = "center.median")

prot <- "P01031ups|CO5_HUMAN_UPS"
data <- pe[["peptideNorm"]][rowData(pe[["peptideNorm"]])$Proteins == prot,
  colData(pe)$lab == "lab3"] %>% assay() %>% as.data.frame() %>% 
  rownames_to_column(var = "peptide") %>% gather(sample, intensity, -peptide) %>% 
  mutate(condition = colData(pe)[sample, "condition"]) %>% na.exclude()
sumPlot <- data %>%
  ggplot(aes(x = peptide, y = intensity, color = condition, group = sample, label = condition), 
         show.legend = FALSE) + geom_text(show.legend = FALSE) + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Peptide") + ylab("Intensity (log2)") + ggtitle(paste0("protein: ", prot))
sumPlot2 <- sumPlot + geom_line(linetype="dashed", alpha=.4)
ggsave("sum_plot.pdf", (sumPlot / sumPlot2), width = 4, height = 8)

# Median summarization
dataHlp <- pe[["peptideNorm"]][rowData(pe[["peptideNorm"]])$Proteins == prot,
  colData(pe)$lab=="lab3"] %>% assay()

sumMedian <- data.frame(intensity = dataHlp %>% colMedians(na.rm = TRUE),
  condition = colnames(dataHlp) %>% substr(12,12) %>% as.factor())

sumMedianPlot <- sumPlot + ggtitle("Median summarization") +
  geom_hline(data = sumMedian, mapping = aes(yintercept = intensity, color = condition))

# Mean summarization
sumMeanMod <- lm(intensity ~ -1 + sample,data)

sumMean <- data.frame(intensity = sumMeanMod$coef[grep("sample", names(sumMeanMod$coef))],
  condition = names(sumMeanMod$coef)[grep("sample", names(sumMeanMod$coef))] %>% 
    substr(18, 18) %>% as.factor())

sumMeanPlot <- sumPlot + ggtitle("Mean summarization") + 
  geom_hline(data = sumMean, mapping = aes(yintercept=intensity,color=condition))

pdf("sum_median_grid.pdf", width = 7, height = 4)
  grid.arrange(sumMedianPlot, sumMeanPlot, ncol=2)
dev.off()

# Model-based summarization
sumMeanPepMod <- lm(intensity ~ -1 + sample + peptide, data)

sumMeanPep <- data.frame(
  intensity = sumMeanPepMod$coef[grep("sample", names(sumMeanPepMod$coef))] + 
    mean(data$intensity) - mean(sumMeanPepMod$coef[grep("sample",names(sumMeanPepMod$coef))]),
  condition = names(sumMeanPepMod$coef)[grep("sample", names(sumMeanPepMod$coef))] %>% 
    substr(18,18) %>% as.factor())

fitLmPlot <-  sumPlot + ggtitle("fit: ~ sample + peptide") + 
  geom_line(data = data %>% mutate(fit = sumMeanPepMod$fitted.values), 
            mapping = aes(x = peptide, y=fit, color = condition, group = sample))
sumLmPlot <- sumPlot + ggtitle("Summarization: sample effect") + 
  geom_hline(data = sumMeanPep, mapping = aes(yintercept = intensity, color = condition))

pdf("sum_model_based.pdf", width = 10, height = 4)  
  grid.arrange(sumMedianPlot, sumMeanPlot, sumLmPlot, nrow=1)
dev.off()

# 
sumMeanPepRobMod <- MASS::rlm(intensity ~ -1 + sample + peptide, data)
resRobPlot <- data %>% mutate(res = sumMeanPepRobMod$residuals, w = sumMeanPepRobMod$w) %>%
  ggplot(aes(x = peptide, y = res, color = condition, label = condition,size=w), show.legend = FALSE) +
  geom_point(shape = 21, size = .2) + geom_point(shape = 21) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Peptide") + ylab("residual") + ylim(c(-1, 1)*max(abs(sumMeanPepRobMod$residuals)))
weightPlot <- qplot(seq(-5, 5, .01), MASS::psi.huber(seq(-5, 5, .01)),
  geom="path") + xlab("standardized residual") + ylab("weight")
ggsave("res_robust_sum.pdf", (resRobPlot + weightPlot), width = 8, height = 5)

# 
sumMeanPepRob <- data.frame(
  intensity = sumMeanPepRobMod$coef[grep("sample", names(sumMeanPepRobMod$coef))] + mean(data$intensity) - 
    mean(sumMeanPepRobMod$coef[grep("sample", names(sumMeanPepRobMod$coef))]),
  condition = names(sumMeanPepRobMod$coef)[grep("sample",names(sumMeanPepRobMod$coef))] %>% 
    substr(18, 18) %>% as.factor())

sumRlmPlot <- sumPlot + ggtitle("Robust") + 
  geom_hline(data = sumMeanPepRob, mapping = aes(yintercept = intensity, color = condition))

pdf("sum_lm_rlm_plot.pdf", width = 8, height = 4)  
  grid.arrange(sumLmPlot + ggtitle("OLS"), sumRlmPlot, nrow = 1)
dev.off()

# Robust regresion results in a better separation between the protein expression values 
# for the different samples according to their spike-in concentration
