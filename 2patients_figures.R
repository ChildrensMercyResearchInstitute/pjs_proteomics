library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(proDA)
library(EnhancedVolcano)
library(ggvenn)

### Proteomic data -----

proteo <- read.csv('inputs/2patients_proteo.csv')

#filter out contaminants
proteo <- proteo[1:4997, ]

proteo_abundance <- proteo[c(3, 4, 38:60)]
#incorporating specimen ID final value from SampleDescriptions excel sheet for colnames
colnames(proteo_abundance) <- c('Accession', 'GeneSymbol', 
                         '0441_Mucosa_02_J', '0441_Mucosa_03_I', '0441_Mucosa_05_A', 
                         '0441_Mucosa_07_A', '0441_Mucosa_09_A', '0441_Mucosa_10_I',
                         '0441_Mucosa_12_I', 
                         '0444_Mucosa_02_J', '0444_Mucosa_03_I', '0444_Mucosa_06_J', 
                         '0444_Mucosa_08_S', '0444_Mucosa_10_I', 
                         '0441_Polyp_01_J', '0441_Polyp_04_A', '0441_Polyp_06_A', 
                         '0441_Polyp_08_A', '0441_Polyp_11_I', 
                         '0444_Polyp_01_J', '0444_Polyp_04_I', '0444_Polyp_05_J', 
                         '0444_Polyp_07_S', '0444_Polyp_09_S', '0444_Polyp_11_I')

proteo_mat <- proteo_abundance[3:25]
rownames(proteo_mat) <- proteo_abundance$Accession

log_proteo <- log2(proteo_mat)

norm_log_proteo <- median_normalization(as.matrix(log_proteo))

## data frame of sample characteristics

proteo_ann_data <- as.data.frame(colnames(proteo_mat))
colnames(proteo_ann_data) <- 'sample'
rownames(proteo_ann_data) <- proteo_ann_data$sample

proteo_ann_data$patient <- c(rep('0441', 7), rep('0444', 5), rep('0441', 5), rep('0444', 6))
proteo_ann_data$tissue <- c(rep('Mucosa', 12), rep('Polyp', 11))
proteo_ann_data$location <- c('Jejunum', 'Ileum', 'Antrum', 'Antrum', 'Antrum', 
                       'Ileum', 'Ileum', 'Jejunum', 'Ileum', 'Jejunum', 
                       'Stomach', 'Ileum', 'Jejunum', 'Antrum', 'Antrum', 
                       'Antrum', 'Ileum', 'Jejunum', 'Ileum', 'Jejunum', 
                       'Stomach', 'Stomach', 'Ileum')
proteo_ann_data$location2 <- case_when(
  proteo_ann_data$location %in% c('Jejunum', 'Ileum') ~ 'SmallIntestine',
  proteo_ann_data$location %in% c('Stomach', 'Antrum') ~ 'Stomach',
  TRUE ~ 'Other'
)
proteo_ann_data$polyp_size_mm <- c(rep(0, 12), 10, 3, 3, 3, 15, 3, 3, 6, 3, 5, 40)

## sample similarity

#correlation matrix
#remove zeros for correlation calculation
norm_log_proteo_nozero <- norm_log_proteo
norm_log_proteo_nozero[is.na(norm_log_proteo_nozero)] <- 0

proteo_cor_mat <- cor(norm_log_proteo_nozero)

pdf('2patients_figures/proteo_cor_mat.pdf', height = 6, width = 6)
pheatmap(proteo_cor_mat, scale = "none", 
         annotation_col = proteo_ann_data[c('patient', 'tissue', 'location2')], 
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

#pca
t_proteo <- as.data.frame(t(norm_log_proteo_nozero))

proteo_pca <- prcomp(t_proteo, center = TRUE)

pdf('2patients_figures/proteo_pca.pdf', height = 4, width = 6)
fviz_pca_ind(proteo_pca, geom = "point", mean.point = FALSE, pointsize = 0, 
             habillage = proteo_ann_data$location, repel = TRUE) + 
  geom_point(aes(shape = proteo_ann_data$tissue, 
                 color = proteo_ann_data$location), size = 3) + 
  scale_shape_manual(name = 'tissue', values = c('Mucosa' = 16, 'Polyp' = 17))
dev.off()


### Phosphoproteomic data -----

phospho <- read.csv('inputs/2patients_phospho.csv')

#filter out contaminants
phospho <- phospho[1:2607, ]

phospho_abundance <- phospho[c(3, 5, 6, 95, 45:67)]
#incorporating specimen ID final value from SampleDescriptions excel sheet for colnames
colnames(phospho_abundance) <- c('Accession', 'GeneSymbol', 'Positions', 'Modifications_all', 
                         '0444_Mucosa_06_J', '0444_Mucosa_10_I', '0444_Mucosa_08_S', 
                         '0444_Mucosa_03_I', '0444_Mucosa_02_J', 
                         '0441_Mucosa_12_I', '0441_Mucosa_10_I', '0441_Mucosa_09_A', 
                         '0441_Mucosa_07_A', '0441_Mucosa_05_A', '0441_Mucosa_02_J', 
                         '0441_Mucosa_03_I', 
                         '0444_Polyp_11_I', '0444_Polyp_05_J', '0444_Polyp_04_I', 
                         '0444_Polyp_09_S', '0444_Polyp_07_S', '0444_Polyp_01_J', 
                         '0441_Polyp_11_I', '0441_Polyp_08_A', '0441_Polyp_06_A', 
                         '0441_Polyp_01_J', '0441_Polyp_04_A')

phospho_mat <- phospho_abundance[5:27]
rownames(phospho_mat) <- make.unique(phospho_abundance$Accession)

log_phospho <- log2(phospho_mat)

norm_log_phospho <- median_normalization(as.matrix(log_phospho))

## data frame of sample characteristics

phospho_ann_data <- as.data.frame(colnames(phospho_mat))
colnames(phospho_ann_data) <- 'sample'
rownames(phospho_ann_data) <- phospho_ann_data$sample

phospho_ann_data$patient <- c(rep('0444', 5), rep('0441', 7), rep('0444', 6), rep('0441', 5))
phospho_ann_data$tissue <- c(rep('Mucosa', 12), rep('Polyp', 11))
phospho_ann_data$location <- c('Jejunum', 'Ileum', 'Stomach', 'Ileum', 'Jejunum', 
                       'Ileum', 'Ileum', 'Antrum', 'Antrum', 'Antrum', 
                       'Jejunum', 'Ileum', 'Ileum', 'Jejunum', 'Ileum', 
                       'Stomach', 'Stomach', 'Jejunum', 'Ileum', 'Antrum', 
                       'Antrum', 'Jejunum', 'Antrum')
phospho_ann_data$location2 <- case_when(
  phospho_ann_data$location %in% c('Jejunum', 'Ileum') ~ 'SmallIntestine',
  phospho_ann_data$location %in% c('Stomach', 'Antrum') ~ 'Stomach',
  TRUE ~ 'Other'
)
phospho_ann_data$polyp_size_mm <- c(rep(0, 12), 40, 6, 3, 5, 3, 3, 15, 3, 3, 10, 3)

## sample similarity

#correlation matrix
#remove zeros for correlation calculation
norm_log_phospho_nozero <- norm_log_phospho
norm_log_phospho_nozero[is.na(norm_log_phospho_nozero)] <- 0

phospho_cor_mat <- cor(norm_log_phospho_nozero)

pdf('2patients_figures/phospho_cor_mat.pdf', height = 6, width = 6)
pheatmap(phospho_cor_mat, scale = "none", 
         annotation_col = phospho_ann_data[c('patient', 'tissue', 'location2')], 
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

#pca
t_phospho <- as.data.frame(t(norm_log_phospho_nozero))

phospho_pca <- prcomp(t_phospho, center = TRUE)

pdf('2patients_figures/phospho_pca.pdf', height = 4, width = 6)
fviz_pca_ind(phospho_pca, geom = "point", mean.point = FALSE, pointsize = 0, 
             habillage = phospho_ann_data$location, repel = TRUE) + 
  geom_point(aes(shape = phospho_ann_data$tissue, 
                 color = phospho_ann_data$location), size = 3) + 
  scale_shape_manual(name = 'tissue', values = c('Mucosa' = 16, 'Polyp' = 17))
dev.off()
