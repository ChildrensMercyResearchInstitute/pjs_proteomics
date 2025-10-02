library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(proDA)
library(EnhancedVolcano)
library(Seurat)
library(ggvenn)
library(limma)
library(qdapRegex)

### Metadata for all samples

ann_data <- read.csv('inputs/sample_metadata.csv')
ann_data$patient <- as.character(ann_data$patient)
rownames(ann_data) <- ann_data$sample

### Proteomic data -----

proteo <- read.csv('inputs/all_patients_proteo.csv')
proteo <- proteo[1:4405, ]

#pull out raw abundance
proteo_abundance <- proteo[c(3, 4, 41:66)]
#purposely ensured that samples in same order in ann_data
colnames(proteo_abundance) <- c('Accession', 'GeneSymbol', ann_data$sample)

proteo_mat <- proteo_abundance[3:28]
rownames(proteo_mat) <- proteo_abundance$Accession

log_proteo <- log2(proteo_mat)

#applying 'conservative' median normalization
norm_log_proteo <- median_normalization(as.matrix(log_proteo))

#replace 0s with NAs
norm_log_proteo_nozero <- norm_log_proteo
norm_log_proteo_nozero[is.na(norm_log_proteo_nozero)] <- 0

## batch correction
proteo_design <- model.matrix(~tissue, ann_data)
norm_log_proteo_nozero_bc <- limma::removeBatchEffect(norm_log_proteo_nozero, 
                                               ann_data$batch, design = proteo_design)

## save batch-corrected counts
write.csv(norm_log_proteo_nozero_bc, 'proteo_norm_log_counts_bc.csv', row.names = TRUE, quote = FALSE)

## correlation
proteo_cor_bc <- cor(norm_log_proteo_nozero_bc)

pdf('all_patients_figures/proteo_cor_mat.pdf', height = 6, width = 6)
pheatmap(proteo_cor_bc, scale = "none", 
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_data[c('patient', 'tissue', 'location')])
dev.off()

## differential expression - polyps vs. mucosa

DEstats_proteo <- proteo[c(3, 4, 9, 10)]
colnames(DEstats_proteo) <- c('Accession', 'GeneSymbol', 'log2FC', 'pValue')

## filter down to statistically significant 
DE_proteo <- DEstats_proteo[abs(DEstats_proteo$log2FC) > 1 & DEstats_proteo$pValue < 0.05, ]

write.csv(DE_proteo, 'all_patients_figures/proteo_DE_sig.csv', row.names = FALSE, quote = FALSE)

pdf('all_patients_figures/proteo_DE_volcano.pdf', height = 8, width = 6)
EnhancedVolcano(DEstats_proteo, 
                lab = DEstats_proteo$GeneSymbol, 
                x = 'log2FC', 
                y = 'pValue', 
                #labSize = 0, 
                pCutoff = 0.05, 
                title = 'Polyp vs. Mucosa - proteomics')
dev.off()


### Phosphoproteomic data -----

phospho <- read.csv('inputs/all_patients_phospho.csv')
phospho <- phospho[1:1516, ]

# match sample order to ann_data/proteo data
phospho_abundance <- phospho[c(3, 5, 6, 103, 
                        53, 54, 52, 51, 50, 49, 47, 48, 55, 56, 57, 58, 
                        64, 63, 62, 61, 60, 59, 70, 65, 66, 67, 68, 69, 71, 72)]

colnames(phospho_abundance) <- c('Accession', 'GeneSymbol', 'Positions', 'Modifications', 
                         ann_data$sample)

phospho_mat <- phospho_abundance[5:30]
rownames(phospho_mat) <- make.unique(phospho_abundance$Accession)

log_phospho <- log2(phospho_mat)

#applying 'conservative' median normalization
norm_log_phospho <- median_normalization(as.matrix(log_phospho))

#replace 0s with NAs
norm_log_phospho_nozero <- norm_log_phospho
norm_log_phospho_nozero[is.na(norm_log_phospho_nozero)] <- 0

## batch correction
phospho_design <- model.matrix(~tissue, ann_data)
norm_log_phospho_nozero_bc <- limma::removeBatchEffect(norm_log_phospho_nozero, 
                                                      ann_data$batch, design = phospho_design)

## save batch-corrected counts
write.csv(norm_log_phospho_nozero_bc, 'phospho_norm_log_counts_bc.csv', row.names = TRUE, quote = FALSE)

## correlation
phospho_cor_bc <- cor(norm_log_phospho_nozero_bc)

pdf('all_patients_figures/phospho_cor_mat.pdf', height = 6, width = 6)
pheatmap(phospho_cor_bc, scale = "none", 
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_data[c('patient', 'tissue', 'location')])
dev.off()

## differential expression - polyps vs. mucosa

DEstats_phospho <- phospho[c(3, 5, 103, 14, 15)]
colnames(DEstats_phospho) <- c('Accession', 'GeneSymbol', 'Modifications', 'log2FC', 'pValue')

## filter down to statistically significant 
DE_phospho <- DEstats_phospho[abs(DEstats_phospho$log2FC) > 1 & DEstats_phospho$pValue < 0.05, ]

write.csv(DE_phospho, 'all_patients_figures/phospho_DE_sig.csv', row.names = FALSE, quote = FALSE)


### Preparing phosphoproteomic data for KSEA App -----

KSEA_stats <- phospho[c('Master.Protein.Accessions', 'Gene.Symbol', 
                         'Modifications.in.Master.Proteins', 
                         'Annotated.Sequence', 
                         'Abundance.Ratio...Polyp.....Mucosa.',                                
                         'Abundance.Ratio..log2....Polyp.....Mucosa.',                         
                         'Abundance.Ratio.P.Value...Polyp.....Mucosa.')]
colnames(KSEA_stats) <- c('Accession', 'GeneSymbol', 'Modifications', 'Sequence', 
                          'AR', 'log2AR', 'pval')

#remove any lines without phospho
KSEA_stats2 <- KSEA_stats[str_detect(KSEA_stats$Modifications, 'Phospho'), ]

#unlist any modifications that have ; between them - specific to "];" to capture distinct ones
KSEA_stats2$Mods2 <- strsplit(as.character(KSEA_stats2$Modifications), '];')
stats3 <- unnest(KSEA_stats2, Mods2)

#remove any lines without 'phospho'
stats3 <- stats3[str_detect(stats3$Mods2, 'Phospho'), ]

#fix any [] that got broken by just adding a ] to all Mods2
stats3$Mods2 <- trimws(stats3$Mods2)
stats3$Mods2 <- paste0(stats3$Mods2, ']')

#pull out accession (since some 'Accession' have two)
stats3$Accession2 <- vapply(strsplit(as.character(stats3$Mods2), ' '), '[', '', 1)
#fix it for those with 'Phospho' in Accession2
stats3$Accession3 <- case_when(
  str_detect(stats3$Accession2, 'Phospho') ~ stats3$Accession, 
  TRUE ~ stats3$Accession2
)

#pull out only content in [] in Mods column
stats3$Mods3 <- rm_between(stats3$Mods2, '[', ']', extract = TRUE)

#get rid of the () in that field
stats3$Mods3 <- gsub('\\([^\\)]+\\)', '', stats3$Mods3)
#remove whitespace
stats3$Mods3 <- gsub(' ', '', stats3$Mods3)

#now fix peptide sequence formatting
#going to just select for things between '.'s
stats3$Sequence2 <- rm_between(stats3$Sequence, '.', '.', extract = TRUE)

#pull out KSEA input
#note that it wants FC NOT log-transformed
KSEA_in <- stats3[c('Accession3', 'GeneSymbol', 'Sequence2', 'Mods3',
                    'pval', 'AR')]
colnames(KSEA_in) <- c('Protein', 'Gene', 'Peptide', 'Residue.Both', 'p', 'FC')

#remove lines with no gene name, or NA FC
KSEA_in <- KSEA_in[KSEA_in$Gene != '', ]
KSEA_in <- KSEA_in[!is.na(KSEA_in$FC), ]

write.csv(KSEA_in, 'KSEA_input.csv', row.names = FALSE, quote = FALSE)

#run KSEA app
shiny::runGitHub('KSEA', 'casecpb')


### Visualize results from KSEA App -----

kinase_results <- read.csv('KSEA_out_kinase_table.csv')

sig_kinases <- kinase_results[kinase_results$p.value < 0.1, ]
sig_kinases <- sig_kinases[order(sig_kinases$z.score, decreasing = FALSE), ]
sig_kinases$Kinase.Gene <- factor(sig_kinases$Kinase.Gene, levels = sig_kinases$Kinase.Gene)

pdf('all_patients_figures/kinase_barchart.pdf', width = 8, height = 3)
ggplot(sig_kinases, aes(x = Kinase.Gene, y = z.score, fill = -log10(p.value))) + 
  scale_fill_continuous(type = 'viridis', limits = c(0, 4.5)) + 
  geom_bar(stat = 'identity') + theme_bw() + labs(fill = '-log10(p-value)') & RotatedAxis()
dev.off()

sig_kinases2 <- sig_kinases[order(sig_kinases$mS, decreasing = FALSE), ]
sig_kinases2$Kinase.Gene <- factor(sig_kinases2$Kinase.Gene, levels = sig_kinases2$Kinase.Gene)

pdf('all_patients_figures/kinase_barchart_mS.pdf', width = 8, height = 3)
ggplot(sig_kinases2, aes(x = Kinase.Gene, y = mS, fill = -log10(p.value))) + 
  scale_fill_continuous(type = 'viridis', limits = c(0, 4.5)) + 
  geom_bar(stat = 'identity') + theme_bw() + labs(fill = '-log10(p-value)') & RotatedAxis()
dev.off()


