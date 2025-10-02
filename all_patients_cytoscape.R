library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(proDA)

### Read in data -----

## read in proteo data
proteo <- read.csv('inputs/all_patients_proteo.csv')
proteo <- proteo[1:4405, ]

proteo_stats <- proteo[c(3, 4, 9, 10)]
colnames(proteo_stats) <- c('Accession', 'Substrate.Gene', 'proteo_log2FC', 'proteo_pValue')

## read in kinase target data (from KSEA)

ksea <- read.csv('KSEA_out_kinase_targets.csv')


### Combining data for network visualization -----

#which kinases are not listed as targets?
kinase_list <- unique(ksea$Kinase.Gene)
length(kinase_list[!(kinase_list %in% ksea$Substrate.Gene)])
#123 kinases

#create df to add kinase genes that are not in Substrate.Gene
add_df <- data.frame(rep('', 123), kinase_list[!(kinase_list %in% ksea$Substrate.Gene)],
                     rep('', 123), rep('', 123), rep(0, 123))
colnames(add_df) <- colnames(ksea)

all_table <- rbind(ksea, add_df)

#add proteomics data
all_table <- merge(all_table, proteo_stats, all.x = TRUE)

write.csv(all_table, 'kinase_target_abundance.csv', row.names = FALSE, quote = FALSE)

#remove NAs in target nodes, which orphans them in the network
all_table0 <- all_table
all_table0$Accession <- case_when(
  is.na(all_table0$Accession) ~ all_table0$Substrate.Gene, 
  TRUE ~ all_table0$Accession
)
all_table0$proteo_log2FC <- case_when(
  is.na(all_table0$proteo_log2FC) ~ 0,
  TRUE ~ all_table0$proteo_log2FC
)
all_table0$proteo_pValue <- case_when(
  is.na(all_table0$proteo_pValue) ~ 1,
  TRUE ~ all_table0$proteo_pValue
)

write.csv(all_table0, 'kinase_target_abundance_noNA.csv', row.names = FALSE, quote = FALSE)

### Proteins of interest -----

polyp_poi <- c('YWHAE', 'SUB1', 'ASAP2', 'PNKP', 'CALU', 'DEK', 'TOP2B', 
               'EIF4G1', 'FIP1L1', 'VIM', 'JUNB', 'SSRP1', 'NOLC1', 
               'FLNA', 'IWS1')

mucosa_poi <- c('BAD', 'YAP1', 'CAV1', 'NPM1', 'RANBP2', 'PRKAR2A', 'CAMK2D')

polyp_network <- all_table0[all_table0$Substrate.Gene %in% polyp_poi, ]

mucosa_network <- all_table0[all_table0$Substrate.Gene %in% mucosa_poi, ]

write.csv(polyp_network, 'cytoscape/polyps.csv', row.names = FALSE, quote = FALSE)

write.csv(mucosa_network, 'cytoscape/mucosa.csv', row.names = FALSE, quote = FALSE)


### Plotting individual proteins and phosphopeptides -----

## read in sample metadata
ann_data <- read.csv('sample_metadata.csv')
ann_data$patient <- as.character(ann_data$patient)
rownames(ann_data) <- ann_data$sample


## read in proteo abundances
proteo_abun <- read.csv('proteo_norm_log_counts_bc.csv', row.names = 1)
colnames(proteo_abun) <- substring(colnames(proteo_abun), 2)

#because order of rows is same between proteo and proteo_abun
#use gene symbols for rownames for plotting
rownames(proteo_abun) <- make.unique(proteo$Gene.Symbol)

#transpose 
t_proteo <- as.data.frame(t(proteo_abun))

#annotate with sample information
ann_proteo <- merge(t_proteo, ann_data, by = 'row.names')
rownames(ann_proteo) <- ann_proteo$Row.names

for(protein in polyp_poi){
  if(protein %in% colnames(ann_proteo)){
    filename <- paste0('polyp_poi/proteo_', protein, '.pdf')
    pdf(filename, height = 4, width = 4)
    print(ggplot(ann_proteo, aes(x = tissue, y = .data[[protein]], fill = tissue)) +
            geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
            xlab(protein) + ylab('log(counts)'))
    dev.off()
  }
  else{
    print(paste0(protein, ' not in dataset'))
  }
}
#ASAP2, FIP1L1, JUNB not in dataset

for(protein in mucosa_poi){
  if(protein %in% colnames(ann_proteo)){
    filename <- paste0('mucosa_poi/proteo_', protein, '.pdf')
    pdf(filename, height = 4, width = 4)
    print(ggplot(ann_proteo, aes(x = tissue, y = .data[[protein]], fill = tissue)) +
            geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
            xlab(protein) + ylab('log(counts)'))
    dev.off()
  }
  else{
    print(paste0(protein, ' not in dataset'))
  }
}
#BAD, YAP1 not in dataset


## read in phospho data

phospho <- read.csv('inputs/all_patients_phospho.csv')
phospho <- phospho[1:1516, ]

## read in phospho abundances
phospho_abun <- read.csv('phospho_norm_log_counts_bc.csv', row.names = 1)
colnames(phospho_abun) <- substring(colnames(phospho_abun), 2)

#because order of rows is same between phospho and phospho_abun
#use gene symbols for rownames for plotting
rownames(phospho_abun) <- make.unique(phospho$Gene.Symbol)

#transpose 
t_phospho <- as.data.frame(t(phospho_abun))

#annotate with sample information
ann_phospho <- merge(t_phospho, ann_data, by = 'row.names')
rownames(ann_phospho) <- ann_phospho$Row.names


## pull out phosphopeptides for all proteins of interest
#because will be multiple per protein

stripped_colnames <- vapply(strsplit(as.character(colnames(ann_phospho)), '\\.'), '[', '', 1)

polyp_ann_phospho <- ann_phospho[ , stripped_colnames %in% c(polyp_poi, colnames(ann_data))]
mucosa_ann_phospho <- ann_phospho[ , stripped_colnames %in% c(mucosa_poi, colnames(ann_data))]

polyp_poi2 <- colnames(polyp_ann_phospho)[!(colnames(polyp_ann_phospho) %in% colnames(ann_data))]
mucosa_poi2 <- colnames(mucosa_ann_phospho)[!(colnames(mucosa_ann_phospho) %in% colnames(ann_data))]

for(protein in polyp_poi2){
  filename <- paste0('polyp_poi/phospho_', protein, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(ann_phospho, aes(x = tissue, y = .data[[protein]], fill = tissue)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(protein) + ylab('log(counts)'))
  dev.off()
}

for(protein in mucosa_poi2){
  filename <- paste0('mucosa_poi/phospho_', protein, '.pdf')
  pdf(filename, height = 4, width = 4)
  print(ggplot(ann_phospho, aes(x = tissue, y = .data[[protein]], fill = tissue)) +
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + theme_classic() + 
          xlab(protein) + ylab('log(counts)'))
  dev.off()
}



