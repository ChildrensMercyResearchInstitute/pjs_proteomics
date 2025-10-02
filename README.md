# Introduction 
Analysis code for global proteomic and phosphoproteomic comparison of polyp and mucosa samples from Children's Mercy PJS patients. 

# Inputs
There are two overlapping datasets explored in this study. The first dataset (presented in Figure 1) consists of data from two patients with samples collected from the stomach as well as small intestine. Inputs files are:
1.	2patients_proteo.csv (proteomic data)
2.	2patients_phospho.csv (phosphoproteomic data)

The rest of the paper focuses on small intestine samples only from all five patients. Inputs for this analysis are:
1.	sample_metadata.csv (manifest)
2.	all_patients_proteo.csv (proteomic data)
3.  all_patients_phospho.csv (phosphoproteomic data)

# Analysis code
There are three files containing R code for analysis of these datasets:
- 2patients_figures.R includes all analysis for the dataset presented in Figure 1
- all_patients_figures.R includes analysis for the 5-patient small intestine dataset and generates the panels in Figure 2
- all_patients_cytoscape.R is used to combine the proteomic and phosphoproteomic data into network structures readble with cytoscape and generates the panels in Figures 3 and 4

Please note that all_patients_cytoscape.R does require files generated in all_patients_figures.R
