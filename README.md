# This script will 
  - fit segmented models using genotype-specific coefficients reported in Schoppach et al. 2017
  - build summary data-sets
  - reproduce plots S1 and S2
  - summarize citation data obtailed from Web of Science (WoS)
  - plot citation data obtained from WoS
#
#
# Instructions for calculating accumulative water use from coefficients published in Schoppach et al. 2017
 1) put all files into a single directory (R script, data):
       - accumulative_t.R, climate_greeley.txt, climate_rf.txt, genotype_coef.csv
 2) open R script "accumulative_t.R" in Rstudio (or similar)
 3) load libraries (top 6 lines of code in script)
 4) Run code 
#
#
# Instruction for plotting data obtained from WoS (citation figures)
 1) put all files into a single directory (R script, data):
       - citation_data_summary_and_figures.R, p_all.rds
 2) open R script "citation_data_summary_and_figures .R" in Rstudio (or similar)
 3) Run code 
#
#
 last tested by SMG on December 5, 2024, Ubuntu Linux 20.04, R 3.6.3, RStudio 2024.04.2 Build #764
