# libraries ----
library(mice)
library(lattice)
library(readxl)  
library(VIM)
library(ggplot2)
library(tidyr)
library(dplyr)

# load ----
  setwd ("C:/Users/jackie/OneDrive - Dartmouth College/Mice_for_data_imputing/") 
  metabolomics <- read.csv("082020RajohnsonMet1_9_.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# clean ----
  colnames(metabolomics) <- c('Baseline_1','Baseline_2','Baseline_3','Baseline_4', 'Baseline_5', 
                              'Dormant_1', 'Dormant_2', 'Dormant_3', 'Dormant_4')
  meta <- metabolomics[-c(1,2),-c(10:12)] # remove unwanted rows, cols
  head(meta)
  
  meta_clean <- as.data.frame(apply(meta,2, as.numeric)) # change str to num 
  str(meta_clean)
  
  meta_names <- row.names(meta) # create vector with meta row.names
  row.names(meta_clean) <- meta_names # re-introduce to meta_clean
  head(meta_clean)
  
  meta_clean_t <- as.data.frame(t(meta_clean)) # transpose 
  head(meta_clean_t)
  str(meta_clean_t)
  class(meta_clean_t)

# Missing data ----
  pMiss <- function(meta_clean_t){sum(is.na(meta_clean_t))/length(meta_clean_t)*100} # % missing data per sample/rows
  apply(meta_clean_t,1,pMiss) # consider removing samples with >30% missing data
  meta_clean_t_filt <- meta_clean_t[ , colSums(is.na(meta_clean_t)) < nrow(meta_clean_t)] # remove metabolites/cols with all NAs 
  aggr_plot <- aggr(meta_clean_t_filt, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, # visualize missing data
                    labels=names(meta_clean_t), cex.axis=.7, gap=3, 
                    ylab=c("Histogram of missing data","Pattern"))
  
  meta_clean_t_filt_IMP <- mice(meta_clean_t_filt, maxit = 10, print=F) # impute missing data using MICE
  colnames(meta_clean_t_filt) <- gsub(" ", "_", colnames(meta_clean_t_filt))

# Transferred back to Anneka to fix metabolite/column names...

