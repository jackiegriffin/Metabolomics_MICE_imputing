## load packages ----


install.packages(c("lattice", "readxl","VIM","ggplot2","tidyr","mice","tidyr","dplyr"))
library(mice)
library(lattice)
library(readxl)  
library(VIM)
library(ggplot2)
library(tidyr)
library(mice)
library(dplyr)
library(VIM)


## data processing ----

# set working directory
setwd ("C:/Users/jackie/OneDrive - Dartmouth College/Mice_for_data_imputing/") 

#Load data into R and name it metabolomics
metabolomics <- read.csv("082020RajohnsonMet1_9_.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# rename columns to match treatment replicates
colnames(metabolomics) <- c('Baseline_1','Baseline_2','Baseline_3','Baseline_4', 'Baseline_5',
                            'Dormant_1', 'Dormant_2', 'Dormant_3', 'Dormant_4')
head(metabolomics)
meta <- metabolomics[-c(1,2), ] # remove rows 1 and 2 
meta <- meta[,-c(10:12) ] # remove columns 10 through 12
head(meta)

# clean data
meta_clean <- as.data.frame(apply(meta,2, as.numeric)) # change all cols from chr to num
str(meta_clean)
head(meta_clean)

meta_names <- row.names(meta) # make a vector with row names from meta df
row.names(meta_clean) <- meta_names # re-introduce row names to clean df
head(meta_clean)

# transpose - samples are row ids, features as column ids
meta_clean_t <- as.data.frame(t(meta_clean))
class(meta_clean_t)
str(meta_clean_t)

# Check for samples (rows) with > 30% missing data, if any, consider removing from analysis
pMiss <- function(meta_clean_t){sum(is.na(meta_clean_t))/length(meta_clean_t)*100}
apply(meta_clean_t,1,pMiss)

# visualize missing data
aggr_plot <- aggr(meta_clean_t, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(meta_clean_t), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

# remove features (columns) with all NAs 
meta_clean_t_filt <- meta_clean_t[ , colSums(is.na(meta_clean_t)) < nrow(meta_clean_t)]

# re-plot
aggr_plot <- aggr(meta_clean_t_filt, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(meta_clean_t), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

## MICE ----
imp <- mice(meta_clean_t_filt, maxit = 10, print=F)

## Need to remove spaces from column names

colnames(meta_clean_t_filt) <- gsub(" ", "_", colnames(meta_clean_t_filt))


## column names need cleaned up..