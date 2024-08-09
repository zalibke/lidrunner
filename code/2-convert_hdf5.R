##################
# Script should be used after 1-pull_hdf5.py to collapse ploidy and convert hdf5 files to rdata
# takes 1 argument POP_NAME, a unique population identifier
# rdata files are written to ../data/

# LOAD THE PACKAGE AND OTHERS
# install.packages("hdf5r")
library(hdf5r)
library(compiler)

#parse commandline args
args <- commandArgs(trailingOnly = TRUE)
POP_NAME <- args[1]
#OUTFOLDER <- args[2]

filename <- paste("../data/", POP_NAME, ".hdf5", sep="")
# MOUNT THE hdf5 FILE TO R, 'r' MEANS READ ONLY

dat<-H5File$new(filename, 'r')
dat

# EXTRACT and collapse genotype and POS vectors for all chr
sample_id<-dat[['sample_id']][]

genotype <- dat[['genotype_3R']][1,,] + dat[['genotype_3R']][2,,] #collapse ploidy
POS<-dat[['POS_3R']][]
outfile <- paste("../data/", "3R_", POP_NAME, ".Rdata", sep="")
save(genotype, POS, sample_id, file=outfile)
rm(POS, genotype, outfile)

genotype <- dat[['genotype_3L']][1,,] + dat[['genotype_3L']][2,,] #collapse ploidy
POS<-dat[['POS_3L']][]
outfile <- paste("../data/", "3L_", POP_NAME, ".Rdata", sep="")
save(genotype, POS, sample_id, file=outfile)
rm(POS, genotype, outfile)

genotype <- dat[['genotype_2R']][1,,] + dat[['genotype_2R']][2,,] #collapse ploidy
POS<-dat[['POS_2R']][]
outfile <- paste("../data/", "2R_", POP_NAME, ".Rdata", sep="")
save(genotype, POS, sample_id, file=outfile)
rm(POS, genotype, outfile)

genotype <- dat[['genotype_2L']][1,,] + dat[['genotype_2L']][2,,] #collapse ploidy
POS<-dat[['POS_2L']][]
outfile <- paste("../data/", "2L_", POP_NAME, ".Rdata", sep="")
save(genotype, POS, sample_id, file=outfile)
rm(POS, genotype, outfile)

genotype <- dat[['genotype_X']][1,,] + dat[['genotype_X']][2,,] #collapse ploidy
POS<-dat[['POS_X']][]
outfile <- paste("../data/", "X_", POP_NAME, ".Rdata", sep="")
save(genotype, POS, sample_id, file=outfile)
rm(POS, genotype, outfile)
