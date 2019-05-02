#! /usr/bin/env Rscript

options(stringsAsFactors = F)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")

library(ShortRead, quietly = T, warn.conflicts = F)
library(BiocParallel, quietly = T, warn.conflicts = F)
library(Rsamtools, quietly = T, warn.conflicts = F)
library(GenomicRanges, quietly = T, warn.conflicts = F)
library(Biostrings, quietly = T, warn.conflicts = F)

print("Packages Loaded")

workingdir = commandArgs(TRUE)

setwd(workingdir)
getwd()

#### load in your bam files #### 
bams = list.files(pattern = ".bam$")

#### make an index (.bai) file for it and define some shit ####
bamidx = indexBam(bams)
bamfiles = BamFileList(bams, bamidx)
name = gsub("-R1.mapping.sorted.bam", "", bams)
nucleotides = seq(1115, 10535, 1) #### 1115 is the start of the first amplicon ####

#### write a csv for the nucleotide counts for each bam file - to get frequencies insteaad, change as.prob to T  #### 
NTFreqs = list()
for(h in 1:length(bamfiles)){
  #### check that the csv doesn't already exist ####
  csvfile = paste0("nucleotideCounts_", name[h], ".csv")
  if(!file.exists(csvfile)){
   print(paste0("Beginning ", name[h]))
    bamfile = bamfiles[[h]]
  
    NTFreqs.count = data.frame()
    for(i in 1:length(nucleotides)){
      data.c = nucleotideFrequencyAt(stackStringsFromBam(bamfile, param = paste0("M33262: ", nucleotides[i],"-", nucleotides[i])), at = 1, as.prob = F)
      data.df.c = as.data.frame(cbind(nucleotides[i], t(data.c), sum(data.c)))
      colnames(data.df.c) = c("Position", "A", "C", "G", "T", "Coverage")
      NTFreqs.count = rbind(NTFreqs.count, data.df.c)
    
     #### add a progress printout ####
      if(i %% 100 == 0){
        print(paste0("Currently at: ", nucleotides[i], " in ", name[h]))
      }
    }
    NTFreqs[[h]] = NTFreqs.count
    print(paste0("Writing csv file"))
    write.csv(NTFreqs.count, paste0("nucleotideCounts_",name[h] , ".csv"))
    print(paste0("CSV done!")) 
  }else{
    print(paste0("CSV for ", name[h], " already exists"))
  }
}

