#! /usr/bin/env Rscript

# Rscript PATH_TO_SCRIPT.R MixingExperiment BaseRoot FileNameBase
options(stringsAsFactors = F)
filesneeded = commandArgs(TRUE)
Experiments = filesneeded[1]
BaseRoot = filesneeded[2]
filename = filesneeded[3]

#### comparing the 
vcfReader = function(Experiments, BaseRoot){
  #### add required packages ####
  require(vcfR)
  require(stringr)
  
  #### get to the fun stuff! #### 
  exps = Experiments
  path = BaseRoot
  vcfdata = list()
  for(a in 1:length(exps)){
    newpath = paste0(path, exps[a])
    print(newpath)
    setwd(newpath)
    tempvcfs = list.files(pattern = "annotated.vcf")
    tempvcfs.names = as.data.frame(paste0(exps[a], "_", tempvcfs))
    sample.df = data.frame()
    for(b in 1:length(tempvcfs)){
      sample = read.vcfR(tempvcfs[b])
      sample.name = gsub("-R1.mapping.annotated.vcf", "", tempvcfs.names[b,])
      tempdf = data.frame()
      for(c in 1:nrow(extract_gt_tidy(sample))){
        sample.snp = as.numeric(getFIX(sample)[c,"POS"])
        sample.freq = as.numeric(gsub("%", "", extract_gt_tidy(sample)[c, "gt_FREQ"]))
        sample.cov = as.numeric(extract_gt_tidy(sample)[c, "gt_SDP"])
        sample.ref.nt = getFIX(sample)[c, "REF"]
        sample.var.nt = getFIX(sample)[c, "ALT"]
        tempdf[c, "sample.name"] = sample.name
        tempdf[c, "sample.snp"] = sample.snp
        tempdf[c, "sample.freq"] = sample.freq
        tempdf[c, "sample.cov"] = sample.cov
        tempdf[c, "sample.ref.nt"] = sample.ref.nt
        tempdf[c, "sample.var.nt"] = sample.var.nt
      }
      sample.df = rbind(sample.df, tempdf)
    }
    vcfdata[[a]] = sample.df
  }
  
  vcfdata.df = data.frame()
  for(i in 1:length(vcfdata)){
    vcfdata.df = rbind(vcfdata[[i]], vcfdata.df)
  }
  return(vcfdata.df)
}

vcffiles = vcfReader(Experiments, BaseRoot)
write.csv(vcffiles, filename)


