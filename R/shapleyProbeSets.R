#!/home/vsochat/bin/Rscript
# Generate probe sets
# This script will generate probe 
# File is path to a set of samples for a given term
# Thresh is the threshold for corrected p values
# Outfile is where to write the filtered raw expression data 
# Outlist is where to write the list of probes (with probeID and p-value)
# b is the number of bootstrap samples (1000)

args <- commandArgs(TRUE)
file = args[1]
thresh = as.numeric(args[2])
outfile = args[3]
outlist = args[4]
b = as.numeric(args[5])

cat("\n",file,":file\n")
cat(thresh,":thresh\n")
cat(outfile,":outfile\n")
cat(outlist,":list\n")
cat(b,":b\n\n")

  # First read in me tables
  me = list.files(path = "/scratch/users/vsochat/DATA/ALLEN/me/", pattern = "*_id.csv")
  
  # Generate column names (combining the files)
  colly = list()
  # Read in colnames for each me file
  for (p in 1:length(me)){
    tmp = read.csv(paste('/scratch/users/vsochat/DATA/ALLEN/me/',me[p],sep=""),sep=',',head=TRUE,nrows=1)
    tmp = list(colnames(tmp))
    colly = c(colly,tmp)
  }  
  
  # For our term, read in file, get gene expression for each sample
  file = read.csv(file,sep='\t',head=TRUE)
  probes = c()
  sampy = as.numeric(gsub("SAMP","",rownames(file[,])))  
  # Get gene expression from each me file (0-11)             
  for (p in 1:length(me)){
    tmp = read.csv(paste('/scratch/users/vsochat/DATA/ALLEN/me/',me[p],sep=""),sep=',',head=TRUE)
    colnames(tmp) = colly[[p]]
    rn = rownames(file)
    # Subset to the samples for the term
    tmp = tmp[sampy,16:dim(tmp)[2]]
    tmp =sapply(tmp,as.numeric)
    probes = cbind(probes,tmp)
  }
  rownames(probes) = rn
  
  # Now we have a matrix of s samples by n probes, and we want to find the "under expressed" and "overexpressed"
  means = apply(probes, 2, mean)
  sd = apply(probes, 2, sd)
  upperlim = means + sd
  lowerlim = means - sd

  # This cuts out the middle between -1 and 1 standard deviation
  B1 <- probes 
  B1[probes < upperlim] <- 0 
  B1[probes >= upperlim] <- 1 
  
  B2 <- probes 
  B2[probes <= lowerlim] <- 1 
  B2[probes > lowerlim] <- 0 

  # Now do Shapley Relative Importance Analysis
  # we have a matrix with probes in columns, samples in rows - let's calculate the shapley matrix, adjust with bootstrap
  # R1 is shapley values for "relatively up" probes
  # R2 is shapley values for "relatively down" probes
  source('/home/vsochat/SCRIPT/R/neurosynth/shapleymat.R')
  R1 <- shapleymat(B1)
  R2 <- shapleymat(B2)
  
  # Run mutliple hypothesis testing for "up" set (R1)
  n = dim(R1)[2]
  k = dim(R1)[1]  
  R <- matrix(NA, n, k)
  R[,1:k] <- R1
  library(multtest)
  OUTPUTup <- MTP(X = R, standardize = FALSE, B = b, get.adjp = TRUE, test = "t.onesamp")

  # Run mutliple hypothesis testing for "down" set (R2)
  n = dim(R2)[2]
  k = dim(R2)[1]  # K and H are equal due to same number of samples
  R <- matrix(NA, n, k)
  R[, 1:k] <- R2
  OUTPUTdown <- MTP(X = R, standardize = FALSE, B = b, get.adjp = TRUE, test = "t.onesamp")

  # Significant UP probes
  # Print raw data to file (probeMatrix), and probes sets (names and pvalues)
  resultup = as.data.frame(cbind(as.character(colnames(probes)[which(OUTPUTup@adjp < thresh)]),OUTPUTup@adjp[which(OUTPUTup@adjp < thresh)]))
  colnames(resultup) = c("probes","pval","rawp")
  
  # Significant DOWN probes
  # Print raw data to file (probeMatrix), and probes sets (names and pvalues)
  resultdown = as.data.frame(cbind(as.character(colnames(probes)[which(OUTPUTdown@adjp < thresh)]),OUTPUTdown@adjp[which(OUTPUTdown@adjp < thresh)]))
  colnames(resultdown) = c("probes","pval")

  # Write UP to output file  
  save(resultup,file=paste(outlist,"_up.Rda",sep=""))
  mat = probes[,which(OUTPUTup@adjp < thresh)]
  save(mat,file=paste(outfile,"_MEdata_up.Rda",sep=""))

  # Write DOWN to output file
  save(resultdown,file=paste(outlist,"_down.Rda",sep=""))
  mat = probes[,which(OUTPUTdown@adjp < thresh)]
  save(mat,file=paste(outfile,"_MEdata_down.Rda",sep=""))
