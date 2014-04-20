# This script will read in gene sets, and create a "database" of gene sets to be used with GSEA.
# Read in file with gene_id lookup - we need to overlap genes in allen brain atlas with this data

# Write to this file
outfile = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/gsea/database/brainTerms3000.gmt",sep="")

# Read in file with probes
probes = read.csv("/scratch/users/vsochat/DATA/ALLEN/Probes.csv",sep=",",header=FALSE)

listdir = "/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/probeSets/9mmsq3000/list/"

# TO DO - create text file with list of 3000 terms
terms = read.csv("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/labels/features3000.txt",head=FALSE)
terms = terms$V1

# UP PROBES
sink(outfile)
for (t in 1:length(terms)) {
  term = terms[t]
    
  if (file.exists(paste(listdir,term,"_probeSet_up.Rda",sep=""))) {
    
    ptmp = probes
    
    # Load the up set
    load(file=paste(listdir,term,"_probeSet_up.Rda",sep=""))
    # Sort by pvalue significance
    resultup$pval = as.numeric(as.character(resultup$pval))
    resultup = resultup[with(resultup, order(pval)),]
    pid = as.numeric(gsub("pid_","",resultup$probes))
    # Find probes in probes table
    ptmp = ptmp[which(ptmp$V1 %in% pid),]
    # Get the gene list
    genes = unique(ptmp$V4)
    label = paste(as.character(term),"_up",sep="")
    cat(label,"Allen Brain Atlas expression up regulated subset for neurosynth fdr .05 corrected brain map sorted by pval more sig first",as.character(genes),"\n",sep="\t") 
    rm(resultup)
  }
  if (file.exists(paste(listdir,term,"_probeSet_down.Rda",sep=""))) {

    # Load the down set
    load(file=paste(listdir,term,"_probeSet_down.Rda",sep=""))
    ptmp = probes
    # Sort by pvalue significance
    resultdown$pval = as.numeric(as.character(resultdown$pval))
    resultdown = resultdown[with(resultdown, order(pval)),]
    pid = as.numeric(gsub("pid_","",resultdown$probes))
    # Find probes in probes table
    ptmp = ptmp[which(ptmp$V1 %in% pid),]
    # Get the gene list
    genes = unique(ptmp$V4)    
    label = paste(as.character(term),"_down",sep="")
    cat(label,"Allen Brain Atlas expression down regulated subset for neurosynth fdr .05 corrected brain map sorted by pval more sig first",as.character(genes),"\n",sep="\t") 
    rm(resultdown)
  }
}
sink()
