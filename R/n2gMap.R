# This script will map a neuro2gene result object to our self organizing map

library(shiny)
library('RColorBrewer')
library('png')

# Load the SOM
load('/home/vanessa/Documents/Dropbox/Code/R/shiny/brainMap/brainMap.Rda')      # This is the som and labels
rbPal <- colorRampPalette(brewer.pal(8,"YlOrRd"))

# Create a lookup table of indices for labels
lookup = c()
for (l in 1:length(brainMap$labels)){
  label = brainMap$labels[l]
  if (nchar(label) > 0){
    pieces = strsplit(label,"\n")[[1]]
    for (p in pieces){
      lookup[p] = l
    }
  }
}

save(lookup,file="/home/vanessa/Documents/Dropbox/Code/R/shiny/brainMap/lookup.Rda")

# Now create a picture for each disorder result
setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults")
load('disorderTermMatrixpt25.Rda')

for (d in 1:nrow(disordermatrix)){
  disorder = rownames(disordermatrix)[d]
  png(paste("/home/vanessa/Desktop/",disorder,".png",sep=""),width=1000,height=800)
  terms = tolower(names(which(disordermatrix[d,]!=0)))
  terms = gsub("_down","",gsub("_up","",terms))
  idx = lookup[terms]
  values = as.numeric(disordermatrix[d,which(disordermatrix[d,]!=0)])
  colors = array(0,dim=506)
  colors[idx] = rbPal(10)[as.numeric(cut(values,breaks = 10))]
  plot(brainMap$som$grid$pts,main=paste("BrainTerms Map for",disorder),xaxt='n', yaxt='n', col=colors,xlab="",ylab="",pch=15,cex=8)
  text(brainMap$som$grid$pts,brainMap$labels,cex=.6)   
  dev.off()
}