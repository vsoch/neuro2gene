# This script aims to:
# Read in file with 3702 unique IDs for each sample
# For each concept map, extract sample sets for specified distance (3mm)

# Get list of .tab files in the directory
files = list.files(path = "/home/vsochat/SCRIPT/python/neurosynth/output", pattern = "*.tab")
# Read in unique IDs for terms
uids = read.csv('/home/vsochat/SCRIPT/python/neurosynth/data/featureUID.txt',sep=',',head=FALSE)
# Read in all files, and save UIDS and distances for
suids = read.csv('/home/vsochat/SCRIPT/python/neurosynth/data/sampleUID.txt',sep=',',head=FALSE)

# We want to keep track of term, term UID, sample ID, and distance
tmatrix = matrix(0,dim(uids)[1],dim(suids)[1])
rownames(tmatrix) = uids[,1]
colnames(tmatrix) = suids[,1]

# HERE IS CODE TO MANUALLY ADD SAMPLE IDS TO TERM FILES with NO THRESHOLDING, AND ALSO WITH THRESHOLDING
for (f in 1:length(files)) {
  file = read.csv(paste("/home/vsochat/SCRIPT/python/neurosynth/output/",files[f],sep=""),sep='\t',head=TRUE)
  file = file[,1:18]

  # Find the term UID
  tidx = strsplit(strsplit(file,"__")[[1]][1],"/")
  tidx = as.character(tidx[[1]][length(tidx[[1]])])
  tidx = which(uids[,2]==tidx)
  termid = uids[tidx,1]
  termy = uids[tidx,2]
  savepoints = c()
  sampids = c()

  for (s in 1:dim(file)[1]) {
    # Find the sample uid
    samp = as.character(paste(f[s,1],f[s,2],f[s,3],f[s,4],f[s,5],f[s,6],f[s,8],f[s,9],f[s,10],f[s,11],f[s,12],f[s,13],f[s,14],sep="-"))
    sidx = which(suids[,2] == samp)
    # Save to result
    sampy = as.character(suids[sidx,1])
    # If squared distance is <= 9 save point
    savepoints = rbind(savepoints,f[s,])   
    sampids = c(sampids,sampy)
  }

  rownames(savepoints) = sampids
  write.table(savepoints, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/thresh/nothresh3000/",termid,"_",termy,".tab",sep=""), append = FALSE, quote = TRUE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

  # Now threshold data, and save to different file
  # VANESSA - NEED TO ADD THIS!
  write.table(savepoints, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/thresh/9mmsq3000/",termid,"_",termy,".tab",sep=""), append = FALSE, quote = TRUE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
  
}

