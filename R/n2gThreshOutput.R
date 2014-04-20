# This script aims to:
# Read in file with 3702 unique IDs for each sample
# For each concept map, extract sample sets for specified distance (3mm)

# Get list of .tab files in the directory
files = list.files(path = "/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/output/", pattern = "*.tab")
# Read in all files, and save UIDS and distances for
suids = read.csv('/home/vsochat/SCRIPT/python/neuro2gene/data/sampleUID.txt',sep=',',head=FALSE)

# HERE IS CODE TO MANUALLY ADD SAMPLE IDS TO TERM FILES with NO THRESHOLDING, AND ALSO WITH THRESHOLDING
for (f in 1:length(files)) {
  
  cat("Processing file",files[f],"of",length(files),"\n")
  file = files[f]
  
  # Find the term
  termid = strsplit(file,"_")[[1]][1]
  savepoints = c()
  sampids = c()
  
  # Read in the file
  file = read.csv(paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/output/",files[f],sep=""),sep='\t',head=TRUE)
  file = file[,1:18]
  
  for (s in 1:dim(file)[1]) {
    # Find the sample uid
    samp = as.character(paste(file[s,1],file[s,2],file[s,3],file[s,4],file[s,5],file[s,6],file[s,8],file[s,9],file[s,10],file[s,11],file[s,12],file[s,13],file[s,14],sep="-"))
    sidx = which(suids[,2] == samp)
    # Save to result
    sampy = as.character(suids[sidx,1])
    # If squared distance is <= 9 save point
    savepoints = rbind(savepoints,file[s,])   
    sampids = c(sampids,sampy)
  }

  rownames(savepoints) = sampids
  write.table(savepoints, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/thresh/nothresh3000/",termid,".tab",sep=""), append = FALSE, quote = TRUE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

  # Now threshold data, and save to different file
  savepoints = savepoints[which(savepoints$sq_distance <= 9),]
  write.table(savepoints, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/thresh/9mmsq3000/",termid,".tab",sep=""), append = FALSE, quote = TRUE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)  
}


