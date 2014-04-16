# n2gGetExpression.R

# This script will read in a file produced by n2gExtractExpression.py
# (a tab separated file with ABA_ID, AFFY_ID, and GENE_ID) and extract
# gene expression from a local Allen Brain Atlas database

# Read in file with genes to get
genes = read.table("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/route66/route66ABA.txt",sep="\t",head=TRUE)

# This is the folder with the Allen Brain Atlas database files
database = "/home/vanessa/Documents/Work/ALLEN/R"

# Read lookup table to get expression file
lookup = read.table("/home/vanessa/Documents/Dropbox/Code/Python/neurosynth/data/MEindexTable.csv",sep=",",head=FALSE)
colnames(lookup) = c("ABA_ID","fid","ABA_ID_PID")

# For each gene, find the correct file, and lookup!
files = c()
for (g in 1:dim(genes)[1]){
  gene = genes[g,]
  idx = which(lookup$ABA_ID == gene$ABA_ID)
  files = c(files,as.character(lookup$fid[idx]))
}

# Now get expression and pacall data!
unis = unique(files)  
expression = array(dim=3702)
pcall = array(dim=3702)
ids = c()
for (u in unis){
  # Read in data
  cat("Processing ",which(unis %in% u),"of",length(unis),"\n")
  tmp = read.table(paste(database,"/",u,".csv",sep=""),sep=",",head=TRUE)
  pacall = read.table(paste(database,"/",gsub("me","pacall",u),".csv",sep=""),sep=",",head=TRUE)
  colnames(tmp) = gsub("pid_","",colnames(tmp))
  colnames(pacall) = gsub("pid_","",colnames(pacall))
  # Get expression for ids in here
  ids = c(ids,colnames(tmp)[which(colnames(tmp) %in% genes$ABA_ID)])
  tmp = tmp[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,which(colnames(tmp) %in% genes$ABA_ID))]
  pacall = pacall[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,which(colnames(pacall) %in% genes$ABA_ID))]
  suids = tmp[,c(1:14)]
  expression = cbind(expression,tmp[,15:dim(tmp)[2]])  
  pcall = cbind(pcall,pacall[,15:dim(pacall)[2]])  
}  

# Add column names
pcall = pcall[,-1]
expression = expression[,-1]
colnames(expression) = ids
colnames(pcall) = ids

# Put variables together
result = list(exp=expression,sig=pcall,rows=suids,genes=genes,lookup=lookup)

# Save expression to file
save(result,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/route66/route66.Rda")

# Then do something with expression data!

}