# This script aims to:
# read in term thresholded results files
# For each set of samples - calculate a mean expression value, and save all to a matrix
#          gene1, gene2, gene3... geneN
# term 1
# term 2
# term 3

# We will then use this matrix to calculate correlations between expression... does that make sense?
datadir = "/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/probeSets/9mmsq525/data"

# Here is how to find a term file
#grep("*audiovisual_probeSet_MEdata_down.Rda",search)

# These are the terms with significant GSEA results - the data across all samples
# Angry Up, Noun Up, Skills Up, Covert Down
files = c("TRM397_angry_probeSet_MEdata_up.Rda","TRM261_noun_probeSet_MEdata_up.Rda","TRM043_colors_probeSet_MEdata_down.Rda","TRM086_audiovisual_probeSet_MEdata_down.Rda","TRM123_rejection_probeSet_MEdata_down.Rda")

# First get the master list of gene probes, for the final matrix
genenames = c()
for (f in 1:length(files)){
  load(paste(datadir,"/",files[f],sep=""))
  genenames = c(genenames,colnames(mat))
}
genenames = unique(genenames)

# STOPPED HERE - need to make this a data frame, and take mean over genes
matrixMeans = array(data=0,dim=c(length(files),length(genenames)))
matrixSums = array(data=0,dim=c(length(files),length(genenames)))
ngenes = c()
nsamps = c()

# Add row and column names
rownames(matrixMeans) = rownames(matrixSums) = files
colnames(matrixMeans) = colnames(matrixSums) = genenames

# Read in each sample by gene table, and get an average expression value
# Put them in our data matrices in the right spot
for (f in 1:length(files)){
  load(paste(datadir,"/",files[f],sep=""))
  numgenes = dim(mat)[2]
  numsamps = dim(mat)[1]
  ngenes = c(ngenes,numgenes)
  nsamps = c(nsamps,numsamps)
  matrixMeans[files[f],colnames(mat)] = colSums(mat) / numgenes
  matrixSums[files[f],colnames(mat)] = colSums(mat)
}

ngenes = as.data.frame(ngenes)
nsamps = as.data.frame(nsamps)

# Look up gene names
probes = read.csv("/scratch/users/vsochat/DATA/ALLEN/Probes.csv",head=FALSE,sep=",")
pids = gsub("pid_","",colnames(matrixMeans))
genes = probes$V4[which(probes$V1 %in% pids)]
tmp = t(matrixMeans)
tmp = as.data.frame(tmp)
for (t in 1:dim(tmp)[2]){
  tmp[,t] = as.numeric(tmp[,t])
}
colnames(tmp) = c("angry","noun","colors","audiovisual","rejection")

# Aggregate over genes
agg = aggregate(tmp,by=list(genes = genes),FUN=mean, na.rm=TRUE)
rownames(agg) = agg$genes
rownames(nsamps) = files
rownames(ngenes) = files
agg = agg[,-1]

# Now save to output file
write.table(agg, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult/termMatrixMeans_ASDSZOALZ.csv"), append = FALSE, quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(ngenes, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult/termNumGenes_ASDSZOALZ.csv"), append = FALSE, quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(nsamps, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult/termNumSamps_ASDSZOALZ.csv"), append = FALSE, quote = FALSE, row.names = TRUE, col.names = TRUE)

# Now calculate correlations for matrixMeans
corMeans = cor(matrixMeans,method="pearson")
write.table(corMeans, file = paste("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult/geneMeansCorrMatrix_n4.csv"), append = FALSE, quote = FALSE, row.names = TRUE, col.names = TRUE)
