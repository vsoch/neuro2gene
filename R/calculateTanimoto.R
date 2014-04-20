#This will read in each of shapley up and down sets, and calculate tanimoto scores.

indir = '/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/probeSets/9mmsq525/list'
downs = list.files(indir,pattern="*_down.Rda")
ups = list.files(indir,pattern="*_up.Rda")

# First read terms into a long list
uplist = list()
downlist = list()

# First do up probes
for (f in 1:length(ups)) {
  file = paste(indir,"/",ups[f],sep="")
  load(file)
  uplist = c(uplist,list(resultup$probes))
}

# Now do down probes
for (f in 1:length(downs)) {
  file = paste(indir,"/",downs[f],sep="")
  load(file)
  downlist = c(downlist,list(resultdown$probes))
}

# NOW TANIMOTO SCORES
TSdown = array(dim=c(length(downs),length(downs)))
TSup = array(dim=c(length(ups),length(ups)))

# Now calculate tanimoto scores for upsets
for (u in 1:length(ups)) {
  cat ("Calculating Tanimotos for Term",u,"of",length(ups),"\n")
  # Now we assess overlap of sets - calculate tanimoto score, intersection / union
  g1 = as.character(uplist[[u]])
  tanimoto = c()  # This will be vector of scores for one value
  for (o in 1:length(ups)) {
    g2 = as.character(uplist[[o]])
    inter = length(which(g1 %in% g2))
    uni = length(unique(c(g1,g2)))
    score = inter/uni
    tanimoto = c(tanimoto,score)
  }
  TSup[u,] = tanimoto
}
  
rownames(TSup) = ups
colnames(TSup) = ups
write.table(TSup,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/up525.csv")
save(TSup,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/up525.Rda")

# Now calculate tanimoto scores for downsets
for (u in 1:length(downs)) {
  cat ("Calculating Tanimotos for Term",u,"of",length(downs),"\n")
  # Now we assess overlap of sets - calculate tanimoto score, intersection / union
  g1 = as.character(downlist[[u]])
  tanimoto = c()  # This will be vector of scores for one value
  for (o in 1:length(downs)) {
    g2 = as.character(downlist[[o]])
    inter = length(which(g1 %in% g2))
    uni = length(unique(c(g1,g2)))
    score = inter/uni
    tanimoto = c(tanimoto,score)
  }
  TSdown[u,] = tanimoto
}

rownames(TSdown) = downs
colnames(TSdown) = downs
write.table(TSdown,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/down525.csv")
save(TSdown,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/down525.Rda")

alls = list.files(indir,pattern="*.Rda") 
alist = list()

for (f in 1:length(alls)) {
  file = paste(indir,"/",alls[f],sep="")
  load(file)
  if (grepl("_down.Rda",file)) {
     alist = c(alist,list(resultdown$probes))
  } else {
    alist = c(alist,list(resultup$probes))
 }
}

TSall = array(dim=c(length(alls),length(alls)))

# Now calculate tanimoto scores for ALL sets
for (u in 1:length(alls)) {
  cat ("Calculating Tanimotos for Term",u,"of",length(alls),"\n")
  # Now we assess overlap of sets - calculate tanimoto score, intersection / union
  g1 = as.character(alist[[u]])
  tanimoto = c()  # This will be vector of scores for one value
  for (o in 1:length(alls)) {
    g2 = as.character(alist[[o]])
    inter = length(which(g1 %in% g2))
    uni = length(unique(c(g1,g2)))
    score = inter/uni
    tanimoto = c(tanimoto,score)
  }
  TSall[u,] = tanimoto
}

rownames(TSall) = alls
colnames(TSall) = alls
write.table(TSall,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/all525.csv")
save(TSall,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/all525.Rda")

# Here we calculate the size of each set for each of up, down, and all, and print to file
asize = c()
for (s in 1:length(alls)){
  asize = c(asize,length(alist[[s]]))
}
dsize = c()
for (s in 1:length(downs)){
  dsize = c(dsize,length(downlist[[s]]))
}
usize = c()
for (s in 1:length(ups)){
  usize = c(usize,length(uplist[[s]]))
}

# Covert to table
asize = as.table(asize)
dsize = as.table(dsize)
usize = as.table(usize)

# Add names of terms
rownames(asize) = alls
rownames(dsize) = downs
rownames(usize) = ups

# Save to file
write.table(asize,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/all_sizes_525.csv",quote=FALSE)
write.table(usize,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/up_sizes_525.csv",quote=FALSE)
write.table(dsize,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/tanimoto/down_sizes_525.csv",quote=FALSE)

