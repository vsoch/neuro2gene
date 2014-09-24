# This script will calculate regional overlap
setwd('/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/thresh/9mmsq')

# Read in file with results at threshold .1
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/disorderGenesN8pt1.Rda")
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/disorderGeneListsN8pt1.Rda")

# For each disorder, save long list of regions
disorders = unique(filter$DISORDER)
samples = read.csv("/home/vanessa/Documents/Work/ALLEN/Samples.csv",sep=",")
# Make a matrix of locations by disorder
regions = as.character(unique(samples$structure_name))
matrix = matrix(0,nrow=length(disorders),ncol=length(regions))
rownames(matrix) = disorders
colnames(matrix) = regions
for (i in 1:length(disorders)){
  d = as.character(disorders[i])
  cat(d,"\n")
  tmp = filter[which(filter$DISORDER==d),]  
  tmp2 = genes[which(genes$DISORDER == d),]
  # Get unique brainterms associated with disorder
  termsO = as.character(tmp2$TERM)
  terms = gsub("_down","",gsub("_up","",tolower(termsO)))
  # Get spatial locations now!
  # We will save a list of the terms and associated genes
  gen = list()
  cat("Generating gene lists for",d,"\n")
  for (t in unique(termsO)){
    # Add list of genes for the term
    tmp3 = tmp2[which(tmp2$TERM == t),]
    gen[t] = list(as.character(unique(tmp3$PROBE))) # Add gene to our list
  }
  
  # Calculate pairwise tanimoto scores 
  genmatrix = array(0,dim=c(length(gen),length(gen)))
  cat("Calculating tanimoto scores for",d,"\n")
  for (g in 1: length(gen)){
    for (t in 1:length(gen)){
     genmatrix[g,t] = length(intersect(gen[g],gen[t])) / length(union(gen[g],gen[t])) 
    }
  }
  diag(genmatrix) = 0
  idx = which(genmatrix!=0,arr.ind=TRUE)
  # This index also applies to the unique term list
  finalterms = unique(termsO)
  indextoremove = c()
  if (length(idx)>0){
    # Remove redundant gene lists from matrix
    idx = idx[1:(dim(idx)[1]/2),1]
    for (j in idx){
      indextoremove = c(indextoremove,j)    
    }
    if (length(!is.na(indextoremove))!=0) {
      finalterms = finalterms[-indextoremove]
    }
  }
  # Generate final structure counts
  for (t in finalterms){
    cat("adding",t,"\n")
    te = gsub("_down","",gsub("_up","",tolower(t)))
    dat = read.csv(list.files(pattern=paste("_",te,".tab",sep="")),sep="\t",head=TRUE)
    structures = as.character(dat$structure_name)
    structures = table(structures)
    matrix[d,names(structures)] = (matrix[d,names(structures)] + structures)    
  }
}
 
# Save matrix of region counts
setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults")
save(matrix,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/regionCountMatrixpt1.Rda")

# Create a matrix of disorder by term with NES scores
terms = sort(unique(genes$TERM))
disordermatrix = array(0,dim=c(length(disorders),length(terms)))
rownames(disordermatrix) = sort(disorders)
colnames(disordermatrix) = terms
for (d in disorders){
  subset = genes[which(genes$DISORDER==as.character(d)),]
  for (t in unique(subset$TERM)){
    disordermatrix[as.character(d),as.character(t)] = mean(subset$NES[which(subset$TERM == t)])
  }
}
save(disordermatrix,file="DisorderTermMatrixNESpt1.Rda")

# Create a matrix of disorder by genes with 0/1 values
gen = sort(unique(genes$PROBE))
genematrix = array(0,dim=c(length(disorders),length(gen)))
rownames(genematrix) = sort(disorders)
colnames(genematrix) = gen
for (d in disorders){
  tmp = unique(as.character(genes$PROBE[which(genes$DISORDER==d)]))
  genematrix[as.character(d),tmp] = 1
}
save(genematrix,file="DisorderGeneMatrixNESpt1.Rda")


# Now we need to get the MNI coordinates for each region
mni = list()
for (c in colnames(matrix)){
  # We will take a middle coordinate for the region
  mni[c] = list(colMeans(samples[which(samples$structure_name == c),c(12,13,14)]))  
}

# Normalize rows between 0 and 1
mnorm = t(apply(matrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
save(mnorm,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/regionNormMatrixpt1.Rda")

# Define range of colors
library('RColorBrewer')
rbPal <- colorRampPalette(brewer.pal(8,"YlOrRd"))

library("scatterplot3d")
options(digits=3)
# Now for each disorder, plot the regions and color by strength
for (d in rownames(mnorm)){
  png(paste(d,"_sigEnrichmentMap.png",sep=""),width=7,height=6,units="in",res=300)
  structures = names(which(mnorm[d,]!=0))
  values = as.numeric(mnorm[d,which(mnorm[d,]!=0)])
  counts = as.numeric(matrix[d,which(matrix[d,]!=0)])
  values = c(0,values,1)
  color = rbPal(10)[as.numeric(cut(values,breaks = 10))]
  color = color[-c(1,length(color))]
  values = values[-c(1,length(values))]
  tmp = as.data.frame(mni[structures])
  scatterplot3d(x=as.numeric(tmp[1,]),y=as.numeric(tmp[2,]),z=as.numeric(tmp[3,]),color=color,pch=19,xlab="mniX",ylab="mniY",zlab="mniZ",main=paste("Regional Centers with Significant Term Enrichment",d)) 
  dev.off()
  # Calculate number of significant regions as percentage of total
  perbraincoverage = rep(length(unique(structures))/length(unique(samples$structure_name)),length(values))
  # Print text file with list of values and coordinates
  result = as.data.frame(cbind(round(as.numeric(tmp[1,]),3),round(as.numeric(tmp[2,]),3),round(as.numeric(tmp[3,]),3),round(values,3),color,structures,counts,perbraincoverage))
  colnames(result) = c("MNIX","MNIY","MNIZ","VALUES","COLOR","STRUCTURES","COUNTS","PERBRAINCOVERAGE")
  save(result,file=paste(d,"_SigRegionEnrichment.Rda",sep=""))
  write.table(result,file=paste(d,"_SigRegionEnrichment.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
}

# Write tables for norm values and counts
write.table(round(mnorm,3),file="/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/regionNormMatrixpt1.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
write.table(matrix,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/regionCountMatrixpt1.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

# Try sparse hierarchical clustering
groups = c()
groups[[1]] = c("ALZ","MS")
groups[[2]] = c("POMPE","PTSD")
groups[[3]] = c("SZO","LUPUS","AIDS","CYTOVIRUS")
groups[[4]] = c("RETT","ASD")
groups[[5]] = c("MD","BRA","CRANIO")
groups[[6]] = c("ALC","PARK")

load("C:\\Users\\Vanessa\\Documents\\Work\\GENE_EXPRESSION\\neurosynth\\sigresult\\regionNormMatrixpt1.Rda")

for (g in 1:length(groups)){
  subset = mnorm[groups[[g]],]
  perm.out <- HierarchicalSparseCluster.permute(as.matrix(mnorm),nperms=100)
  plot(perm.out)
  sparsehc = HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw, method="complete")
  plot(sparsehc$hc, labels=rep("", nrow(mnorm)))
  print(sparsehc)

# Read in nifti image of MNI template
library(Rniftilib)
library("AnalyzeFMRI")
nii = nifti.image.read("/var/www/gbm/data/X.nii",read_data=1)
f.read.nifti.header(system.file("/var/www/gbm/data/X.nii", package="AnalyzeFMRI"))
save(mni,file="mniCoordinatesAllenStructures.Rda")

# Convert from coordinate space to MNI space
mnidf = as.data.frame(mni)
hdr = f.read.header("X.hdr")
mnicoord = xyz2ijk(xyz=mnidf,method=2,hdr)
mnicoord = mnicoord$ijk

# Read and empty the nifti file, fill coordinates with values
mr = f.read.nifti.volume("X.nii")
mr[,,,1] = 0

# Try writing region centers to an MNI image file
for (d in rownames(mnorm)){
  structures = names(which(mnorm[d,]!=0))
  values = as.numeric(mnorm[d,which(mnorm[d,]!=0)])
  counts = as.numeric(matrix[d,which(matrix[d,]!=0)])
  regions = names(mni[structures])
  # Here we fill empty image with values
  for (i in 1:length(regions)){
    r = regions[i]
    idx = which(names(mni)==r)
    coord = round(mnicoord[,idx])
    mr[coord[1],coord[2],coord[3],1] = values[i]
  }
  # Print new image to file
  f.write.nifti(mr,paste(d,"_sigRegions",sep=""),size="float",nii=TRUE)
}

# Finally, create a plot of unique regions by disorder
library(gplots)
library(ggplot2)
library(reshape)
# Get rid of regions with 0 significant
zerosig = names(which(colSums(mnorm)==0))
cat(zerosig,file="ZeroSigRegionspt1.txt",sep="\n")
zerosig = which(colSums(mnorm)==0)
hm = heatmap.2(mnorm[,-as.numeric(zerosig)],dendrogram=c("none"),key=FALSE,col="heat.colors",labCol=NULL)
tmp = mnorm[,-as.numeric(zerosig)]
tmp = tmp[hm$rowInd,hm$colInd]
tmp = melt(tmp)


# Make ggplot
png(file="disorderRegionHeatmap.png",width=18,height=6,units="in",res=300)
(p = ggplot(data = tmp, aes(X2, X1)) 
      + geom_tile(aes(fill = value),colour = "white") 
      + scale_fill_gradient(low = "white",high = "steelblue")
 + scale_x_discrete(expand = c(0, 0))
 + scale_y_discrete(expand = c(0, 0))
      + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 9 *0.8, angle = 90, hjust = 1, colour = "grey50"),axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 90, hjust = 1)))
dev.off()

# Attempt Number 2 at ggplot
library(reshape2)
library(ggdendro)

x = mnorm[,-as.numeric(zerosig)]
dd.col <- as.dendrogram(hclust(dist(x)))
col.ord <- order.dendrogram(dd.col)
dd.row <- as.dendrogram(hclust(dist(t(x))))
row.ord <- order.dendrogram(dd.row)
xx <- x[col.ord, row.ord]
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]
df$disorder <- xx_names[[1]]
df$disorder = with(df, factor(disorder, levels=disorder, ordered=TRUE))

mdf = melt(df, id.vars="disorder")
ddata_x <- dendro_data(dd.row)
ddata_y <- dendro_data(dd.col)

### Set up a blank theme
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

### Create plot components ###    
# Heatmap
png(file="disorderRegionHeatmapTall.png",width=6,height=18,units="in",res=300)
p1 <- ggplot(mdf, aes(x=disorder, y=variable)) + 
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low="red",mid="white",high="blue") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 9 *0.8, angle = 90, hjust = 1, colour = "grey50"),axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p1
dev.off()


# TAKE TWO - make plot of tanimoto scores by disorder
setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/")
load('DisorderRegionTanimotoMatrixpt1.Rda')

library(gplots)
library(ggplot2)
library(reshape)
# Get rid of regions with 0 significant
hm = heatmap.2(tanimotos,dendrogram=c("none"),key=FALSE,col="heat.colors",labCol=NULL)

library(ggdendro)
library(gplots)
library(ggplot2)
library(reshape)

x = tanimotogenes
dd.col <- as.dendrogram(hclust(dist(x)))
col.ord <- order.dendrogram(dd.col)
dd.row <- as.dendrogram(hclust(dist(t(x))))
row.ord <- order.dendrogram(dd.row)
xx <- x[col.ord, row.ord]
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]
df$disorder <- xx_names[[1]]
df$disorder = with(df, factor(disorder, levels=disorder, ordered=TRUE))

mdf = melt(df, id.vars="disorder")
ddata_x <- dendro_data(dd.row)
ddata_y <- dendro_data(dd.col)

### Set up a blank theme
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

### Create plot components ###    
library('RColorBrewer')

# Heatmap
png(file="disorderGenesTanimotos.png",width=8,height=7,units="in",res=300)
p1 <- ggplot(mdf, aes(x=disorder, y=variable)) + 
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low="red",mid="white",high="blue") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 9 *0.8, angle = 90, hjust = 1, colour = "grey50"),axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p1
dev.off()

# Define groups based on tree clustering
load("DisorderRegionHeatmappt1.Rda")
plot(hclust(dist(t(hmresult$matrix))))
gid = c()
#gid["LUPUS"] = 3
gid["SZO"] = 2
#gid["ALZ"] = 4
#gid["MS"] = 5
gid["PTSD"] = 2
gid["ASD"] = 1
gid["PARK"] = 2
gid["RETT"] = 1
matrix = mnorm
matrix = cbind(matrix,gid)
tree = rpart(gid ~.,data=as.data.frame(matrix))
plot(tree)

load("regionNormMatrixpt1.Rda")
matrix = mnorm
# Can we just look at regions in each group ~= 0?
for (g in unique(gid)){
 subset = matrix[names(which(gid==g)),]
 cat("\n\n",names(which(gid==g)),"\n")
 cat(names(which(colSums(subset)!=0)),sep="\n")
}

# Build tree with sparse Hierarchical clustering to do feature selection
library('sparcl')
# Here is the filtered matrix with region features, 1 == has, 0 == doesn't
load("regionNormMatrixpt1.Rda")
matrix = mnorm
perm.out <- HierarchicalSparseCluster.permute(matrix, wbounds=c(1.5,2:6),nperms=5)
print(perm.out)
plot(perm.out)
sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw, method="complete")
print(sparsehc)
par(mfrow=c(1,1))
plot(sparsehc$hc, labels=rownames(matrix))

# Let's try clustering by regions
load("regionNormMatrixpt1.Rda")
filt = mnorm[,-which(colSums(mnorm) == 0)]
disty = dist(filt)
hc = hclust(disty)
plot(hc)

# Can we go back to brain samples?
load("DisorderSampleTanimotoMatrixpt1.Rda")
sampsim
hms = heatmap.2(sampsim,dendrogram=c("none"),key=FALSE,col="heat.colors",labCol=NULL)

