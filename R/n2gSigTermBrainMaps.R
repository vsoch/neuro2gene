setwd('/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult')
setwd('/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults')
setwd("C:\\Users\\Vanessa\\Documents\\Work\\GENE_EXPRESSION\\neurosynth\\sigresult")

load('DisorderTermMatrixNESpt1.Rda')

samplelist = list()
all = c()
# For each disorder
for (d in rownames(disordermatrix)){
  #mapdir = "/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/thresh/9mmsq525"
  mapdir = "/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/thresh/9mmsq"
  # What I'm trying:
  # 1) load all significant brain maps for a disease (first trying threshold .1)
  # 2) plot sample locations (and compare to other diseases)
  # 3) Do some kind of subset selection for locations (if too big)
  # 4) Compare against others again
  
  # 5) Save unique regions for all disorders
  # 6) Also do subset selection on regions
  # 7) Come up with final list?
  
  # 8) Also look at shared regions between disorders?
  
  # Find each term file and read in sample locations
  allsamples = c()
  
  maps = names(which(disordermatrix[d,] != 0))
  
  for (m in maps){
    map = strsplit(tolower(m),"_")[[1]][1]
    filey = list.files(mapdir,pattern=paste("_",map,".tab",sep=""),full.names=TRUE)            
    samples = read.csv(filey,sep="\t",head=TRUE)
    samples = cbind(samples,rep(m,nrow(samples)))
    allsamples = rbind(allsamples,samples)
  }
  sampids = rownames(allsamples)
  samplelist = c(samplelist,(list(allsamples)))
  all = rbind(all,cbind(sampids,allsamples,rep(d,nrow(allsamples))))
}

colnames(all)[1] = "SID"
colnames(all)[20] = "TERM"
colnames(all)[21] = "DISORDER"

# Now plot with a color for each disorder
library(scatterplot3d)

s3d = scatterplot3d(all$mni_x,all$mni_y,all$mni_z,pch=20,color=as.numeric(all$DISORDER),main="Overlap of Sample Locations from Brainterms Thresh .1",xlab="MNIX",ylab="MNIY",zlab="MNIZ")
legend(5,-1, col= unique(as.numeric(all$DISORDER)), bg="white", lty=c(1,1), lwd=2, yjust=0, legend = unique(all$DISORDER), cex = .7)

# Now let's use each location as a feature - calculate similarity
sampsim = array(0,dim=c(length(unique(all$SID)),length(unique(all$DISORDER))))
rownames(sampsim) = sort(unique(all$SID))
colnames(sampsim) = unique(all$DISORDER)
for (d in unique(all$DISORDER)){
  sampsim[as.character(all$SID[which(all$DISORDER == d)]),d] = sampsim[as.character(all$SID[which(all$DISORDER == d)]),d] + 1
}  

save(sampsim,file="DisorderSampleTanimotoMatrixpt1.Rda")
disty = dist(t(sampsim))
hc = hclust(disty)
plot(hc,main="Clustering Based on Sample Location, thresh .1")

#pdf('/home/vsochat/sampleheatmap.pdf',width=60,height=4)
library("gplots")
hm = heatmap.2(t(sampsim),dendrogram=c("none"),key=FALSE,col="heat.colors",labCol=NULL,main="Neuro2Gene Term Enrichment")
#dev.off()

# Also try lower resolution - region names
# Now let's use each location as a feature - calculate similarity
regsim = array(0,dim=c(length(unique(all$structure_name)),length(unique(all$DISORDER))))
rownames(regsim) = sort(unique(all$structure_name))
colnames(regsim) = unique(all$DISORDER)
for (d in unique(all$DISORDER)){
  regsim[all$structure_name[which(all$DISORDER == d)],d] = regsim[all$structure_name[which(all$DISORDER == d)],d] + 1
}  

png(file="DisorderRegionMatrixpt1.png",width=14,height=3,units="in",res=300)
save(regsim,file="DisorderRegionMatrixpt1.Rda")
load("DisorderRegionMatrixpt1.Rda") # binary values
load("regionCountMatrixpt1.Rda")  # counts
load("regionNormMatrixpt1.Rda") # normalized values
disty = dist(t(regsim))
hc = hclust(disty)
plot(hc,main="Disorder Region Similarity")
dev.off()

#pdf('/home/vsochat/sampleheatmap.pdf',width=60,height=4)
# STOPPED HERE - FIGURE OUT HOW TO GET MEANINGFUL FEATURES FOR GROUPS
hm = heatmap.2(t(regsim),dendrogram=c("none"),key=FALSE,col="heat.colors",labCol=NULL,main="Disorders Clustered by BrainTerm Regions")
hmresult = list(matrix=regsim,hm=hm,all=all)
save(hmresult,file="DisorderRegionHeatmappt1.Rda")
#dev.off()

# Calculte Tanimoto scores between disorders based on region overlap
tanimotos = array(0,dim=c(8,8))
disorders = sort(unique(all$DISORDER))
rownames(tanimotos) = sort(unique(all$DISORDER))
colnames(tanimotos) = sort(unique(all$DISORDER))
for (i in 1:length(disorders)){
  for (j in 1:length(disorders)){
    region1 = unique(as.character(all$structure_name[which(all$DISORDER == disorders[i])]))
    region2 = unique(as.character(all$structure_name[which(all$DISORDER == disorders[j])]))
    tanimotos[i,j] = length(intersect(region1,region2)) / length(union(region1,region2))
  }
}
# Save tanimoto distance matrix
png(file="DisorderRegionTanimotoMatrixpt1.png",width=10,height=8,units="in",res=300)
# Create matrix of these values for plot!
save(tanimotos,file="DisorderRegionTanimotoMatrixpt1.Rda")
disty = dist(as.dist(tanimotos))
hc = hclust(disty)
plot(hc,main="Disorder Region Similarity")
dev.off()


# First define groups based on region names
# Define groups based on similar sample locations
groups=list()
groups[[1]] = c("BRAINT")
groups[[2]] = c("ALZ","MS")
groups[[3]] = c("POMPE","PTSD")
groups[[4]] = c("FRX","RART")
groups[[5]] = c("SZO","LUPUS","AIDS","CYTOVIRUS")
groups[[6]] = c("RETT","ASD")
groups[[7]] = c("MD","BRA","CRANIO")
groups[[8]] = c("ALC","PARK")

gid = c()
gid["CYTOVIRUS"] = 5
gid["AIDS"] = 5
gid["LUPUS"] =5
gid["SZO"] = 5
gid["MD"] = 7
gid["BRA"] = 7
gid["CRANIO"] = 7
gid["ALZ"] = 2
gid["MS"] = 2
gid["FRX"] = 4
gid["RART"] = 4
gid["BRAINT"] = 1
gid["POMPE"] = 3
gid["PTSD"] = 3
gid["ASD"] = 6
gid["PARK"] = 8
gid["ALC"] = 8
gid["RETT"] = 6

# Append to matrix
tmp=c()
for (d in 1:length(all$DISORDER)){
  tmp = c(tmp,gid[all$DISORDER[d]])
}

all = cbind(all,tmp)
colnames(all)[23] = "REGIONGROUP"


# Define groups based on similar sample locations
groups=list()
groups[[1]] = c("CYTOVIRUS","AIDS","LUPUS","SZO")
groups[[2]] = c("MD","BRA","CRANIO")
groups[[3]] = c("ALZ","MS")
groups[[4]] = c("FRX","RART")
groups[[5]] = c("BRAINT")
groups[[6]] = c("POMPE")
groups[[7]] = c("PTSD")
groups[[8]] = c("ASD")
groups[[9]] = c("ALC")
groups[[10]] = c("RETT")
groups[[11]] = c("PARK")

# Let's append a group variable to the matrix
gid = c()
gid["CYTOVIRUS"] = 1
gid["AIDS"] = 1
gid["LUPUS"] = 1
gid["SZO"] = 1
gid["MD"] = 2
gid["BRA"] = 2
gid["CRANIO"] = 2
gid["ALZ"] = 3
gid["MS"] = 3
gid["FRX"] = 4
gid["RART"] = 4
gid["BRAINT"] = 5
gid["POMPE"] = 6
gid["PTSD"] = 7
gid["ASD"] = 8
gid["PARK"] = 9
gid["ALC"] = 10
gid["RETT"] = 11

tmp=c()
for (d in 1:length(all$DISORDER)){
  tmp = c(tmp,gid[all$DISORDER[d]])
}
gid[all$DISORDER]

all = cbind(all,tmp)
colnames(all)[22] = "GROUP"

save(all,file="SampleLocationsFDRpt1.Rda")

library(gplots)
heatmap.2(sampsim)

# For each group
overlappingRegions = c()
for (f in 1:length(groups)){
  members = groups[[f]]
  regions = all$structure_name[which(all$DISORDER==members[1])]
  if (length(members) > 1) {
    for (m in 2:length(members)){
      regions = intersect(regions,all$structure_name[which(all$DISORDER==members[m])])
    }
  }
  overlappingRegions = c(overlappingRegions,list(as.character(sort(unique(regions)))))
}  

# Calculate unique regions for each, as well as pairwise overlapping regions.
for (d in 1:ncol(regsim)){
  regions = names(which(regsim[,d] != 0))
  cat(regions,sep="\n",file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/overlap-regions/",colnames(regsim)[d],".txt",sep=""))  
}

# Now calculate pairwise overlap
for (f in 1:ncol(regsim)){
  for (g in 1:ncol(regsim)){
    d1 = colnames(regsim)[f]
    d2 = colnames(regsim)[g]
    overlap = intersect(names(which(regsim[,f] != 0)),names(which(regsim[,g] != 0)))
    name = paste(sort(c(d1,d2)),collapse="-")
    cat(overlap,sep="\n",file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/overlap-regions/",name,".txt",sep=""))    
  }  
}

# Save unique region names to file
reg = colnames(matrix)
save(reg,file="AllenUniqueRegions.Rda")
cat(reg,file="AllenUniqueRegions.txt",sep="\n")

# Save to file
for (r in 1:length(overlappingRegions)){
  cat("\n",groups[[r]],"\n",file="/home/vsochat/GroupedOverlappingRegionNameGroups.txt",append=TRUE)
  cat(overlappingRegions[[r]],file="/home/vsochat/GroupedOverlappingRegionNameGroups.txt",append=TRUE)
}

# Now plot data based on group
pdf("/home/vsochat/GroupSampleLocations.pdf")
for (f in 1:length(groups)){
  members = groups[[f]]
  subset = all[which(all$DISORDER==members[1]),]
  label = paste(members,collapse="-")
  s3d = scatterplot3d(all$mni_x,all$mni_y,all$mni_z,pch=20,color=sample(1:10,1),main=paste(label," Group Sample Locations Thresh .1"),xlab="MNIX",ylab="MNIY",zlab="MNIZ")
}
dev.off()

cat(intersect(all$structure_name[which(all$DISORDER=="RART")],all$structure_name[which(all$DISORDER=="FRX")]),sep="\n",file="/home/vsochat/overlapRART-FRX.txt")
allsub = subset(all,DISORDER %in% c("RART","FRX"))

# How do we select a subst of brain spatial locations for each disorder?
load("regionNormMatrixpt1.Rda")
load("DisorderRegionMatrixpt1.Rda")
#need a similarity metric that takes into account:
# variance of features in group
# scores of group members (higher == more important feature)
positives = as.vector(mnorm)
negatives = positives * -1
hist(c(positives,negatives))
normdist = c(positives,negatives)
std = sd(normdist)  

# Three sd is:
# FIRST calculate what would the range need to be to be in top 3 standard deviations of scores for ALL features? (eg, look at distribution of all scores, calculate Z scores, threshold Z score at value of 3, and then look at range of those values above threshold).  This is defined as $RANGE

# Then, to determine subsets of regions for each group of disorders:
# calculate average difference between pairwise scores
# Define the groups of disorders:
# Let's append a group variable to the matrix
gid = c()
gid["CYTOVIRUS"] = 3
gid["AIDS"] = 3
gid["LUPUS"] = 3
gid["SZO"] = 1
gid["MD"] = 3
gid["BRA"] = 4
gid["CRANIO"] = 4
gid["ALZ"] = 4
gid["MS"] = 4
gid["FRX"] = 3
gid["RART"] = 4
gid["BRAINT"] = 4
gid["POMPE"] = 3
gid["PTSD"] = 1
gid["ASD"] = 2
gid["PARK"] = 1
gid["ALC"] = 2
gid["RETT"] = 2

for (g in c(1,2,3,4)){
  # Select the members of the group
  subset = mnorm[names(which(gid==g)),]
  cat(names(which(gid==g)),"\n")
  # Get rid of zero features
  subset = subset[,-which(colSums(subset)==0)]
  # Calculate variance of column features (regions)
  vars = sort(apply(subset, 2, var))
  # Just look at nonzero
  subset[subset!=0] =1
  # Which regions are shared by all maps?
  cat(names(which(colSums(subset) == nrow(subset))),sep="\n")
  # Sort by variance
  subset = subset[,names(vars)]
  # Plot each
  dis = paste(names(which(gid==g)),collapse=" ")
  colors = c(1,2,3,4,5,6)
  plot(subset[1,],col="blue",pch=19,main=paste("Normalized Scores for",dis,"Across 100 Non-Zero Regions",ylab="Normalized Region Score"))
  for (s in 2:nrow(subset)){
    points(subset[s,],col=colors[s],pch=19)
  }
}
  
# Forget about regions - just try genes
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/DisorderGeneMatrixNESpt1.Rda")
disty = dist(genematrix)
hc = hclust(disty,method="ward")

# What about heatmap?
hm = heatmap.2(genematrix,dendrogram=c("none"),key=FALSE,col="heat.colors",labCol=NULL,main="Neuro2Gene Term Enrichment")

# Calculate just overlap for the neuro2gene interface
# Calculte Tanimoto scores between disorders based on region overlap
overlaps = array(0,dim=c(8,8))
disorders = sort(rownames(genematrix))
rownames(overlaps) = sort(rownames(genematrix))
colnames(overlaps) = sort(rownames(genematrix))
for (i in 1:length(disorders)){
  for (j in 1:length(disorders)){
    set1 = names(genematrix[disorders[i],which(genematrix[disorders[i],] != 0)]) 
    set2 = names(genematrix[disorders[j],which(genematrix[disorders[j],] != 0)])
    overlaps[i,j] = length(intersect(set1,set2))
  }
}
# THERE ARE NO OVERLAP!
save(overlaps,file="OverlapGeneCountN8pt1.Rda")

# VANESSA - HERE IS THE FINAL DECISION - TANIMOTO WITH GENES, and DIST MATRIX
# What about tanimoto based on genes?
# Calculte Tanimoto scores between disorders based on region overlap
tanimotogenes = array(0,dim=c(8,8))
disorders = sort(rownames(genematrix))
rownames(tanimotogenes) = sort(rownames(genematrix))
colnames(tanimotogenes) = sort(rownames(genematrix))
for (i in 1:length(disorders)){
  for (j in 1:length(disorders)){
    set1 = names(genematrix[disorders[i],which(genematrix[disorders[i],] != 0)]) 
    set2 = names(genematrix[disorders[j],which(genematrix[disorders[j],] != 0)])
    tanimotogenes[i,j] = length(intersect(set1,set2)) / length(union(set1,set2))
  }
}
save(tanimotogenes,file="TanimotoGenesN8pt1.Rda")

# Make tree with tanimoto as distance 
png(file="tanimotoTreeN8Wardpt1.png",width=9,height=8,units="in",res=300)
hc = hclust(as.dist(tanimotogenes),method="ward")
plot(hc,main="Tanimoto Matrix, Complete Criterion, for GENES")

# Heatmap
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

library("ggplot2")
setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults")
png(file="disorderGeneTanimotospt1.png",width=9,height=8,units="in",res=300)
p1 <- ggplot(mdf, aes(x=disorder, y=variable)) + 
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low="red",mid="white",high="blue") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 9 *0.8, angle = 90, hjust = 1, colour = "grey50"),axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p1
dev.off()