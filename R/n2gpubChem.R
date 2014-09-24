# This script will look up disorder genes in pubChem, specifically finding
# bioassays that are "active" for the gene - the idea is that we
# want to look at "druggable genes" - enriched genes in a disorder 
# that can be impacted by drugs!

# Set working directory to where we have data
library("RCurl")
setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults")
load("disorderGenesN8pt1.Rda")
load("disorderGeneListsN8pt1.Rda")

# For each disorder, create a list of genes with active
# compounds
disorders = as.character(unique(genes$DISORDER))

# First get list of gene IDs
library("org.Hs.eg.db")
# Convert the object to a list
lookup = as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
lookup = lookup[!is.na(lookup)]

disorderaids = list()
for (i in 1:length(disorders)){
  dis = disorders[i]
  cat("Processing",dis,"\n")
  gen = as.character(unique(genes$PROBE[which(genes$DISORDER == dis)]))
  # This list will hold aids for each gene
  aidlist = list()
  for (g in 1:length(gen)){
    cat("Processing",g,"of",length(gen),"\n")
    # get number
    querygene = as.character(gen[g])    
    geneid = lookup[[querygene]]
    # Look up assay ids
    query = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol/",querygene,"/aids/TXT",sep="")
    result = getURL(query)
    aids = strsplit(result,"\n")[[1]]
    for (k in geneid){
      aidlist[[k]] = aids
    }
  }
  disorderaids[[dis]] = aidlist
}

# Save object with AIDS
save(disorderaids,file="/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult/disorderGeneAssayLists8.Rda")

# Now we should look up each of the SIDS, and find the SIDS for which the gene is active
bioassayFolder = "/scratch/PI/dpwall/DATA/PUBCHEM/bioassay"
subfolders = list.files(bioassayFolder,pattern="^[0-9]{7}_[0-9]{7}$")
# Convert to ranges
starts = c()
for (s in 1:length(subfolders)){
   folder = subfolders[s]
   start = as.numeric(strsplit(folder,"_")[[1]][1])
   end = as.numeric(strsplit(folder,"_")[[1]][2])
   starts = c(starts,start)
}
names(starts) = subfolders

# Save to file for later!
save(starts,file="pubChemfolderStartIndex.Rda")

# Now for each of our disorderaids, read in the files, save active sid info
disordersids = list()
for (b in 1:length(disorderaids)){
   sids = c()
   dis = names(disorderaids[b])
   cat("Starting",dis,":",b,"of",length(disorderaids),"\n")
   genetocheck = disorderaids[[b]]
   for (a in 1:length(genetocheck)){
     geneid = names(genetocheck)[a]
     aids = genetocheck[[a]]
     for (c in 1:length(aids)){
     # Find the file to lookup the aid
       aid = as.numeric(aids[c])
       # Get the index for the right folder
       idx = as.numeric(which(starts<=aid)[length(which(starts<=aid))])
       folder = names(starts)[idx]
       # Read in file
       filey = paste(bioassayFolder,"/",folder,"/",aid,".csv.gz",sep="")
       filey = read.csv(filey,head=TRUE,sep=",")
       idx = which(colnames(filey)=="X3.Gene.Target.ID.STRING....gene.target.id")
       idx = which(as.character(filey[,idx])==geneid)
       sids = rbind(sids,filey[idx,])
     }
   }
   disordersids[[dis]] = sids
}
