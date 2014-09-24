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

for (i in 1:length(disorders)){
  dis = disorders[i]
  gen = as.character(genes$PROBE[which(genes$DISORDER == dis)])
  # This list will hold aids for each gene
  aidlist = list()
  for (g in 1:length(gen)){
    # get number
    querygene = as.character(gen[g])    
    geneid = lookup[[querygene]]
    # Look up assay ids
    query = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol/",querygene,"/aids/TXT",sep="")
    result = getURL(query)
    aids = strsplit(result,"\n")[[1]]
    aidlist[[geneid]] = aids
    # We will be looking up these results in our data on sherlock
  }
}

# Run this on Sherlock
4) in data file - find the gene, if active, save SID for it!
  query = paste("http://www.ncbi.nlm.nih.gov/gene/?term=",querygene,sep="")
# This will get assay ids
