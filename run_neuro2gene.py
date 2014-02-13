#!/usr/bin/python

import neuro2gene
import datetime

# HERE IS RUNNING FOR ONE TERM

SEARCH_TERM = "anxiety"
THRESHOLD = 0.001
print "Starting analysis for term " + SEARCH_TERM
# First query neurosynth for relevant activation
dataset = neuro2gene.neurosynthInit()
# Return MNI coordinates where there is FDR corrected, significant activation
coords = neuro2gene.neurosynthQuery(SEARCH_TERM,THRESHOLD,dataset)
# Read in data file to get specimen IDs and structures to query
keepers,colnames = neuro2gene.ReadAllenBrainSample(coords)
# Get the analysis time for all output files
analysisTime = datetime.datetime.now().strftime('%b-%d-%Y-%I-%M-%S')
# Print table to file of matches
neuro2gene.printSampleMatches(keepers,SEARCH_TERM,analysisTime,colnames)
# Now query to download expression data for our list
# neuro2gene.AllenQuery(keepers,SEARCH_TERM)        


# HERE IS RUNNING FOR ALL TERMS
dataset = neuro2gene.neurosynthInit()
features = neuro2gene.getFeatures(dataset)
THRESHOLD = 0.001
sids,xyz,colnames,aba = neuro2gene.loadABA()

for f in features:
  print "Saving sample subset for " + f + "..."
  coords = neuro2gene.neurosynthQuery(f,THRESHOLD,dataset)
  # Read in data file to get specimen IDs and structures to query
  keepers,colnames = neuro2gene.ReadAllenBrainSample(coords,sids,xyz,colnames,aba)
  neuro2gene.printSampleMatches(keepers,f,'_2-12-2014',colnames)


