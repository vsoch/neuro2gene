#!/usr/bin/python

import neuro2gene
import datetime

# HERE IS RUNNING FOR ONE TERM

SEARCH_TERM = "associative"
THRESHOLD = 0.001
print "Starting analysis for term " + SEARCH_TERM
# First query neurosynth for relevant activation
dataset = neuro2gene.neurosynthInit()
# Return MNI coordinates where there is FDR corrected, significant activation
coords = neuro2gene.neurosynthQuery(SEARCH_TERM,THRESHOLD,dataset)
# Read in data file to get specimen IDs and structures to query
keepers,colnames = neuro2gene.ReadAllenBrainSample(coords)
# Get the analysis time for all output files
#analysisTime = datetime.datetime.now().strftime('%b-%d-%Y-%I-%M-%S')
analysisTime = '_4-4-2014'
# Print table to file of matches
neuro2gene.printSampleMatches(keepers,SEARCH_TERM,analysisTime,colnames)
# Now query to download expression data for our list
# neuro2gene.AllenQuery(keepers,SEARCH_TERM)        


# HERE IS RUNNING FOR ALL TERMS
import neuro2gene
dataset = neuro2gene.neurosynthInit('525')
features = neuro2gene.getFeatures(dataset)
THRESHOLD = 0.05
sids,xyz,colnames,aba = neuro2gene.loadABA()

#len(features)

for i in range(0,100):
  f = features[i]
  print "Saving sample subset for " + f + "..."
  coords = neuro2gene.neurosynthQuery(f,THRESHOLD,dataset)
  # Read in data file to get specimen IDs and structures to query
  keepers,colnames = neuro2gene.ReadAllenBrainSample(coords,sids,xyz,colnames,aba)
  neuro2gene.printSampleMatches(keepers,f,'_4-5-2014',colnames)
