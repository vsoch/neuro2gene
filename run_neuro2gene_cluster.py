#!/usr/bin/python

import sys
import neuro2gene

SEARCH_TERM = sys.argv[1]
dataset = neuro2gene.neurosynthInit()
features = neuro2gene.getFeatures(dataset)
THRESHOLD = 0.001
sids,xyz,colnames,aba = neuro2gene.loadABA()
analysisTime = '2-12-2014'

print "Starting analysis for term " + SEARCH_TERM
# First query neurosynth for relevant activation
dataset = neuro2gene.neurosynthInit()
# Return MNI coordinates where there is FDR corrected, significant activation
coords = neuro2gene.neurosynthQuery(SEARCH_TERM,THRESHOLD,dataset)
# Read in data file to get specimen IDs and structures to query
keepers,colnames = neuro2gene.ReadAllenBrainSample(coords)
# Get the analysis time for all output files   
# Print table to file of matches
neuro2gene.printSampleMatches(keepers,f,analysisTime,colnames)
    

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


