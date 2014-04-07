#!/home/vanessa/venv-python/bin/python

import sys
import neuro2gene

SEARCH_TERM = sys.argv[1]
THRESHOLD = sys.argv[2] # 0.001
analysisTime = sys.argv[3] #'2-12-2014'
dataset = neuro2gene.neurosynthInit()
features = neuro2gene.getFeatures(dataset)
sids,xyz,colnames,aba = neuro2gene.loadABA()

print "Starting analysis for term " + SEARCH_TERM
# First query neurosynth for relevant activation
#dataset = neuro2gene.neurosynthInit()
>>>>>>> aa25dac2d3dcd68c9015fc7aaea7c95b3f425ec6
# Return MNI coordinates where there is FDR corrected, significant activation
coords = neuro2gene.neurosynthQuery(SEARCH_TERM,THRESHOLD,dataset)
# Read in data file to get specimen IDs and structures to query
keepers,colnames = neuro2gene.ReadAllenBrainSample(coords)
# Get the analysis time for all output files   
# Print table to file of matches
neuro2gene.printSampleMatches(keepers,SEARCH_TERM,analysisTime,colnames)
