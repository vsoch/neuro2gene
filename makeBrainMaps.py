#!/usr/bin/python

import neuro2gene

# First query neurosynth for relevant activation
dataset = neuro2gene.neurosynthInit()
# Get the feature names
features = neuro2gene.getFeatures(dataset)

# OPTION 1
# Just create similarity matrix without images
sims = neuro2gene.getSimilarity(features,THRESHOLD,dataset,outdir)

# OPTION 2 make images
# Set the output directory for images
outdir = '/home/vanessa/Documents/Work/ATLAS/NEUROSYNTH'

# For each feature, create spatial map
for f in features:
  SEARCH_TERM = f
  THRESHOLD = 0.001
  # Return MNI coordinates where there is FDR corrected, significant activation
  coords = neuro2gene.neurosynthQuery(SEARCH_TERM,THRESHOLD,dataset,outdir)
  

# Now we need to run analysis for all the terms!
