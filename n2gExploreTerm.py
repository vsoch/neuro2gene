#!/bin/python

# This script will explore the derivation of a particular feature
# It uses neuro2gene module
# We would want to be able to take a term map, and understand what the articles "were getting at" from which the term is derived.  We will use the neurosynth API to retrieve information about the articles, download them from Pubmed, read the full text, and create a document vector for each.  We will then cluster similar articles, and then visualize the similar articles with a word cloud (wordle!)

# vsochat 4/12/2014

import neuro2gene
import os
THRESHOLD = 0.001

# Read in features from dataset
dataset = neuro2gene.neurosynthInit()
queryfeatures = ['angry','noun','skills','covert']
features = neuro2gene.getFeatures(dataset)
email = 'vsochat@stanford.edu'

# Make sure user selected features are in our list
if len([e for e in queryfeatures if e in '\n'.join(features)]) != len(queryfeatures):
  print "Error: Please check that query features are in feature list!"

word_dicts = list()

# For each feature, get the articles
for f in queryfeatures:
  word_dicts.append(neuro2gene.getWordCounts(dataset,f,THRESHOLD,email))

# Print to file
for wrd in range(0,len(word_dicts)):
  w = word_dicts[wrd]
  keys = w.keys()
  vals = w.values()
  filey = open("/home/vanessa/Desktop/" + queryfeatures[wrd] + "_wordcounts.txt",'w')
  filey.writelines('word\tcount\n')
  for i in range(0,len(w)):
    filey.writelines(keys[i] + '\t' + str(w[keys[i]]) + '\n')
  filey.close()


