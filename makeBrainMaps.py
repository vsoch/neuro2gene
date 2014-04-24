import neuro2gene
thresh = 0.001
import os.path

# Read in features from dataset
dataset = neuro2gene.neurosynthInit('525')
features = neuro2gene.getFeatures(dataset)

count = 1
for f in features:
   fname = "/scratch/users/vsochat/DATA/BRAINMAP/nsynth525/" + f + ".nii.gz"
   if not os.path.isfile(fname):
     print "Processing " + str(count) + " of " + str(len(features)) + " : " + f
     tmp = neuro2gene.neurosynthQuery(f,thresh,dataset,outdir="/scratch/users/vsochat/DATA/BRAINMAP/nsynth525/")
     count = count + 1

