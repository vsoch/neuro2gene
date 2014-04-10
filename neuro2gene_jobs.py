# This script will load neurosynth features, and create job files to submit
# on the sherlock cluster
#
# vsochat 4/4/2014

import neuro2gene
import os

THRESHOLD = 0.001
analysisTime = '4-9-2014'

# Read in features from dataset
dataset = neuro2gene.neurosynthInit()
features = neuro2gene.getFeatures(dataset)
sids,xyz,colnames,aba = neuro2gene.loadABA()

# For each feature, print job file for running on cluster
for f in features:
  print "Creating job file for " + f + "..."
  jobfile = open(".jobs/" + f + ".job",'w')
  jobfile.writelines("#!/bin/bash\n")
  jobfile.writelines("#SBATCH --job-name=" + f + "_genes.job\n")
  jobfile.writelines("#SBATCH --output=.out/"+ f + "_genes.out\n")
  jobfile.writelines("#SBATCH --error=.out/"+ f + "_genes.err\n") 
  jobfile.writelines("#SBATCH --time=2-00:00\n") 
  jobfile.writelines("#SBATCH --mem=12000\n")
  jobfile.writelines("source /home/vsochat/python-lapack-blas/bin/activate\n")   
  jobfile.writelines("/home/vsochat/python-lapack-blas/bin/python /home/vsochat/SCRIPT/python/neuro2gene/run_neuro2gene_cluster.py " + f + " " + str(THRESHOLD) + " " + analysisTime + "\n")  
  jobfile.close()

# submit jobs
for i in range(2050,2050+512):
  f = features[i]
  os.system('sbatch .jobs/'+ f + ".job")



