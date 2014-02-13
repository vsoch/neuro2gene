# This script queries neurosynth for a term of interest, does
# a meta analysis to get a spatial map of voxels from the
# literature associated with that term, and then queries
# the allen brain atlas to extract probes of interest, and
# write their pa call (if they are significantly above background)
# to a table in big query.  This table can be downloaded to
# find a subset of genes that have expression above background 
# (I will likely use R for this analysis)
#
# vsochat@stanford.edu
# Wall-Lab, February 2014

import numpy as np
import json
from neurosynth.base.dataset import Dataset
from neurosynth.base import imageutils
from neurosynth.analysis import meta
import nibabel as nb
from nibabel import nifti1
from copy import deepcopy
import subprocess
import os
import matplotlib.pyplot as plt
import pickle

# -- NEUROSYNTH FUNCTIONS --------------------------------------------------------------

def neurosynthInit():
    print "Initializing Neurosynth database..."
    dataset = Dataset('data/database.txt')
    dataset.add_features('data/features.txt')    

    #print "Loading standard space brain..."
    #img = nb.load("data/MNI152_T1_2mm_brain.nii.gz")
    #standard = img.get_data()
    return dataset

def neurosynthQuery(searchTerm,thresh,dataset,outdir=None):
    query = dataset.get_ids_by_features('*' + searchTerm + '*',threshold=thresh)
    ma = meta.MetaAnalysis(dataset,query)
    # This gets the absolute value FDR corrected at threshold
    data = ma.images[ma.images.keys()[4]]
    # Print this image to file, to look at later
    if not outdir:
      imageutils.save_img(data, '%s.nii.gz' % (searchTerm), dataset.volume)
    else:
      imageutils.save_img(data, '%s/%s.nii.gz' % (outdir, searchTerm), dataset.volume)
    img = dataset.volume.unmask(data)
    # These are x,y,z coordinates of nonnzero voxels    
    idx = np.nonzero(img)
    affine =  dataset.volume.volume.get_affine()
    coords = np.dot(affine,[ idx[0],idx[1],idx[2] ,np.ones(len(idx[0]))])
    coords = np.transpose(coords)
    # Get rid of 4th column
    coords = np.delete(coords,-1,1)
    # These are MNI coordinates for non-zero voxels
    return coords

def getFeatures(dataset):
    return dataset.get_feature_names()


# -- SIMILARITY METRICS ----------------------------------------------------------------

# Take a list of features, and calculate similarity of maps
# Similarity is based on shared voxels - non-shared voxels.
# difference_score = average_activation_per_voxel_shared - average_activation_per_voxel_notshared
def getSimilarity(features,thresh,dataset,outdir=None):

    # Here we hold our distance scores, smaller = more similar
    sims = nu.zeros(shape=(len(features),len(features)))

    for i in range(1,len(features)):
      print "Starting feature #: " + str(i) + "/" + str(len(features)) + " " + f1
      f1 = features[i]
      query = dataset.get_ids_by_features('*' + f1 + '*',threshold=thresh)
      ma = meta.MetaAnalysis(dataset,query)
      # Absolute value for FDR corrected at threshold
      data = ma.images[ma.images.keys()[4]]
      for j in range(1,len(features)):
        f2 = features[j]
        query = dataset.get_ids_by_features('*' + f2 + '*',threshold=thresh)
        ma = meta.MetaAnalysis(dataset,query)
        data2 = ma.images[ma.images.keys()[4]]
        # Count number of shared voxels
        overlap = float(np.sum((data != 0) & (data2 != 0)))
        # Count number of nonzero voxels in both
        total = float(np.sum(data !=0) + np.sum(data2 !=0))
        # As a percentage of total number of voxels for both maps
        sims[i,j] = np.divide(overlap,(total - overlap))
     
    # Print similarity matrix to file
    np.savetxt('featureSimMatrix.dat', sims, fmt='%.18e', delimiter=',', newline='\n', header=",".join(features))
    # How to read in data file
    #np.loadtxt('featureSimMatrix.dat', comments='#', delimiter=',')
    return sims
 

# -- ALLEN BRAIN ATLAS FUNCTIONS -----------------------------------------------------------

# Output text file and pickle dictionary with unique ID for each sample
def AllenBrainSampleUID():
    # Load list of sample IDs and MNI coordinates, 3702 across 6 specimens
    filey = open('data/AllenBrainSamplesSortedMNITab.csv','r')
    aba = []
    for f in filey:
      aba.append(f.strip('\n').split('\t'))
    colnames = aba.pop(0)
    colnames = list(["uid"]) + colnames
    uids = dict()
    print "Assigning unique id to " + str(len(aba)) + " samples in Allen Brain Atlas..."
    filey = open('sampleUID.txt','w')
    filey.writelines(colnames + list(["\n"]))
    for a in range(0,len(aba)):
      f = aba[a]
      if a < 10:
        filey.writelines("SAMP000 " + str(a) + "," + "-".join(f[0:6] + f[7:14]) + "\n")
        uids["-".join(f[0:6] + f[7:14])] = "SAMP000 " + str(a)
      elif a < 100:
        filey.writelines("SAMP00 " + str(a) + "," + "-".join(f[0:6] + f[7:14]) + "\n")
        uids["-".join(f[0:6] + f[7:14])] = "SAMP00 " + str(a)
      elif a < 1000:
        filey.writelines("SAMP0 " + str(a) + "," + "-".join(f[0:6] + f[7:14]) + "\n")
        uids["-".join(f[0:6] + f[7:14])] = "SAMP0 " + str(a)
      else:
        filey.writelines("SAMP " + str(a) + "," + "-".join(f[0:6] + f[7:14]) + "\n")
        uids["-".join(f[0:6] + f[7:14])] = "SAMP " + str(a)
    filey.close()

    # Save dictionary to pickle file as well
    with open('sampleUID.pkl', 'wb') as f:
      pickle.dump(uids, f, pickle.HIGHEST_PROTOCOL)

# Return 
def loadABA():
    # Load list of sample IDs and MNI coordinates, 3702 across 6 specimens
    filey = open('data/AllenBrainSamplesSortedMNITab.csv','r')
    aba = []
    for f in filey:
      aba.append(f.strip('\n').split('\t'))
    colnames = aba.pop(0)
    xyz = []   # To find all xyz coordinates
    sids = []  # to find list of unique specimens
    for f in aba:
      xyz.append([float(f[11]),float(f[12]),float(f[13])])
      sids.append(f[0])

    # Find the number of unique specimens
    sids = np.unique(sids)
    return sids,xyz,colnames,aba
    

def ReadAllenBrainSample(coords,sids=None,xyz=None,colnames=None,aba=None):

    if not aba:
      sids,xyz,colnames,aba = loadABA()

    print "Finding close sample locations in Allen Brain Atlas..." 
    # This will be a list of structures, specimens, etc to keep
    keepers = []
    keepers_string = []  # for checking if we have an entry

    # The column names will have original header, plus point
    # MNI coordinate, plus the distance
    colnames = colnames + ['nsyn_mni_x','nsyn_mni_y','nsyn_mni_z','sq_distance']

    # First, find the closest len(sids) points for each coordinate
    # the idea being that we might have one close sample / specimen
    for c in coords:
      aba_temp = deepcopy(aba)
      xyz_temp = deepcopy(xyz)  # Hold a temporary list of samples
      for i in range(0,1): #len(sids)):   # Find closest point for each
        idx = np.argmin(np.sum((xyz_temp - c)**2, axis=1))
        point = aba[idx] + [ str(c[0]),str(c[1]),str(c[2]),str(np.min(np.sum((xyz_temp - c)**2, axis=1)))]
        if "".join(aba[idx]) not in keepers_string: #rownum,structure,parent_structure,mnix,mniy,mniz,nsynx,nsyny,nsynz,distance
          keepers.append(point)
        keepers_string.append("".join(aba[idx]))
        del aba_temp[idx]
        del xyz_temp[idx]

      # TODO: if we want to threshold based on a distance, add a "contenders" variable
      # In the loop, and add subset to keepers
    print "Found " + str(len(keepers)) + " close points."

    # Return list of hitsh
    return keepers,colnames


def printSampleMatches(keepers,searchTerm,analysisTime,colnames):
    print "Writing Allen Brain Atlas Matches for " + searchTerm + " to file..."
    filey = open('output/' + searchTerm + "_" + analysisTime + '_samples.tab','w')
    filey.writelines("\t".join(colnames) + "\n")
    for k in keepers:
      filey.writelines("\t".join(k) + "\n")
    filey.close()



# THIS FUNCTION IS NOT COMPLETE - BIG QUERYING IS TOO SLOW!
def AllenQuery(keepers,searchTerm):
    print "Running queries for " + str(len(keepers)) + " samples..."

    # First we need to organize the data but the subject ID - we will query for each
    # subject separately
    # First we need to find unique IDs for our subjects, structures, and coordinates
    sids = np.unique([x[0] for x in keepers])
    data = dict()
    
    for s in sids:
      data[str(s)] = list()

    for k in keepers:
      data[k[0]].append(k)

 
    for key,val in data.iteritems():
      subid = key
      struc = np.unique([x[1] for x in keepers])
      #mnix = np.unique(["\'" + x[11] + "\'" for x in val])
      #mniy = np.unique(["\'" + x[12] + "\'" for x in val])
      #mniz = np.unique(["\'" + x[13] + "\'" for x in val])
      mnix = np.unique([x[11] for x in val])
      mniy = np.unique([x[12] for x in val])
      mniz = np.unique([x[13] for x in val])
      #nx = np.unique(["\'" + x[8] + "\'" for x in val])
      #ny = np.unique(["\'" + x[9] + "\'" for x in val])
      #nz = np.unique(["\'" + x[10] + "\'" for x in val])

      # Put in comma separated list
      struc = ",".join(struc)
      mnix = ",".join(mnix)
      mniy = ",".join(mniy)
      mniz = ",".join(mniz)

      for t in range(0,12):
        bqCommand = "bq query --destination_table=allen_brain_result." + searchTerm +"_" + str(t) + " --append_table \"SELECT * FROM allen_brain_human.pacall" + str(t) + " WHERE (specimen_id in (" + subid + ")) AND (structure_id in (" + struc + ")) AND (mni_x in (" + mnix + ")) AND (mni_y in (" + mniy + ")) AND (mni_z in (" + mniz + "))\""
        filey = open(searchTerm + "-query-" + str(subid) + "_" + str(t),'w')
        filey.writelines('#!/bin/sh\n')
        filey.writelines(bqCommand)
        filey.close()
        os.system('chmod u+x ' + searchTerm + "-query-" + str(subid) + "_" + str(t))
        #os.system(bqCommand)

    print "Querying is complete! Results can be found in tables allen_brain_result." + searchTerm + "0..11 in Google Big Query."
          
if __name__ == "__main__":
  print "Please import as a module"
