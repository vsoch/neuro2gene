# This script queries the allen brain atlas to get a list of
# structure ids, and prings to file
#
# vsochat@stanford.edu

import api
import numpy as np
import json
import nibabel as nb
    
def ReadAllenBrainStructures():

    # This file has a list of sample IDs and MNI coordinates
    filey = open('data/AllenBrainAtlasSampleMNISorted.csv','r')
    strucID = []
    temp = filey.readlines()

    for f in range(1,len(temp)):
      strucID.append(int(temp[f].split(',')[2]))

    strucID = np.unique(strucID)

    structures = api.getStructures(strucID)
    
    # Return specimen
    return structures    
    

def printStructures(structures):
    filey = open('data/AllenBrainAtlasStructures.csv','w')
    filey.writelines("id\tname\tweight\tacronym\tparent_structure_id\tsphinx_id\themisphere_id\tcolor_hex_triplet\tstructure_id_path\tdepth\tneuro_name_structure_id_path\tneuro_name_structure_id\tontology_id\tatlas_id\n")
    for s in structures:
      filey.write(s[0] + "\n")
    filey.close()


if __name__ == "__main__":

    # Read in data file to get specimen IDs and structures to query
    structures = ReadAllenBrainStructures()
    printStructures(structures)
