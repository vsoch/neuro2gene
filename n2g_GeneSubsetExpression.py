#!/usr/bin/python

# This script will take in a list of gene IDs, lookup probe IDs that we have 
# in the Allen Brain Atlas, and return expression data for all samples for each
# probe

import neuro2gene

# Define list of genes to lookup in ABA
# These are Leticia's Core enrichment genes for route 66
GENES = ["SULF2","ARID4B","SVIL","AGTPBP1","UPF2","SORL1","RBM25","PPM1B","ACSL4","RBM39","SIRPA","TM9SF3","CNOT4","ZNF644","NFE2L2","ITGB1","UBE2D3","SEPT2","OCIAD1","PCNP","TIMP2","COL4A3BP","WIPF1","CPD","RAB24","PHF20L1","WLS","SLC44A2","ASAH1","CUL4A","BAX","MAP1LC3B","LYST","VPS29","FAM198B","LAMP2","CSF2RA","DSE","NUDT16","CNBP","BCL6","PDE4B","HSPD1","TDP2","RHEB","CAPZA2","RNF141","PANK2","RAB2A"]

# Lookup pids for genes we have expression for
print "Looking up " + str(len(GENES)) + " in Allen Brain Atlas..." 
[abaid,pids,genes] = neuro2gene.lookup(GENES)

# Print to file for import into R
filey = open('/home/vanessa/Desktop/route66ABA.txt','w')
filey.writelines('ABA_ID\tAFFY_ID\tGENE_ID\n')
for p in range(0,len(pids)):
  filey.writelines(abaid[p] + "\t" + pids[p] + "\t" + genes[p] + "\n")


