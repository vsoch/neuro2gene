#!/usr/bin/python

# This script will take in a list of gene IDs, lookup probe IDs that we have 
# in the Allen Brain Atlas, and return expression data for all samples for each
# probe

# This top portion is for route 66 genes

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


# This bottom portion is for generating a matrix of genes from significant GSEA results - we will calculate
# mean expression values for different brain maps of the significant results, but we need the CORE genes from all sets

# These genes are from angry, color, audiovisual, noun, rejection brain terms sets
GENES =["A4GNT","ACAP1","ACSL3","ADH6","ADM2","AGXT2","ANKK1","ANXA10","ARG1","ARHGAP11A","ARL13B","ART4","ASB4","ASF1B","ATP2A1","BAALC","BBS5","BEX2","BHLHE41","C12orf76","C16orf3","C16orf57","C17orf53","C1orf116","C22orf29","C6orf201","C6orf89","C8orf45","CA5B","CACNA1I","CALM1","CALM3","CD4","CD86","CD99L2","CDH29","CDKAL1","CEBPD","CFTR","CHMP4B","CIITA","CLLU1","CLTC","CPE","CST1","CXorf21","DNAJC8","EEF1D","EFHA2","EFNB1","ENAM","ESF1","EXOC3L2","FAHD2B","FAM166A","FAM174A","FAM181B","FAM24A","FBRSL1","FBXL16","FER1L4","FHIT","FLJ37396","FOXQ1","FUT6","GATM","GFRA4","GLRA4","GLRX2","GPIHBP1","GPX1","GPX4","GRB2","GRM2","HAL","HDGF","HFM1","HMGCS2","HOXA3","HOXC4","HSDL2","HSP90B1","ICT1","IDO2","IL23A","ILDR1","ITGB1BP3","IYD","JUND","KIF5A","KLC3","KRT81","KRTAP10-10","LAMB4","LAMP1","LCE1E","LDHB","LGALS1","LHFPL4","LHPP","LOC391722","LOC642031","LOC729866","LOC730256","LONRF2","LSMD1","MAPK4","MDP1","MEG3","MGEA5","MIER2","MKNK2","MOGAT2","MTPN","MYEOV2","MYO7A","NACA2","NDRG2","NDUFA11","NDUFS7","NEB","NEDD4","NF2","NIN","NME7","NODAL","NPAS3","NSF","ODF1","OLAH","OR10G3","OR1J1","OR2T8","OR52I1","OR52I2","OR5H1","OR5I1","OR6K3","ORC1L","OTOR","PALM","PARP10","PCNP","PCNXL2","PDC","PDDC1","PEBP1","PGM2L1","PLSCR2","POFUT1","PPP3CB","PRKCZ","PRR11","PTGER2","PTGER4","PXN","RNF208","RNF26","RPL22","RPL23A","RPL26","RPL35","RPL36","RPL37","RPL3L","RPP25","RPS12","RPS6KA2","RPS6KA3","RPSA","RSL1D1","SCGB1D4","SEMG1","SENP5","SERPINA4","SERTAD3","SHISA4","SI","SIK1","SIRPA","SLAMF1","SLC25A22","SLC29A3","SLC35E4","SLC35F2","SLC9A3R1","SMARCC2","SMYD4","SOAT2","SOHLH2","SPATA2L","SPRR2B","STMN1","STMN2","STRA13","STX11","STXBP2","SYT1","TAS2R38","TBX22","TBX5","TCHHL1","THOC7","TM4SF4","TMX4","TPI1","TPSD1","TRIM29","TRIM38","TRIM48","TSPAN10","TSPAN3","TSPYL2","TTC3","TTC33","TUBB","UBC","UBE2O","UGT2B17","UQCRQ","VIT","VNN2","VNN3","VPS4A","WDR87","WNT6","YBX2","ZC3H11A","ZCCHC16","ZFAND5","ZFPM1","ZNF24","ZNF493","ZNF579","ZNF746","ZNF92","ZYG11A","ZZEF1"]

# Lookup pids for genes we have expression for
print "Looking up " + str(len(GENES)) + " in Allen Brain Atlas..." 
[abaid,pids,genes] = neuro2gene.lookup(GENES)

# Print to file for import into R
filey = open('/home/vanessa/Desktop/brainterm6GenesExpression.txt','w')
filey.writelines('ABA_ID\tAFFY_ID\tGENE_ID\n')
for p in range(0,len(pids)):
  filey.writelines(abaid[p] + "\t" + pids[p] + "\t" + genes[p] + "\n")

