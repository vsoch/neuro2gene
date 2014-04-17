# runShapleyProbeSets will submit jobs to run on cluster, for each of our behavioral terms

indir = "/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/thresh/9mmsq3000"
files = list.files(path = indir, pattern = "*.tab")
#uids = read.csv('/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/labels/featureUID.txt',sep=',',head=FALSE)
suids = read.csv('/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/labels/sampleUID.txt',sep=',',head=FALSE)
thresh = 0.05
b = 1000 # Number bootstrap samples
outdir = "/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/probeSets/9mmsq3000/"

for (f in 1:length(files)){
  file = paste(indir,"/",files[f],sep="")
  # List of probanda1es and pvalues
  outlist = paste(outdir,"list/",gsub(".tab","_probeSet",files[f]),sep="")
  # Expression data for the subset
  outfile = paste(outdir,"data/",gsub(".tab","_probeSet",files[f]),sep="")
  # We need to print a job file for each term to run:
  if (!file.exists(paste(outfile,"_MEdata_up.Rda",sep=""))) {
    jobby = paste(files[f],".job",sep="")
    sink(paste(".jobs/",jobby,sep=""))
    cat("#!/bin/bash\n")
    cat("#SBATCH --job-name=",jobby,"\n",sep="")  
    cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
    cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
    cat("#SBATCH --time=2-00:00\n",sep="")
    cat("#SBATCH --mem=12000\n",sep="")
    cat("Rscript /home/vsochat/SCRIPT/R/neurosynth/shapleyProbeSets.R",file,thresh,outfile,outlist,b,"\n")
    sink()
  
    # SUBMIT R SCRIPT TO RUN ON CLUSTER  
    system(paste("sbatch",paste(".jobs/",jobby,sep="")))
  }
}
