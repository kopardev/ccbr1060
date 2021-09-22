#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("sets"))

## functions
checkfile <- function(filename) {
  if( file.access(filename) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", filename))
  }
}

readfile <- function(filename,indexcols,datacol,samplename) {
  d=read.csv(filename,header=TRUE,sep="\t")
  cols=colnames(d)
  reqdcols=c(indexcols,datacol)
  if(set_is_subset(as.set(reqdcols),as.set(cols))){
    d=d[,reqdcols]
    cols=colnames(d)
    cols=gsub(pattern = datacol,replacement = samplename,cols)
    colnames(d)=cols
    return(d)
  } else {
    stop(sprintf("Required columns: %s are missing in file %s", reqdcols,filename))
  }
}

debug=0



# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-l", "--filelist", 
                    type="character", 
                    help="comma separated list of files",
                    required=TRUE)
parser$add_argument("-i", "--indexcols", 
                    type="character",
                    help="comma separated list of columns to use as index while merging. These need to exist in all files provided in the -l option.",
                    required=TRUE)
parser$add_argument("-c","--datacol", 
                    type="character",
                    help="column to be extracted from all samples and aggregated in the outputfile. Should exist in all files provided to -l. Its name will be replaced by corresponding value in -s argument.",
                    required=TRUE)
parser$add_argument("-s","--samplenames", 
                    type="character",
                    help="comma separated list of sample names. Need to be unique. Will replace datacol (-c) when reported in the output file. ",
                    required=TRUE)
parser$add_argument("-o","--outfile", 
                    type="character",
                    help="aggregated outfile",
                    required=TRUE)


args <- parser$parse_args()

filelist=unlist(strsplit(args$filelist,","))
indexcols=unlist(strsplit(args$indexcols,","))
datacol=args$datacol
samplenames=unlist(strsplit(args$samplenames,","))
outfile=args$outfile

if(debug==1){
  filelist=unlist(strsplit("/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis/results/596-7-2/STAR/withChimericJunctions/a.tsv,/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis/results/596-7-2/STAR/withChimericJunctions/b.tsv",","))
  indexcols=unlist(strsplit("ensemblID,gene_name,mRNA_length",","))
  datacol="tpm"
  samplenames=unlist(strsplit("A,B",","))
  outfile="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis/results/596-7-2/STAR/withChimericJunctions/test.tpm.tsv"
}

# 2 or more files required
if (length(filelist)<=1) {
  stop(sprintf("Two are more comma separated files needed"))
}
# nfiles and nsamplenames should match
if (length(unique(filelist))!=length(unique(samplenames))){
  stop(sprintf("Number of files need to match number of samplenames. Duplicate files and samplenames not allowed!"))
}

for (f in filelist) {checkfile(f)}

m=readfile(filename = filelist[1],
           indexcols = indexcols,
           datacol = datacol,
           samplename = samplenames[1])

for (i in 2:length(filelist)){
  n=readfile(filename = filelist[i],
             indexcols = indexcols,
             datacol = datacol,
             samplename = samplenames[i])
  m=merge(m,n,by=indexcols)
}

write.table(m,file=outfile,quote = FALSE,sep="\t",row.names = FALSE)
