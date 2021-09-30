#!/usr/bin/env Rscript
## functions
checkfile <- function(filename) {
if( file.access(filename) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", filename))
}
}

debug=0

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", 
                    action="store_true", 
                    default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-r", "--readspergene", 
                    type="character", 
                    help="STAR reads per gene file",
                    required=TRUE)
parser$add_argument("-s", "--strandedness", 
                    type="character",
                    choices=c('no', 'yes', 'reverse'),
                    help="yes/no/reverse default [\"%(default)s\"]",
                    default="reverse")
parser$add_argument("-g","--genelengthlookup", 
                    type="character",
                    help="Gene length lookup table with 3 columns (ensemblID,HUGOname,mRNAlength)",
                    required=TRUE)
parser$add_argument("-o","--outfile", 
                    type="character",
                    help="STAR TPM/RPKM/rawcounts/HUGOname file",
                    required=TRUE)


args <- parser$parse_args()

readspergene=args$readspergene
genelengthlookup=args$genelengthlookup
strandedness=args$strandedness
outfile=args$outfile

if (debug==1){
###DEBUGGING#########
readspergene="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis/results/596-7-2/STAR/withChimericJunctions/596-7-2_p1.ReadsPerGene.out.tab"
genelengthlookup="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/genelist.lookup.tsv"
strandedness="reverse"
outfile="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis/results/596-7-2/STAR/withChimericJunctions/596-7-2_p1.counts.tsv"
#####################
}

checkfile(readspergene)
checkfile(genelengthlookup)

# read in the rawcounts from STAR
d=read.csv(readspergene,sep="\t",header = FALSE)
if (strandedness=="no"){
  d=d[!grepl("^N_",d$V1),c(1,2)]
} else if (strandedness=="yes"){
  d=d[!grepl("^N_",d$V1),c(1,3)]
} else if (strandedness=="reverse"){
  d=d[!grepl("^N_",d$V1),c(1,4)]
}
colnames(d)=c("ensemblID","raw_counts")

# read in the lookuptable
l=read.csv(genelengthlookup,sep="\t",header=FALSE)
colnames(l)=c("ensemblID","gene_name","mRNA_length")
m=merge(l,d,by="ensemblID")
rate=m$raw_counts/m$mRNA_length
m$rpkm=rate / sum(m$raw_counts) * 1e6
m$tpm=rate / sum(rate) * 1e6

write.table(m,file=outfile,row.names = FALSE,quote=FALSE,sep="\t")
