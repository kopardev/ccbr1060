#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-r", "--rawcountsmatrix", 
                    type="character", 
                    help="file with raw counts matrix",
                    required=TRUE)
parser$add_argument("-c", "--coldata", 
                    type="character", 
                    help="two tab delimited columns.. sample_name and condition",
                    required=TRUE)
parser$add_argument("-i", "--indexcols", 
                    type="character",
                    help="comma separated list of columns that do not contain any counts eg. ensemblID, geneName, etc., ie., columns to be excluded from normalization by included in the output file.",
                    required=TRUE)
parser$add_argument("-x","--excludecols", 
                    type="character",
                    help="comma separated list of columns in the input that should be excluded from the output file.",
                    required=FALSE)
parser$add_argument("-t","--cpmthreshold", 
                    type="character", default="1",
                    help="cpm threshold (Default=1.0). Genes will cpm less than threshold are filtered out.",
                    required=FALSE)
parser$add_argument("-f","--mingroupfraction", 
                    type="character", default="0.5",
                    help="Fraction of samples per group that should meet the CPM threshold",
                    required=FALSE)
parser$add_argument("-o","--outfile", 
                    type="character",
                    help="name of outfile",
                    required=TRUE)


args <- parser$parse_args()


suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("tidyverse"))
debug=0

rawcountsmatrix=args$rawcountsmatrix
coldata=args$coldata
indexcols=unlist(strsplit(args$indexcols,","))
#print(indexcols)
if (length(args$excludecols)==0){
excludecols=c()
} else {
excludecols=unlist(strsplit(args$excludecols,","))
}
#print(excludecols)
outfile=args$outfil
cpmthreshold=as.numeric(args$cpmthreshold)
min_group_fraction=as.numeric(args$mingroupfraction)
outfile2=paste0(outfile,".antilog")


if(debug==1){
  rawcountsmatrix="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_raw_counts_counts.tsv"
  coldata="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_raw_counts_counts.coldata"
  indexcols=unlist(strsplit("ensemblID,gene_name,mRNA_length",","))
  excludecols=unlist(strsplit("596-7-2_p1",","))
  cpmthreshold=as.numeric("1.0")
  min_group_fraction=as.numeric("0.5")
  outfile="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_limmavoom_normalized_counts.tsv"
  outfile2="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_limmavoom_normalized_counts.antilog.tsv"
}


# read in raw counts
d=read.csv(rawcountsmatrix,header=TRUE,sep="\t",check.names = FALSE)

# remove excludecols, concate includecols into a single column and use it as index
d %>% select(-all_of(excludecols)) %>% 
  unite("geneID",all_of(indexcols),sep="##",remove=TRUE) %>% 
  column_to_rownames(.,var="geneID") -> e

# load coldata
as.data.frame(read.csv(coldata,header = TRUE,sep="\t")) %>%
  select(c("sample_name","condition")) -> cdata
  # column_to_rownames(.,var="sample_name") -> cdata
# change hyphen to underscore in conditions
cdata$condition = as.factor(gsub("-","_",cdata$condition))
design=model.matrix(~0+cdata$condition)

d0 <- DGEList(as.matrix(e))
d0 <- calcNormFactors(d0)

# apply cpm filters
conditions <- as.vector(unique(cdata$condition))
cpmd0 <- cpm(d0)
group = conditions[1]
print(group)
cpmsubset = cpmd0[,cdata[cdata$condition==group,]$sample_name]
nsamples = ncol(cpmsubset)
keep = !(rowSums(cpmsubset<cpmthreshold)/nsamples>min_group_fraction)
for (i in 2:length(conditions)){
  group=conditions[i]
  print(group)
  cpmsubset = cpmd0[,cdata[cdata$condition==group,]$sample_name]
  nsamples = ncol(cpmsubset)
  k = !(rowSums(cpmsubset<cpmthreshold)/nsamples>min_group_fraction)
  keep = (keep|k)
}
d <- d0[keep,]

# apply voom
v <- voom(as.matrix(d),design,plot=FALSE,normalize="quantile")


# Plot column sums according to size factor
# plot(d0$samples$norm.factors, d0$samples$lib.size)
# abline(lm(d0$samples$norm.factors ~ d0$samples$lib.size + 0))

# get normalized counts
logcounts <- v$E
as.data.frame(logcounts) %>% rownames_to_column(.,var="geneID") %>%
  separate(col="geneID",into=indexcols,sep="##",remove = TRUE) -> outdf

# write output
write.table(outdf,file=outfile,sep="\t",quote = FALSE,row.names = FALSE)

antilogcounts <- 2^logcounts
as.data.frame(antilogcounts) %>% rownames_to_column(.,var="geneID") %>%
  separate(col="geneID",into=indexcols,sep="##",remove = TRUE) -> outdf2
write.table(outdf2,file=outfile2,sep="\t",quote = FALSE,row.names = FALSE)
