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
                    help="two tab delimited columns.. samplename and condition",
                    required=FALSE)
parser$add_argument("-i", "--indexcols", 
                    type="character",
                    help="comma separated list of columns that do not contain any counts eg. ensemblID, geneName, etc., ie., columns to be excluded from normalization by included in the output file.",
                    required=TRUE)
parser$add_argument("-x","--excludecols", 
                    type="character",
                    help="comma separated list of columns in the input that should be excluded from the output file.",
                    required=FALSE)
parser$add_argument("-o","--outfile", 
                    type="character",
                    help="name of outfile",
                    required=TRUE)


args <- parser$parse_args()


suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("tidyverse"))
debug=0

rawcountsmatrix=args$rawcountsmatrix
coldata=args$coldata
indexcols=unlist(strsplit(args$indexcols,","))
excludecols=unlist(strsplit(args$excludecols,","))
outfile=args$outfil


if(debug==1){
  rawcountsmatrix="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_raw_counts_counts.tsv"
  coldata="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_raw_counts_counts.coldata"
  indexcols=unlist(strsplit("ensemblID,gene_name,mRNA_length",","))
  excludecols=unlist(strsplit("596-7-2_p1",","))
  outfile="/Volumes/CCBR/projects/ccbr1060/Hg38_shRNA_hybrid/HGHY2DRXY_analysis_v2/results/all_DESeq2_normalized_counts.tsv"
}


# read in raw counts
d=read.csv(rawcountsmatrix,header=TRUE,sep="\t",check.names = FALSE)

# remove excludecols, concate includecols into a single column and use it as index
d %>% select(-all_of(excludecols)) %>% 
  unite("geneID",all_of(indexcols),sep="##",remove=TRUE) %>% 
  column_to_rownames(.,var="geneID") -> e

if ( is.null(coldata) ){
  # convert count matrix df to SE object
  se <- SummarizedExperiment(list(counts=as.matrix(e)))
  # head(assay(se))
  
  # convert SE object to DESeqDataSet
  dds <- DESeqDataSet( se, design = ~ 1 )
  
} else {
  as.data.frame(read.csv(coldata,header = TRUE,sep="\t")) %>%
    select(c("sample_name","condition")) %>%
    column_to_rownames(.,var="sample_name") -> cdata
  # change hyphen to underscore in conditions
  cdata$condition = as.factor(gsub("-","_",cdata$condition))
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(e),
                                colData = cdata,
                                design = ~ condition)
}

#Estimate size factors
dds <- estimateSizeFactors( dds )
# sizeFactors(dds)

# Plot column sums according to size factor
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

# get normalized counts
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
as.data.frame(logcounts) %>% rownames_to_column(.,var="geneID") %>%
  separate(col="geneID",into=indexcols,sep="##",remove = TRUE) -> outdf

# write output
write.table(outdf,file=outfile,sep="\t",quote = FALSE,row.names = FALSE)
