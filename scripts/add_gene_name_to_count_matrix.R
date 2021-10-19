#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-r", "--rawcountsmatrix", 
                    type="character", 
                    help="file with raw counts matrix with gene_id column",
                    required=TRUE)
parser$add_argument("-g", "--gtf", 
                    type="character", 
                    help="GTF file with gene_name and gene_id",
                    required=TRUE)
parser$add_argument("-o","--outfile", 
                    type="character",
                    help="count matrix with gene_name column included",
                    required=TRUE)


args <- parser$parse_args()

library("rtracklayer")
library("tidyverse")
gtf2<-rtracklayer::import(args$gtf)
gtf2<-as.data.frame(gtf2)
unique(data.frame(gene_id=gtf2$gene_id,gene_name=gtf2$gene_name)) %>% drop_na() -> lookuptable

in_df=read.csv(args$rawcountsmatrix,sep="\t",header=TRUE)

out_df=merge(in_df,lookuptable,by=c("gene_id"))

write.table(out_df,file=args$outfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
