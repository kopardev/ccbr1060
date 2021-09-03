library(ggplot2)
library(rlist)
library(tidyverse)
library(reshape2)
# setwd("~/Documents/Projects/ccbr1060/test")
options <- commandArgs(trailingOnly = TRUE)

# options <- c("596-11-1.detailed_counts_matrix.tsv",
#              "596-11-2.detailed_counts_matrix.tsv",
#              "596-11-Dox-1.detailed_counts_matrix.tsv",
#              "596-11-Dox-2.detailed_counts_matrix.tsv")

loadfile <- function(fn){
  sn=unlist(strsplit(basename(fn),"\\."))[1]
  df=read.csv(fn,sep="\t",header=TRUE,check.names = FALSE)
  coln=colnames(df)
  coln[1]="integration_window"
  colnames(df)=coln
  df$sample_name=sn
  return(df)
}

d=list.rbind(lapply(options[2:length(options)],loadfile))
# d=d[colSums(d[,2:(ncol(d)-1)])>10,]

wid=0
newcolnames=c()
for (i in colnames(d)){
  if(grepl("chr_shRNA",i)){
    wid=wid+1
    swid=sprintf("%02d", wid)
    newname=gsub(pattern = "chr_shRNA:", replacement = paste0("sw_",swid,"|"),i)
  }else{
    newname=i
  }
  newcolnames=c(newcolnames,newname)
}
colnames(d)=newcolnames

samples=unique(d$sample_name)
integration_windows=unique(d$integration_window)
sums=c()
# f=function(x) paste(as.character(strsplit(strsplit(as.character(x),split = "|",fixed = TRUE)[[1]][1],":")[[1]][2]),strsplit(as.character(x),split="|",fixed=TRUE)[[1]][2],sep="|")
iw=integration_windows[9]

plot_list = list()
i=0
pdf(options[1])
for (iw in integration_windows){
  d2=d[d$integration_window==iw,]
  d3=melt(d2,id.vars = c("sample_name","integration_window"))
  colnames(d3)=c("sample","integration_window","shRNA_window","counts")
  s=sum(d3$counts)
  sums=c(sums,s)
  
  if (s>0){
    i=i+1
    d3$shRNA_coordinate=d3$shRNA_window
    # d3$shRNA_coordinate = sapply(d3$shRNA_window,f)
    #avoid reordering
    # d3$shRNA_coordinate = factor(unique(d3$shRNA_coordinate),levels=unique(d3$shRNA_coordinate))
    p=ggplot(data=d3,aes(x=shRNA_coordinate,y=counts,fill=shRNA_coordinate)) +
      geom_col()+facet_wrap(~sample,nrow = 2) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      theme(plot.title = element_text(size=10))+ggtitle(iw)
    
    plot_list[[i]] = p
    # fn=paste0(iw,"_shRNA_integration_distribution.pdf")
    # pdf(fn)
    print(plot_list[[i]])
    # dev.off()
  }
}
dev.off()

