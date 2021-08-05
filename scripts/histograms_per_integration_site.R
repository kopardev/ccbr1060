rm(list=ls())
setwd("~/Projects/ccbr1060")


library(ggplot2)
library(reshape2)

d=read.csv("detailed_counts.txt",sep="\t",header=TRUE,check.names = FALSE)
# View(d)
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

samples=unique(d$sample)
integration_windows=unique(d$integration_window)
sums=c()
f=function(x) paste(as.character(strsplit(strsplit(as.character(x),split = "|",fixed = TRUE)[[1]][1],":")[[1]][2]),strsplit(as.character(x),split="|",fixed=TRUE)[[1]][2],sep="|")
iw=integration_windows[3]

plot_list = list()
i=0
for (iw in integration_windows){
  d2=d[d$integration_window==iw,]
  d3=melt(d2,id.vars = c("sample","integration_window"))
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
  fn=paste0(iw,"_shRNA_integration_distribution.pdf")
  pdf(fn)
  print(plot_list[[i]])
  dev.off()
  }
}


iw=integration_windows[71]
d3=melt(d2,id.vars = c("sample","integration_window"))
colnames(d3)=c("sample","integration_window","shRNA_window","counts")

f=function(x) paste(as.character(sum(as.integer(strsplit(strsplit(strsplit(as.character(x),split = "|",fixed = TRUE)[[1]][1],":")[[1]][2],"-")[[1]]))/2),strsplit(as.character(x),split="|",fixed=TRUE)[[1]][2],sep="|")





d3$shRNA_coordinate = sapply(d3$shRNA_window,f)
ggplot(data=d3,aes(x=shRNA_coordinate,y=counts,fill=sample)) +
  geom_col()+facet_wrap(~sample,nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(iw)
  

x=d3$shRNA_window[1]
paste(as.character(sum(as.integer(strsplit(strsplit(strsplit(x,split = "|",fixed = TRUE)[[1]][1],":")[[1]][2],"-")[[1]]))/2),strsplit(x,split="|",fixed=TRUE)[[1]][2],sep="|")

df2 <- data.frame(supp=rep(c("VC", "OJ"), each=3),
                  dose=rep(c("D0.5", "D1", "D2"),2),
                  len=c(6.8, 15, 33, 4.2, 10, 29.5))
View(df2)



