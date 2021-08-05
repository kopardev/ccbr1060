rm(list=ls())
setwd("~/Projects/ccbr1060")


library(ggplot2)
library(reshape2)




shRNA_coordinates=read.csv("shRNA_coordinates.txt",header = FALSE)
# fviz_nbclust(as.data.frame(shRNA_coordinates$V1), kmeans, method = "wss",k.max = 10)
# 
# 
# gap_stat <- clusGap(shRNA_coordinates,
#                     FUN = kmeans,
#                     nstart = 25,
#                     K.max = 10,
#                     B = 50)
# fviz_gap_stat(gap_stat)
km=list()
for (i in seq(1,15)){
  km[[i]]=kmeans(shRNA_coordinates$V1,centers=i,nstart=20,iter.max = 100000)
}

get_count <- function(y) 
  c(label=length(y), y=median(y))


for (i in seq(1,15)){
  for (j in seq(1,i)){
    fn=as.integer(fivenum(shRNA_coordinates[km[[i]]$cluster==j,1]))
    x=sprintf("%d %d %d %d %d %d %d",i,j,fn[[1]],fn[[2]],fn[[3]],fn[[4]],fn[[5]])
    print(x)
  }
  shRNA_coordinates$cluster=paste0("C",km[[i]]$cluster)
  p=ggplot(shRNA_coordinates,aes(x=cluster,y=V1))+geom_boxplot(aes(col=cluster))+stat_summary(fun.data=get_count, geom="text", vjust=-0.5, col="blue", hjust=0.5)+coord_flip()+ylab("shRNA_coordinate")+theme(legend.position = "none")+theme_classic()
  pdf(paste0("shRNA_coordinate_nclusters_",as.character(i),".pdf"))
  print(p)
  dev.off()
}


# 13 clusters give best coverage accross shRNA
# save sprintf output to tmp file and 
# awk -F"\"" '{print $2}' tmp|grep ^13|awk '{print $2,$6-$4,$4,$6}'|sort -k3,3n                                       
# 11 0 636 636
# 1 111 2192 2303
# 8 73 3002 3075
# 10 64 3176 3240
# 5 74 3366 3440
# 6 131 3548 3679
# 3 43 3803 3846
# 13 80 4474 4554
# 12 178 5416 5594
# 2 254 5780 6034
# 9 307 6349 6656
# 4 195 7446 7641
# 7 367 8040 8407
# picking quantile values for coordinates
