rm(list=ls())
library(Seurat)
library(dplyr)
library(patchwork)
library(devtools)
library(SingleR)
library(mindr)
library(Matrix)
library(magrittr)
library(hdf5r)
library(SeuratObject)

setwd('F:\\HLQ\\cell_cell_communication')
mypvals <- read.table("out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
mymeans <- read.table("out/means.txt",header = T,sep = "\t",stringsAsFactors = F) 
kp = grepl(pattern = "B_cell", colnames(mypvals)) 

table(kp)
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,2,5,6,8,9),pos)]
choose_means <- mymeans[,c(c(1,2,5,6,8,9),pos)]

logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.05, 1, sum) 
logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]

choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

library(tidyverse)
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = meansdf$interacting_pair,
                      CC = meansdf$variable,
                      means = meansdf$value)

pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = pvalsdf$interacting_pair,
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)

pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >0))$means)
head(pldf)
pcc =  pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 1  )+ 
  theme_bw()+ 
  # scale_color_manual(values = rainbow(100))+
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8))
pcc

mypvals <- read.table("out/pvalues.txt",
                      header = T,sep = "\t",stringsAsFactors = F)
mymeans <- read.table("out/means.txt",
                      header = T,sep = "\t",stringsAsFactors = F) 

colnames(mypvals)[12:ncol(mypvals)] 
sm = as.data.frame(
  do.call(rbind,
          lapply( 12:ncol(mypvals) , function(i){
            return(c( strsplit(colnames(mypvals)[i],'\\.')[[1]],
                      sum(mypvals[,i] <0.05)))
          }))
)
head(sm)
colnames(sm)=c('SOURCE' ,'TARGET' ,'count')
sm$count = as.numeric( sm$count )
sm
write.table(sm,file = 'count_network.txt',
            sep = '\t',
            quote = F,row.names = F)



library(reshape2)
sm_df =dcast(as.data.frame(sm),SOURCE~TARGET )
sm_df[is.na(sm_df)]=0
sm_df

rownames(sm_df) = sm_df[,1]
sm_df = sm_df[,-1]#删除第一列
p1=pheatmap::pheatmap(sm_df,display_numbers = T)
p2=pheatmap::pheatmap(log(sm_df+1),display_numbers = T)
library(patchwork)
library(cowplot)
library(ggplotify)
as.ggplot(p1) + as.ggplot(p2)

ggplot2::ggsave('heatmap.pdf',width = 10,height = 5)

rm(list = ls())
library(psych)
library(qgraph)
library(igraph)
library(purrr) # map function

netf<- read.table("count_network.txt", header = T)
mynet <- netf
head(mynet)

net<- graph_from_data_frame(mynet)
plot(net)

allcolour=c('#a6cee3','#1f78b4','#b2df8a','#33a02c' )

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))

E(net)$width  <- E(net)$count/40 
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 

net2 <- net  # 复制一份备用

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
} 

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.size = 30,
     vertex.label.cex=.7) 


dev.off()
length(unique(mynet$SOURCE)) 
par(mfrow=c(2,3), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1,
       edge.curved=0.4,
       edge.width=3,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1,
       vertex.size = 30)
  
} 


for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count 
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.5, 
       edge.curved=0.5,
       edge.width=3,
       edge.label = E(net1)$count, 
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1.2
  ) 
}
