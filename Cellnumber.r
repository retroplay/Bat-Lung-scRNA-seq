library(ggthemes)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(tidyr)
setwd("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\7.immunecells\\")

metapath<-c("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\7.immunecells\\metadata\\")

meta<-c()
metadata<-c()
for (i in c("Bat","Cat","Pangolin","Tiger")) {
  data<-read.table(paste0(metapath,"metadata_",i,".txt"),sep="\t",header = T)
  data<-data[which(data$NewCelltype== "ATI" | data$NewCelltype == "ATII" | data$NewCelltype == "Ciliated cells" |
                     data$NewCelltype == "Secretory cells" | data$NewCelltype == "Fibroblasts" |data$NewCelltype == "Endothelial cells" |
                     data$NewCelltype == "Macrophages" | data$NewCelltype == "B cells" | data$NewCelltype == "T cells"),]
  
  data<-data[,c("NewCelltype","species")]
  meta<-rbind(meta,data)
  data1<-mutate(as.data.frame(table(data.frame(id=data$NewCelltype))))
  colnames(data1)<-c("NewCelltype","Number")
  data1$species<-i
  metadata<-rbind(metadata,data1)
  
}

write.table(metadata,"cellnumber_4species.txt",sep = "\t",row.names = F)

colors<-c(brewer.pal(9,"Paired"))
color4<-c("#DE2D26","#C51B8A","#31A354","#3182BD")

######celltype
p1<-ggplot(meta,aes(NewCelltype,fill=NewCelltype))+geom_bar(position="dodge")+
  theme_set(theme_bw())+
  scale_fill_manual(values=colors)+
  guides(fill=guide_legend(title=NULL))+
  ggtitle("Celltypes statistics of each species")+ 
  theme(axis.title = element_blank(), plot.title = element_text(size = 20))+
  facet_grid(rows = vars(species),scales = "free_y")
  
p2<-ggplot(meta,aes(NewCelltype,fill=species))+geom_bar(position="dodge")+
  theme_set(theme_bw())+
  scale_fill_manual(values=color4)+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title = element_blank(), panel.grid.major=element_line(colour=NA))


p3<-ggplot(meta,aes(species,fill=NewCelltype))+geom_bar(position="fill")+
  theme_set(theme_bw())+
  scale_fill_manual(values=colors)+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title = element_blank(), panel.grid.major=element_line(colour=NA))

p4<-ggplot(meta,aes(NewCelltype,fill=NewCelltype))+geom_bar(position="dodge")+
  theme_set(theme_bw())+
  scale_fill_manual(values=colors)+
  guides(fill=guide_legend(title=NULL))+
  facet_wrap(~species) +
  theme(axis.title = element_blank(), panel.grid.major=element_line(colour=NA))


pdf("cellnumber_4species.pdf", width = 25, height = 20)
print(p1+p2+p3+p4)
dev.off()









