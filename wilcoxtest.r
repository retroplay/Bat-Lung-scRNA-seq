library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(dplyr)
setwd("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\10.pvalue")

###receptor
datasource<-c("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\10.pvalue\\Receptor\\")
datalist<-list()
i=1
S<-c("Bat")
datalist[[i]] <- read.table(paste0(datasource,S,"_expmatrix_receptors.txt"),sep = "\t",stringsAsFactors = F)
alldata<-as.data.frame(t(datalist[[i]]))
alldata$key<-rownames(alldata)

Species <- c("Cat","Pangolin","Tiger")
for (S in Species) {
  i=i+1
  datalist[[i]] <- read.table(paste0(datasource,S,"_expmatrix_receptors.txt"),sep = "\t",stringsAsFactors = F)
  df<-as.data.frame(t(datalist[[i]]))
  df$key<-rownames(df)
  alldata<-full_join(alldata,df,by = "key")
}

alldata<-as.data.frame(t(alldata))
colnames(alldata)<-alldata["key",]
alldata<-alldata[-which(rownames(alldata)=="key"),]

write.table(alldata,"expmatrix_receptors.txt",sep = "\t")

alldata<-read.table("expmatrix_receptors.txt",sep = "\t",stringsAsFactors = F)

#setwd("Boxplot")
alldata<-alldata[which(alldata$NewCelltype == "ATI" | alldata$NewCelltype == "ATII" | alldata$NewCelltype == "Ciliated cells" |
                 alldata$NewCelltype == "Secretory cells" | alldata$NewCelltype == "Fibroblasts" |alldata$NewCelltype == "Endothelial cells" |
                 alldata$NewCelltype == "Macrophages" | alldata$NewCelltype == "B cells" | alldata$NewCelltype == "T cells"),]
alldata$NewCelltype<-factor(alldata$NewCelltype,levels = c("ATI","ATII","Ciliated cells","Secretory cells","Fibroblasts",
                                                     "Endothelial cells", "Macrophages","B cells", "T cells"))

alldata$species<-factor(alldata$species,levels = c("Bat","Cat","Pangolin","Tiger"))
compaired <- list(c("Bat", "Cat"),c("Bat","Pangolin"),c("Bat","Tiger"))
color4<-c("#DE2D26","#C51B8A","#31A354","#3182BD")

receptors<-c("Itgb5","Anpep","Cd86","Cdhr3","Rpsa","Ace2","Scarb1","Nrp1","Axl")
if(!dir.exists("Boxplot_Receptors")){dir.create("Boxplot_Receptors")}
setwd("Boxplot_Receptors")

for (i in unique(alldata$NewCelltype)) {
  tmp<-alldata[which(alldata$NewCelltype==i),]
  for (j in receptors) {
    tmp1<-data.frame(gne=j,exp=tmp[,which(colnames(tmp) == j)],species=tmp$species,celltype=tmp$NewCelltype)
   p<-ggplot(data = tmp1,mapping = aes(x = species,y = exp,fill=species)) + geom_boxplot(outlier.colour = 1)+
      scale_fill_manual(values=color4)+
      ggtitle(paste(i,j,sep = "_"))+theme_set(theme_bw())+
      theme(plot.background = element_rect(colour = NA),panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA),axis.line = element_line(colour = "black"))+
      geom_signif(comparisons = compaired,step_increase = 0.05,map_signif_level = T,test = wilcox.test)
    pdf(paste0("Boxplot_",i,"_",j,".pdf"),width = 4,height = 4)
    print(p)
    dev.off()
  }
}


plot<-list()
pcell<-list()
x=1
for (i in levels(alldata$NewCelltype)) {
  tmp<-alldata[which(alldata$NewCelltype==i),]
  plot<-list()
  y=1
  for (j in receptors) {
    tmp1<-data.frame(gne=j,exp=tmp[,which(colnames(tmp) == j)],species=tmp$species,celltype=tmp$NewCelltype)
plot[[y]]<-ggplot(data = tmp1,mapping = aes(x = species,y = exp,fill=species)) + geom_boxplot(outlier.colour = 1)+
      scale_fill_manual(values=color4)+theme_set(theme_bw())+
      theme(plot.background = element_rect(colour = NA),panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA),axis.line = element_line(colour = "black"))+
      geom_signif(comparisons = compaired,step_increase = 0.05,map_signif_level = T,test = wilcox.test)+
      theme(legend.position = 'none', 
           axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_text(face = "italic"),
           plot.title =element_blank())+
      labs(y=j)
    y=y+1
  }
  pcell[[x]]<-ggarrange(plot[[1]],plot[[2]],plot[[3]],
            plot[[4]],plot[[5]],plot[[6]],
            plot[[7]],plot[[8]],plot[[9]],
            ncol = 1, nrow = 9,labels = i)
  x=x+1
}
pall<-ggarrange(pcell[[1]],pcell[[2]],pcell[[3]],
                pcell[[4]],pcell[[5]],pcell[[6]],
                pcell[[7]],pcell[[8]],pcell[[9]],
                ncol = 9, nrow = 1)
pdf("Boxplot_9receptor.pdf",width = 20,height = 22)
print(pall)
dev.off()





###Cytokines

setwd("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\10.pvalue")
datasource<-c("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\10.pvalue\\Cytokines\\")
datalist<-list()
i=1
S<-c("Bat")
datalist[[i]] <- read.table(paste0(datasource,S,"_expmatrix_cytokines.txt"),sep = "\t",stringsAsFactors = F)
alldata<-as.data.frame(t(datalist[[i]]))
alldata$key<-rownames(alldata)

Species <- c("Cat","Pangolin","Tiger")
for (S in Species) {
  i=i+1
  datalist[[i]] <- read.table(paste0(datasource,S,"_expmatrix_cytokines.txt"),sep = "\t",stringsAsFactors = F)
  df<-as.data.frame(t(datalist[[i]]))
  df$key<-rownames(df)
  alldata<-full_join(alldata,df,by = "key")
}

alldata<-as.data.frame(t(alldata))
colnames(alldata)<-alldata["key",]
alldata<-alldata[-which(rownames(alldata)=="key"),]

write.table(alldata,"expmatrix_cytokines.txt",sep = "\t")

alldata<-read.table("expmatrix_cytokines.txt",sep = "\t",stringsAsFactors = F)

#setwd("Boxplot")
alldata<-alldata[which(alldata$NewCelltype == "ATI" | alldata$NewCelltype == "ATII" | alldata$NewCelltype == "Ciliated cells" |
                         alldata$NewCelltype == "Secretory cells" | alldata$NewCelltype == "Fibroblasts" |alldata$NewCelltype == "Endothelial cells" |
                         alldata$NewCelltype == "Macrophages" | alldata$NewCelltype == "B cells" | alldata$NewCelltype == "T cells"),]
alldata$NewCelltype<-factor(alldata$NewCelltype,levels = c("ATI","ATII","Ciliated cells","Secretory cells","Fibroblasts",
                                                           "Endothelial cells", "Macrophages","B cells", "T cells"))

alldata$species<-factor(alldata$species,levels = c("Bat","Cat","Pangolin","Tiger"))
compaired <- list(c("Bat", "Cat"),c("Bat","Pangolin"),c("Bat","Tiger"))
color4<-c("#DE2D26","#C51B8A","#31A354","#3182BD")

cytokines<-c("Osmr", "Lif","Lifr","Ifngr2", "Ifnar1", "Ifnar2")
if(!dir.exists("Boxplot_Cytokines")){dir.create("Boxplot_Cytokines")}
setwd("Boxplot_Cytokines")
for (i in unique(alldata$NewCelltype)) {
  tmp<-alldata[which(alldata$NewCelltype==i),]
  for (j in cytokines) {
    tmp1<-data.frame(gne=j,exp=tmp[,which(colnames(tmp) == j)],species=tmp$species,celltype=tmp$NewCelltype)
    p<-ggplot(data = tmp1,mapping = aes(x = species,y = exp,fill=species)) + geom_boxplot(outlier.colour = 1)+
      scale_fill_manual(values=color4)+
      ggtitle(paste(i,j,sep = "_"))+theme_set(theme_bw())+
      theme(plot.background = element_rect(colour = NA),panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA),axis.line = element_line(colour = "black"))+
      geom_signif(comparisons = compaired,step_increase = 0.05,map_signif_level = T,test = wilcox.test)
    pdf(paste0("Boxplot_",i,"_",j,".pdf"),width = 4,height = 4)
    print(p)
    dev.off()
  }
}

plot<-list()
pcell<-list()
x=1
for (i in levels(alldata$NewCelltype)) {
  tmp<-alldata[which(alldata$NewCelltype==i),]
  plot<-list()
  y=1
  for (j in cytokines) {
    tmp1<-data.frame(gne=j,exp=tmp[,which(colnames(tmp) == j)],species=tmp$species,celltype=tmp$NewCelltype)
    plot[[y]]<-ggplot(data = tmp1,mapping = aes(x = species,y = exp,fill=species)) + geom_boxplot(outlier.colour = 1)+
      scale_fill_manual(values=color4)+theme_set(theme_bw())+
      theme(plot.background = element_rect(colour = NA),panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA),axis.line = element_line(colour = "black"))+
      geom_signif(comparisons = compaired,step_increase = 0.05,map_signif_level = T,test = wilcox.test)+
      theme(legend.position = 'none', 
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "italic"),
            plot.title =element_blank())+
      labs(y=j)
    y=y+1
  }
  pcell[[x]]<-ggarrange(plot[[1]],plot[[2]],plot[[3]],
                        plot[[4]],plot[[5]],plot[[6]],
                        ncol = 1, nrow = 9,labels = i)
  x=x+1
}
pall<-ggarrange(pcell[[1]],pcell[[2]],pcell[[3]],
                pcell[[4]],pcell[[5]],pcell[[6]],
                pcell[[7]],pcell[[8]],pcell[[9]],
                ncol = 9, nrow = 1)
pdf("Boxplot_6cytokine.pdf",width = 20,height = 22)
print(pall)
dev.off()









alldata$color<-as.character(alldata$species) 
alldata$color[which(alldata$color == "Bat")]<-c("#DE2D26")
alldata$color[which(alldata$color == "Cat")]<-c("#C51B8A")
alldata$color[which(alldata$color == "Pangolin")]<-c("#31A354")
alldata$color[which(alldata$color == "Tiger")]<-c("#3182BD")
alldata$color<-factor(alldata$color,levels = c("#DE2D26","#C51B8A","#31A354","#3182BD"))
