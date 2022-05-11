###载入需要的R包
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

setwd("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\9.figure\\fig1c")
data<-readRDS("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\BT_NewCelltype.rds")

DimPlot(data,reduction = "umap",group.by = "NewCelltype")
data$NewCelltype<-as.character(data$NewCelltype)
data$NewCelltype[which(data$NewCelltype %in% c("Brush cells","Mesothelial cells","SMCs","Unknown"))] <- "Others"
data$NewCelltype<-factor(data$NewCelltype,levels = c("ATI","ATII","Ciliated cells","Secretory cells","Fibroblasts","B cells",
                                                     "T cells","Macrophages","Endothelial cells","Others"))

colors<-as.character(brewer.pal(10,"Paired"))
pdf("fig1b.pdf",width = 8,height = 7)
DimPlot(data,reduction = "umap",group.by = "NewCelltype",cols = colors,label = T)
dev.off()

data$NewCelltype<-fct_rev(data$NewCelltype)
colors<-c("#6A3D9A","#CAB2D6","#FF7F00","#FDBF6F","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3")
marker<-c("Ager","Pdpn","Scnn1g","Etv5","Slc34a2","Ccdc39","Dnah6","Itga8","Pdgfra","Scgb3a2","Ebf1","Il7r","Itk","Mrc1","Msr1","Pecam1","Vwf","F8","Ace","Ace2")
DefaultAssay(data)<-"RNA"

p<-list()
j=1
for (i in marker) {
  p[[j]]<- VlnPlot(data, features = i,cols = colors,group.by = "NewCelltype",pt.size = 0,y.max = 6)+
                  coord_flip()+
                  theme(legend.position = 'none', 
                  axis.text.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text(hjust = 0.5,face= "italic",angle=45))+
                  guides(scale = "none")
  j=j+1
}

pdf("fig1c.pdf",width = 18,height = 10)
ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],
          p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],
          p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],
          p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],
          ncol = 20, nrow = 1)
#ggarrange(plotlist = p,ncol = 20, nrow = 1)
dev.off()

