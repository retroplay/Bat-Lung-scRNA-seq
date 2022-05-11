library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)
setwd("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\4.dotplot\\002cytokine")
data<-read.table("Dotplotdata_all_Cytokine.txt",sep = "\t",header = T,stringsAsFactors = F)
cytokine <-read.table("002Cytokine.txt",sep = "\t",header = T,stringsAsFactors = F)

data<-data[which(data$id== "ATI" | data$id == "ATII" | data$id == "Ciliated cells" |
                   data$id == "Secretory cells" | data$id == "Fibroblasts" |data$id == "Endothelial cells" |
                   data$id == "Macrophages" | data$id == "B cells" | data$id == "T cells"),]
Data<-merge(data,cytokine,by.x = "features.plot",by.y = "Gene")

no0Data<-Data[which(Data$pct.exp != "0"),]

p<-ggplot(no0Data,aes(species,features.plot,fill =avg.exp.scaled, size=pct.exp))+
  geom_point(shape = 21,colour="black")+
  theme(axis.text.x = element_text(angle=45,vjust=0.6, hjust = 0.5),
        axis.text.y = element_text(face = "italic"),
        strip.background.x = element_rect(color="black",  size=1.5, linetype="solid"),
        strip.background.y = element_rect(color="black",  size=1.5, linetype="solid"),
        #panel.border = element_rect(color="black",fill =),
        #panel.background = element_rect(fill="white"),
        axis.line=element_line(color="black",size=1),
        strip.text.y.left=element_text(angle=0,face="bold"),
        panel.spacing.x = unit(0.5, "cm"),
        panel.spacing.y = unit(0.5, "cm")
  )+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red")+
  facet_grid(Function~id,scales = "free",space="free",switch = "y")+
  labs(x = "Species", y = "Cytokines")+
  scale_y_discrete(position = "right")


pdf("002Cytokines.pdf",width = 30,height = 25)
print(p)
dev.off()

write.table(no0Data,"002Cytokines_plotdata.txt",sep = "\t",row.names = F)


