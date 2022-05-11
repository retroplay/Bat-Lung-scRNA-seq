library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)
setwd("D:\\1_precision_M_B\\2_project_end\\bat_mink_atlas_plan\\bat_lung\\Revision\\4.dotplot\\001receptor")
data<-read.table("Dotplotdata_all_Receptor.txt",sep = "\t",header = T,stringsAsFactors = F)
virus_receptor <-read.table("001Virus_receptor.txt",sep = "\t",header = T,stringsAsFactors = F)

data<-data[which(data$id== "ATI" | data$id == "ATII" | data$id == "Ciliated cells" |
                   data$id == "Secretory cells" | data$id == "Fibroblasts" |data$id == "Endothelial cells" |
                   data$id == "Macrophages" | data$id == "B cells" | data$id == "T cells"),]
Data<-merge(data,virus_receptor,by.x = "features.plot",by.y = "Receptor")

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
        strip.text.y.left=element_text(angle=0,face="bold")
  )+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red")+
  facet_grid(Respiratory~id,scales = "free",space="free",switch = "y")+
  labs(x = "Species", y = "Receptors")+
  scale_y_discrete(position = "right")
 
pdf("001Virus_receptor.pdf",width = 25,height = 15)
print(p)
dev.off()

write.table(no0Data,"001Virus_receptor_plotdata.txt",sep = "\t",row.names = F)


###DPW
p <- ggplot(data=data,mapping=aes_string(x='species',y='features.plot'))+
  geom_point(mapping=aes_string(size='pct.exp',color='avg.exp'),show.legend = TRUE)+
  #  scale_size(range=c(2,15),limits = c(min(data$pct.exp),max(data$pct.exp)),breaks = NULL)+
  scale_color_gradient(low="blue",high = "red",breaks=NULL)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent Expressed'),color=guide_legend(title = 'Average Expression')) +
  labs(x = 'features.plot', y = 'species',color='Average Expression',size="Percent Expressed",title='')+
  theme_bw()+
  coord_flip()+theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust = 0.5))+
  facet_grid(id~species,shrink = F,scales = "free_y",space = "free_y",switch = "y",as.table =F)+
  theme(strip.text.x = element_text(size = 12,face="bold"),strip.text.y.left=element_text(angle=0,size=15,face="bold",color="red"))

