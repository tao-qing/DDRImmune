setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_overlap/")
### dependencies ###
source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../load_somatic.R")

library(dplyr)
detach(package:plyr)
#core DDR pathway
ddr=as.data.frame(read.csv("../../Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/DDR_Pathways.csv",h=T))

ddr_list=list()
ddr_list$Base_Excision_Repair=unique(as.character(ddr$Base.Excision.Repair..BER.))[-1]
ddr_list$Nucleotide_Excision_Repair=unique(as.character(ddr$Nucleotide.Excision.Repair..NER..including.TC.NER.and.GC.NER..))[-1]
ddr_list$Mismatch_Repair=unique(as.character(ddr$Mismatch.Repair..MMR.))[-1]
ddr_list$Fanconi_Anemia=unique(as.character(ddr$Fanconi.Anemia..FA.))[-1]
ddr_list$Homologous_Recomination=unique(as.character(ddr$Homologous.Recomination..HR.))[-1]
ddr_list$Nonhomologous_End_Joining=unique(as.character(ddr$Non.homologous.End.Joining..NHEJ.))[-1]
ddr_list$Direct_Repair=unique(as.character(ddr$Direct.Repair..DR.))[-1]
ddr_list$Translesion_Synthesis=unique(as.character(ddr$Translesion.Synthesis..TLS.))[-1]
ddr_list$Damage_Sensor=unique(as.character(ddr$Damage.Sensor.etc.))[-1]

#cancer samples
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]

pathVarP=pathVarP%>%filter(bcr_patient_barcode%in%samples)

somatic_likelyfunctional_driver=somatic_likelyfunctional_driver%>%filter(bcr_patient_barcode%in%samples)

clin=clin%>%filter(bcr_patient_barcode%in%samples)
cancerType=table(clin$type)

somatic_likelyfunctional_driver<-somatic_likelyfunctional_driver%>%left_join(.,clin[,c("bcr_patient_barcode","type")],by="bcr_patient_barcode")%>%rename(cancer=type)


mat=NULL
#germline somatic mutation by sample
for(ca in unique(clin$type)){
  for(path in unique(names(ddr_list))){
    gnum=length(unique(pathVarP$bcr_patient_barcode[pathVarP$cancer==ca & pathVarP$HUGO_Symbol%in%ddr_list[[path]]]))
    snum=length(unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$cancer==ca & somatic_likelyfunctional_driver$Hugo_Symbol%in%ddr_list[[path]]]))
    gfreq=round(length(unique(pathVarP$bcr_patient_barcode[pathVarP$cancer==ca & pathVarP$HUGO_Symbol%in%ddr_list[[path]]]))/as.numeric(cancerType[ca]),digits = 4)
    sfreq=round(length(unique(somatic_likelyfunctional_driver$bcr_patient_barcode[ somatic_likelyfunctional_driver$cancer==ca & somatic_likelyfunctional_driver$Hugo_Symbol%in%ddr_list[[path]]]))/as.numeric(cancerType[ca]),digits = 4)
    mat=cbind(mat,rbind(Cancer=ca,Pathway=path,germline=gnum,somatic=snum,gFreq=gfreq,sFreq=sfreq))
  }
}

mat=as.data.frame(t(mat))
mat$germline=as.numeric(as.character(mat$germline))
mat$somatic=as.numeric(as.character(mat$somatic))
mat$gFreq=as.numeric(as.character(mat$gFreq))*100
mat$sFreq=as.numeric(as.character(mat$sFreq))*100


#spath=tapply(mat$sFreq,mat$Pathway,mean)
#rmpath=intersect(names(gpath)[gpath==0],names(spath)[spath==0])

#mat=mat[-which(mat$Pathway%in%rmpath),]

library(ggplot2)
col=brewer.pal(n = 8, name = "RdYlBu")
#heatmap germline
mat$Pathway=gsub("_"," ",as.character(mat$Pathway))

gpathmean=tapply(mat$gFreq,mat$Pathway,mean)
pathIndex=names(gpathmean)[order(gpathmean,decreasing = T)]
mat$Pathway=factor(mat$Pathway,levels=pathIndex)

gcancermean=tapply(mat$gFreq,mat$Cancer,mean)
cancerIndex=names(gcancermean)[order(gcancermean,decreasing = F)]
mat$Cancer=factor(mat$Cancer,levels=cancerIndex)

gpath=tapply(mat$gFreq,mat$Pathway,mean)
gmat=mat[-which(mat$Pathway%in%names(gpath)[gpath==0]),]

p = ggplot(data=gmat)
p = p + geom_tile(data=gmat,aes(y=Cancer, x=Pathway, fill= as.numeric(as.matrix(gFreq))), linetype="blank") + scale_fill_gradientn(name= "% of samples with mutation\n", colours=c("white",col[3],col[1]), na.value=NA, limit=c(0,20))
p = p + geom_text(data=gmat,aes(y=Cancer, x=Pathway, label = round(gFreq,digits = 1)), color="black", size=3)
p = p  + theme_bw() + theme_nogrid() +
  theme(legend.position = "none",axis.text.x = element_text(colour="black", size=12, angle=30, hjust = 0.95), axis.text.y = element_text(colour="black", size=10,hjust = 0.95),axis.title=element_text(size=14,face="bold"),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + labs(title="Germline",y="",x = "")+guides(fill = guide_legend(title.position = "top"))
p
fn = './out/germlineDDR_heatmap.pdf'
ggsave(fn,height=6,width=3.5)


#heatmap somatic
spathmean=tapply(mat$sFreq,mat$Pathway,mean)
pathIndex=names(spathmean)[order(spathmean,decreasing = T)]
mat$Pathway=factor(mat$Pathway,levels=pathIndex)

scancermean=tapply(mat$sFreq,mat$Cancer,mean)
cancerIndex=names(scancermean)[order(scancermean,decreasing = F)]
mat$Cancer=factor(mat$Cancer,levels=cancerIndex)

spath=tapply(mat$sFreq,mat$Pathway,mean)
smat=mat[-which(mat$Pathway%in%names(spath)[spath==0]),]

#smat$Pathway=gsub("_"," ",as.character(smat$Pathway))

p = ggplot(data=smat)
p = p + geom_tile(data=smat,aes(y=Cancer, x=Pathway, fill= as.numeric(as.matrix(sFreq))), linetype="blank") + scale_fill_gradientn(name= "% of samples\nwith mutation\n", colours=c("white",col[3],col[1]), na.value=NA, limit=c(0,20))
p = p + geom_text(data=smat,aes(y=Cancer, x=Pathway, label = round(sFreq,digits = 1)), color="black", size=3)
p = p  + theme_bw() + theme_nogrid() +
  theme(legend.position = "right",axis.text.x = element_text(colour="black", size=12, angle=30, hjust = 0.95), axis.text.y = element_text(colour="black", size=10,hjust = 0.95),axis.title=element_text(size=14,face="bold"),axis.ticks = element_blank())
p = p + labs(title="Somatic",y="",x = "")+guides(fill = guide_legend(title.position = "top"))
p
fn = './out/somaticDDR_heatmap.pdf'
ggsave(fn,height=6,width=5)

