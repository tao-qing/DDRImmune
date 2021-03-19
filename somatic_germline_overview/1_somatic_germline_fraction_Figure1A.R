##### somatic_germline_overlap.R #####
# Kuan-lin Huang @ WashU 2018
# Find overlap of genes/variants for somatic/germline variants
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_overlap/")

overlapsplrmhyper<-read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",stringsAsFactors = F)[,1]

#all core DDR genes
ddr=read.csv("../../Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/DDR_Pathways.csv",h=T)
ddrgene=unique(as.character(unlist(ddr)))[-1]

source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../load_somatic.R")
rmhypermutator=TRUE
pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%overlapsplrmhyper),]
somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[which(somatic_likelyfunctional_driver$bcr_patient_barcode%in%overlapsplrmhyper),]
somatic_likelyfunctional_driver$cancer=clin$type[sapply(somatic_likelyfunctional_driver$bcr_patient_barcode,function(x)which(clin$bcr_patient_barcode==x))]
###################################################
####Figure 1A Distribution of DDR mutated genes####
###################################################
allspl=table(clin$type[clin$bcr_patient_barcode%in%overlapsplrmhyper])
freqMat=as.data.frame(cbind(Cancer=names(allspl),SampleSize=allspl))
freqMat$Cancer=as.character(freqMat$Cancer)
freqMat$SampleSize=as.numeric(as.matrix(freqMat$SampleSize))

tmp0=pathVarP[which(pathVarP$HUGO_Symbol%in%ddrgene),]
tmp0=tmp0[-which(duplicated(tmp0$bcr_patient_barcode)),]
gfreq=table(tmp0$cancer)

#percentage of germline DDR carriers across all sample
length(unique(tmp0$bcr_patient_barcode))/9738
length(unique(tmp0$bcr_patient_barcode[tmp0$cancer=="OV"]))/386

tmp1=somatic_likelyfunctional_driver[which(somatic_likelyfunctional_driver$Hugo_Symbol%in%ddrgene),]
tmp1=tmp1[-which(duplicated(tmp1$bcr_patient_barcode)),]
sfreq=table(tmp1$cancer)

#percentage of somatic DDR carriers across all sample
length(unique(tmp1$bcr_patient_barcode))/9738

length(unique(tmp1$bcr_patient_barcode[tmp1$cancer=="UCEC"]))/478

freqMat$Germline=as.numeric((gfreq[freqMat$Cancer]/freqMat$SampleSize)*100)
freqMat$Somatic=as.numeric((sfreq[freqMat$Cancer]/freqMat$SampleSize)*100)

library(reshape2)
plotMat=melt(freqMat,id=c("Cancer","SampleSize"))
plotMat$Label=paste0(plotMat$Cancer," (",plotMat$SampleSize,")")
tmp3=plotMat[plotMat$variable=="Germline",]
or=tmp3$Label[order(tmp3$value,decreasing=T)]
plotMat$Label=factor(plotMat$Label,levels=or)
colnames(plotMat)=gsub("variable","Type",colnames(plotMat))

p= ggplot(plotMat,aes(y=as.numeric(as.matrix(value)),x=Label,fill=Type)) + geom_bar(stat="identity",position='dodge')
#p= p+ geom_text(aes(label=Counts), hjust=1, size=3)
p = p  + theme_bw() 
p = p  + theme(legend.position = c(0.85, 0.85),axis.text.x = element_text(colour="black", size=14, angle=90, vjust = 0.5,hjust = 0.95), axis.text.y = element_text(colour="black", size=14,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0,size=16,face="bold"),axis.title=element_text(size=14,face="bold"),panel.border = element_blank(),axis.line= element_line(color='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p = p + scale_fill_manual("Type", values = c("Germline" = "blue","Somatic"="orange"))
p = p + labs(title="Samples affected by DDR mutations",y="Percentage",x="")
p

fn = "./out/Figure1E_somaticgermline_frequency_of_DDR_affected_samples.pdf"
ggsave(file=fn, width=8, h =3.5, useDingbats=FALSE)

