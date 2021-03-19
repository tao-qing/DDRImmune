#we try to analyze how combined DDR mutation and neoantigens cause immune response 
#Neoantigen hotspots, from driver genes, 
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")

#################
####Load data####
#################
library(data.table)
library(dplyr)
source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../load_somatic_rmhypermut.R")
rmhypermutator=TRUE
label="_rmhypermutator_"

library(ggplot2)
library(dplyr)
#short name
xlabel=scale_y_discrete(position = "left",labels=c("Aneuploidy_Score" = "Aneuploidy", "BCR_Evenness" = "BCR_E","BCR_Richness" = "BCR_R","BCR_Shannon"="BCR_S","CTA_Score"="CTA","Homologous_Recombination_Defects"="HRD","TIL_Regional_Fraction"="TIL","Leukocyte_Fraction"="Leukocyte","Macrophage_Regulation"="Macrophage","Lymphocyte_Infiltration_Signature_Score"="Lymphocyte","IFN_gamma_Response"="IFN","TGF_beta_Response"="TGF_beta","Th1_Cells"="Th1","Th2_Cells"="Th2","Th17_Cells"="Th17"))

th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black", size=12, angle = 60,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

#cancer samples
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]
pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%samples),]
clin=clin[clin$bcr_patient_barcode%in%samples,]

#Indel neoantigen samples
indelneoantigen<-fread("~/Box Sync/Others/controledData/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.indel_neoantigens_10062017.tsv",data.table = F,stringsAsFactors = F)

##########DDR mutation#############
mmrgene=c("MLH1","MLH3","MSH2","MSH3","MSH6","PMS1","PMS2","BRCA1","BRCA2","PALB2","POLE","POLQ","ATM","ATR","CHEK2")
gddrspl=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%mmrgene])
sddrspl=unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%mmrgene])
DDRmut=unique(c(gddrspl,sddrspl))

##########Neoantigen hotspots################
indelhot=read.table(sep="\t",header=T,file="/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/out/mutation_Indel_neoantigen_association_byCancer_annotated.txt", stringsAsFactors=FALSE)
hotspot=indelhot$Indel[which(indelhot$P<0.05)]
hotspotsample=unique(indelneoantigen$sample[which(indelneoantigen$indel%in%hotspot)])

########Neoantigen level in germline and somatic mutated sample#######
immuneProfile=read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]
immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]

immuneProfile$DDRmut=ifelse(immuneProfile$bcr_patient_barcode%in%DDRmut,"mut","none")
immuneProfile$hotspot=ifelse(immuneProfile$bcr_patient_barcode%in%hotspotsample,"hotspot","none")

immuneProfile$type=paste0(immuneProfile$DDRmut,"_",immuneProfile$hotspot)


##########################################################
####immune signature hotspots vs none (by cancer type)####
##########################################################

immuneNeo<-function(label,cancer,pheno,neo,short_name,ylim,col_value,genelist,log=TRUE,th,w,h,comparisions){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]
 
  if(any(is.na(subimmnue$value))){
    subimmnue=subimmnue[-which(is.na(subimmnue$value)),]
  }
  
  subimmnue$value=as.numeric(as.matrix(subimmnue$value))
  subimmnue$label=paste0(subimmnue$cancer,subimmnue$hotspot)
  if(log){
    subimmnue$value=log2(as.numeric(as.matrix(subimmnue$value))+0.1)
  }else{
    subimmnue$value=as.numeric(as.matrix(subimmnue$value))
  }
  
      ann_text=NULL
      pvalueMat=NULL
      for(ca in as.character(unique(subimmnue$TCGA_Study))){
        subMat=subimmnue[subimmnue$TCGA_Study==ca,]
                tmpPvalue=wilcox.test(subMat$value[subMat$hotspot=="hotspot"],subMat$value[subMat$hotspot=="none"])$p.value
                if(tmpPvalue<0.05){
                  plab="*"
                }
                if(tmpPvalue<0.01){
                  plab="**"
                }
                if(tmpPvalue<0.001){
                  plab="***"
                }
                if(tmpPvalue>0.05){
                  plab=""
                }
                
              ann_text=rbind(ann_text,t(cbind(c(cancer=ca,mut="hotspot",pheno=pheno,lab=plab))))
              pvalueMat=rbind(pvalueMat,t(cbind(c(cancer=ca,mut="hotspot",pheno=pheno,lab=as.numeric(tmpPvalue)))))
      }
      
      ann_text=as.data.frame(ann_text)
      pvalueMat=as.data.frame(pvalueMat)
      pvalueMat$lab=as.numeric(as.matrix(pvalueMat$lab))
      
      ann_text$mut=as.character(ann_text$mut)
      ann_text$cancer=as.character(ann_text$cancer)
      #ann_text$type=paste0(ann_text$cancer,"_",ann_text$mut)
 
      #within facet ordering
      subimmnue$TCGA_Study=as.character(subimmnue$TCGA_Study)
      subimmnue$hotspot=factor(as.character(subimmnue$hotspot),levels=c("hotspot","none"))
      
      sigCancer=intersect(unique(ann_text$cancer[ann_text$lab!=""]),cancer)
      subimmnue=subimmnue[subimmnue$TCGA_Study%in%sigCancer,]
      
      subimmnue$TCGA_Study=as.character(subimmnue$TCGA_Study)
      ann_text=ann_text[ann_text$cancer%in%subimmnue$TCGA_Study,]
      if(dim(subimmnue)[1]==0){next}
      
      eval(parse(text=paste0("p=ggplot(subimmnue, aes(x=hotspot, y=as.numeric(as.matrix(value))))+geom_violin()+geom_boxplot(width=0.1,color=\"grey\",outlier.size=0.2)")))
      p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
      p=p+geom_signif(comparisons = list(c("hotspot", "none")),step_increase = 0.1)
      p=p+facet_grid(.~TCGA_Study, space="free",scale="free",drop=T)
      p=p+xlab("")+ylim(ylim)+ylab(short_name)
      p=p+scale_fill_manual(values=col_value)+scale_color_manual(values=col_value)
      #eval(parse(text=paste0("p=p+geom_text(data = ann_text,aes(x=mut,y=labpos,label=ann_text$lab),colour=\"#56B4E9\")")))
      p=p+th
      print(p)
      fn = paste0("out/",label,"_",pheno,"_distribution.pdf")
      ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
 }
  
col_value=c("hotspots"="#56B4E9","none"="#56B4E9")
#PDCD1
short_name="PDCD1"
cancer=c("UCEC")#unique(clin$type)
pheno="PDCD1"

immuneNeo("immunePlusNeoHotspot_",cancer=cancer,pheno=pheno,short_name="PDCD1",ylim=c(0,13),log=FALSE,th=th,col_value=col_value,w=1.5,h=3)


#CYTScore
short_name="CYTScore"
cancer=c("UCEC")
#cancer=unique(clin$type)
pheno="CYTScore"
labpos=5

immuneNeo("immunePlusNeoCate_",cancer=cancer,pheno=pheno,short_name="CYTScore",ylim=c(0,15),log=FALSE,th=th,col_value=col_value,w=1.5,h=3)


#####################################################
####immune signature hotspots vs none (pancancer)####
#####################################################
indelhot=read.table(sep="\t",header=T,file="/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/out/DDRsignificant_enriched_Indel_neoantigens_rmhypermutator_pancancerLevel.txt", stringsAsFactors=FALSE)
#hotspot=indelhot$Indel[which(indelhot$P<0.05)]
#hotspot="17_56435160_56435161_AC_A"
hotspot=unique(indelhot$Indel)
hotspotsample=unique(indelneoantigen$sample[which(indelneoantigen$indel%in%hotspot)])

#immune score
immuneProfile=read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]
immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]

immuneProfile$DDRmut=ifelse(immuneProfile$bcr_patient_barcode%in%DDRmut,"mut","none")
immuneProfile$hotspot=ifelse(immuneProfile$bcr_patient_barcode%in%hotspotsample,"hotspot","none")

immuneProfile$type=paste0(immuneProfile$DDRmut,"_",immuneProfile$hotspot)


immuneNeoAdj<-function(label,cancer,pheno,neo,short_name,ylim,col_value,genelist,log=TRUE,th,w,h,comparisions){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]

  #adjust immune
  medarray=tapply(subimmnue$value,subimmnue$TCGA_Study,median)
  subimmnue$valueadj=apply(subimmnue,1,function(x)as.numeric(as.matrix(x["value"]))/as.numeric(as.matrix(medarray[as.character(subimmnue[1,]["TCGA_Study"])])))
  
  if(any(is.na(subimmnue$value))){
    subimmnue=subimmnue[-which(is.na(subimmnue$value)),]
  }
  
  subimmnue$valueadj=as.numeric(as.matrix(subimmnue$valueadj))
  subimmnue$label=paste0(subimmnue$cancer,subimmnue$hotspot)
  if(log){
    subimmnue$value=log2(as.numeric(as.matrix(subimmnue$value))+0.1)
  }else{
    subimmnue$value=as.numeric(as.matrix(subimmnue$value))
  }
  

  #within facet ordering
  subimmnue$hotspot=factor(as.character(subimmnue$hotspot),levels=c("hotspot","none"))
  if(dim(subimmnue)[1]==0){next}
  
  require(ggsignif)
  eval(parse(text=paste0("p=ggplot(subimmnue, aes(x=hotspot, y=as.numeric(as.matrix(valueadj))))+geom_violin()+geom_boxplot(width=0.1,color=\"grey\",outlier.size=0.2)")))
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+geom_signif(comparisons = list(c("hotspot", "none")),step_increase = 0.1)
  #p=p+facet_grid(.~TCGA_Study, space="free",scale="free",drop=T)
  p=p+xlab("")+ylim(ylim)+ylab(short_name)
  p=p+scale_fill_manual(values=col_value)+scale_color_manual(values=col_value)
  p=p+th
  print(p)
  fn = paste0("out/",label,"_",pheno,"_distribution.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}


col_value=c("hotspots"="#56B4E9","none"="#56B4E9")
short_name="TIL"
cancer=unique(clin$type)
pheno="Lymphocyte_Infiltration_Signature_Score"

immuneNeoAdj("PanCancerimmunePlusNeoHotspot_",cancer=cancer,pheno=pheno,short_name="Adjusted TIL",ylim=c(-4,4),log=FALSE,th=th,col_value=col_value,w=1.5,h=3)


cancer=unique(clin$type)
pheno="PDCD1"
immuneNeoAdj("PanCancerimmunePlusNeoHotspot_",cancer=cancer,pheno=pheno,short_name="Adjusted PD1",ylim=c(0,5),log=FALSE,th=th,col_value=col_value,w=1.5,h=3)

cancer=unique(clin$type)
pheno="CD274"
immuneNeoAdj("PanCancerimmunePlusNeoHotspot_",cancer=cancer,pheno=pheno,short_name="Adjusted PD-L1",ylim=c(-5,15),log=FALSE,th=th,col_value=col_value,w=1.5,h=3)



cancer=unique(clin$type)
pheno="CYTScore"
immuneNeoAdj("PanCancerimmunePlusNeoHotspot_",cancer=cancer,pheno=pheno,short_name="Adjusted CYTScore",ylim=c(0,3),log=FALSE,th=th,col_value=col_value,w=1.5,h=3)


#############################################################
####immune signature specific variant vs none (by cancer)####
#############################################################
hs=c("DOCK3|p.P1852Qfs*78","INPPL1|p.R1156Gfs*46","RNF43|p.G659Vfs*41")
indelhot=read.table(sep="\t",header=T,file="/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/out/mutation_Indel_neoantigen_association_annotated.txt", stringsAsFactors=FALSE)
indelhot$hs=paste0(indelhot$gene.1,"|",indelhot$aa)


hs1=hs[3]

hotspot=unique(indelhot$Indel[which(indelhot$hs==hs1)])
hotspotsample=unique(indelneoantigen$sample[which(indelneoantigen$indel%in%hotspot)])

#immune score
immuneProfile=read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]
immuneProfile$hotspot=ifelse(immuneProfile$bcr_patient_barcode%in%hotspotsample,"hotspot","none")

th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black", size=12, angle = 60,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(size=8,hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))


immuneHotspot<-function(label,cancer,hs1,pheno,neo,short_name,ylim,col_value,genelist,log=TRUE,th,w,h,comparisions){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]
  
  if(any(is.na(subimmnue$value))){
    subimmnue=subimmnue[-which(is.na(subimmnue$value)),]
  }
  
  subimmnue$value=as.numeric(as.matrix(subimmnue$value))
  subimmnue$label=paste0(subimmnue$cancer,subimmnue$hotspot)
  if(log){
    subimmnue$value=log2(as.numeric(as.matrix(subimmnue$value))+0.1)
  }else{
    subimmnue$value=as.numeric(as.matrix(subimmnue$value))
  }
  
#  ann_text=NULL
#  pvalueMat=NULL
#  for(ca in as.character(unique(subimmnue$TCGA_Study))){
#    subMat=subimmnue[subimmnue$TCGA_Study==ca,]
#    tmpPvalue=wilcox.test(subMat$value[subMat$hotspot=='hotspot'],subMat$value[subMat$hotspot=="none"])$p.value
#    if(tmpPvalue<0.05){
#      plab="*"
#    }
#    if(tmpPvalue<0.01){
#      plab="**"
#    }
#    if(tmpPvalue<0.001){
#      plab="***"
#    }
#    if(tmpPvalue>0.05){
#      plab=""
#    }
    
#    ann_text=rbind(ann_text,t(cbind(c(cancer=ca,mut='hotspot',pheno=pheno,lab=plab))))
#    pvalueMat=rbind(pvalueMat,t(cbind(c(cancer=ca,mut='hotspot',pheno=pheno,lab=as.numeric(tmpPvalue)))))
#  }
  
#  ann_text=as.data.frame(ann_text)
#  pvalueMat=as.data.frame(pvalueMat)
#  pvalueMat$lab=as.numeric(as.matrix(pvalueMat$lab))
  
#  ann_text$mut=as.character(ann_text$mut)
#  ann_text$cancer=as.character(ann_text$cancer)
  #ann_text$type=paste0(ann_text$cancer,"_",ann_text$mut)
  
  #within facet ordering
#  subimmnue$TCGA_Study=as.character(subimmnue$TCGA_Study)
#  subimmnue$hotspot=factor(as.character(subimmnue$hotspot),levels=c("hotspot","none"))
  
#  sigCancer=intersect(unique(ann_text$cancer[ann_text$lab!=""]),cancer)
#  subimmnue=subimmnue[subimmnue$TCGA_Study%in%sigCancer,]
  
#  subimmnue$TCGA_Study=as.character(subimmnue$TCGA_Study)
#  ann_text=ann_text[ann_text$cancer%in%subimmnue$TCGA_Study,]
#  if(dim(subimmnue)[1]==0){next}
  
  eval(parse(text=paste0("p=ggplot(subimmnue, aes(x=hotspot, y=as.numeric(as.matrix(value))))+geom_violin()+geom_boxplot(width=0.1,color=\"grey\",outlier.size=0.2)")))
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+geom_signif(comparisons = list(c("hotspot", "none")),step_increase = 0.1)
  p=p+facet_grid(.~TCGA_Study, space="free",scale="free",drop=T)
  p=p+xlab("")+ylim(ylim)+ylab(short_name)+ggtitle(hs1)
  p=p+scale_color_manual(values=col_value)#+scale_fill_manual(values=col_lab)
  #eval(parse(text=paste0("p=p+geom_text(data = ann_text,aes(x=mut,y=labpos,label=ann_text$lab),colour=\"#56B4E9\")")))
  p=p+th
  print(p)
  fn = paste0("./out/",cancer,"_",hs1,"_",label,"_",pheno,"_distribution.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}

#eval(parse(text=paste0("col_value=c(\"",hs1,"\"=\"#56B4E9\",\"none\"=\"#56B4E9\")")))
col_value=c("hotspot"="#56B4E9","none"="#56B4E9")
cancer="COADREAD"
pheno="PDCD1"
immuneHotspot("SpecificHotspot_",hs1=hs1,cancer=cancer,pheno=pheno,short_name="PDCD1",ylim=c(0,15),log=FALSE,th=th,col_value=col_value,w=2,h=3)

pheno="Lymphocyte_Infiltration_Signature_Score"
immuneHotspot("SpecificHotspot_",hs1=hs1,cancer=cancer,pheno=pheno,short_name="TIL",ylim=c(0,15),log=FALSE,th=th,col_value=col_value,w=2,h=3)

pheno="CYTScore"
immuneHotspot("SpecificHotspot_",hs1=hs1,cancer=cancer,pheno=pheno,short_name="CYTScore",ylim=c(0,15),log=FALSE,th=th,col_value=col_value,w=2,h=3)

#UCEC
col_value=c("hotspot"="#56B4E9","none"="#56B4E9")
cancer="UCEC"
pheno="PDCD1"
immuneHotspot("SpecificHotspot_",hs1=hs1,cancer=cancer,pheno=pheno,short_name="PDCD1",ylim=c(0,15),log=FALSE,th=th,col_value=col_value,w=2,h=3)

pheno="Lymphocyte_Infiltration_Signature_Score"
immuneHotspot("SpecificHotspot_",hs1=hs1,cancer=cancer,pheno=pheno,short_name="TIL",ylim=c(0,5),log=FALSE,th=th,col_value=col_value,w=2,h=3)

pheno="CYTScore"
immuneHotspot("SpecificHotspot_",hs1=hs1,cancer=cancer,pheno=pheno,short_name="CYTScore",ylim=c(0,15),log=FALSE,th=th,col_value=col_value,w=2,h=3)

#############################################################
####immune signature specific variant vs none (pancancer)####
#############################################################
immuneHotspot<-function(label,cancer,hs,pheno,neo,short_name,ylim,col_value,genelist,log=TRUE,th,w,h,comparisions){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]
  
  #adjust immune
  medarray=tapply(subimmnue$value,subimmnue$TCGA_Study,median)
  subimmnue$valueadj=apply(subimmnue,1,function(x)as.numeric(as.matrix(x["value"]))/as.numeric(as.matrix(medarray[as.character(subimmnue[1,]["TCGA_Study"])])))
  
  if(any(is.na(subimmnue$value))){
    subimmnue=subimmnue[-which(is.na(subimmnue$value)),]
  }
  
  subimmnue$valueadj=as.numeric(as.matrix(subimmnue$valueadj))
  subimmnue$label=paste0(subimmnue$cancer,subimmnue$hotspot)
  if(log){
    subimmnue$value=log2(as.numeric(as.matrix(subimmnue$value))+0.1)
  }else{
    subimmnue$value=as.numeric(as.matrix(subimmnue$value))
  }
  
  
  #within facet ordering
  subimmnue$hotspot=factor(as.character(subimmnue$hotspot),levels=c("hotspot","none"))
  if(dim(subimmnue)[1]==0){next}
  
  require(ggsignif)
  eval(parse(text=paste0("p=ggplot(subimmnue, aes(x=hotspot, y=as.numeric(as.matrix(valueadj))))+geom_violin()+geom_boxplot(width=0.1,color=\"grey\",outlier.size=0.2)")))
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+geom_signif(comparisons = list(c("hotspot", "none")),step_increase = 0.1)
  #p=p+facet_grid(.~TCGA_Study, space="free",scale="free",drop=T)
  p=p+xlab("")+ylim(ylim)+ylab(short_name)
  p=p+scale_fill_manual(values=c("hotspot"=hs))+scale_color_manual(values=col_value)
  p=p+th
  print(p)
  fn = paste0("./out/",label,"_",pheno,"_distribution.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}


