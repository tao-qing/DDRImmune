#######################
####Part1 Load data####
#######################
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")
library(data.table)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)

#short name
xlabel=scale_y_discrete(position = "left",labels=c("Aneuploidy_Score" = "Aneuploidy", "BCR_Evenness" = "BCR_E","BCR_Richness" = "BCR_R","BCR_Shannon"="BCR_S","CTA_Score"="CTA","Homologous_Recombination_Defects"="HRD","TIL_Regional_Fraction"="TIL","Leukocyte_Fraction"="Leukocyte","Macrophage_Regulation"="Macrophage","Lymphocyte_Infiltration_Signature_Score"="Lymphocyte","IFN_gamma_Response"="IFN","TGF_beta_Response"="TGF_beta","Th1_Cells"="Th1","Th2_Cells"="Th2","Th17_Cells"="Th17"))

th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black", size=12, angle = 60,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

##########Neoantigen level in germline and somatic mutated sample######### 
source("../dependency_files_tq.R")
source("../load_somatic_rmhypermut.R")
rmhypermutator=TRUE
label="_rmhypermutator_"

#immune score
immuneProfile=read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]

#cancer samples
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]

  immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]
  pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%samples),]
  somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$bcr_patient_barcode%in%samples,]

##########################################################
####Figure1 germline somatic genes to neoantigen level####
##########################################################
immuneFeature<-function(label,cancer,pheno,short_name,labpos,col_value,genelist,log=TRUE,ord,cancerord=cancerord,th,w,h,comparisions){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  cancers=clin$type[pmatch(samples,clin$bcr_patient_barcode)]
  mut_status=sapply(samples,function(x){
    g=pathVarP$HUGO_Symbol[pathVarP$bcr_patient_barcode==x]
    s=somatic_likelyfunctional_driver$Hugo_Symbol[somatic_likelyfunctional_driver$bcr_patient_barcode==x]
    gindex=which(unlist(lapply(genelist,function(x)any(x%in%g))))
    sindex=which(unlist(lapply(genelist,function(x)any(x%in%s))))
    label=paste0(paste0(ifelse(length(gindex)==0,"",paste0("g",names(genelist)[gindex])),collapse = "-"),paste0(ifelse(length(sindex)==0,"",paste0("s",names(genelist)[sindex])),collapse = "-"),collapse = "-")
    return(label)
  })
  
  mutMat=cbind(samples,pheno=subimmnue$value[pmatch(samples,subimmnue$bcr_patient_barcode)],TCGA_Study=cancers,mut=mut_status) %>% as.data.frame() %>% mutate(mut=replace(as.character(mut),mut==" ","")) %>% mutate(mut=replace(as.character(mut),mut=="","WT"))%>%mutate(pheno=as.numeric(as.matrix(pheno)))%>%filter(!is.na(mut))%>%filter(!is.na(pheno))%>%mutate(label=paste0(TCGA_Study,"_",mut))%>%group_by(label)%>%filter(n() > 3)
  
  if(log){
    mutMat$pheno=log2(as.numeric(as.matrix(mutMat$pheno))+0.1)
  }

  neo_comps <- compare_means(pheno ~ mut,data=mutMat, group.by = "TCGA_Study")%>%filter(group1=="WT")%>%mutate(p.adj=p.adjust(p,method = "fdr")) %>%mutate(p.adj.signif=ifelse(p.adj<0.001,"***",ifelse(p.adj<0.01,"**",ifelse(p.adj<0.05,"*",""))))%>%filter(TCGA_Study%in% TCGA_Study[which(p.adj.signif!="")])%>%mutate(label=as.character(paste0(TCGA_Study,"_",group2)))%>%mutate(TCGA_Study=factor(TCGA_Study,levels=cancerord))
  
  plotMat<-mutMat%>%filter(label%in%neo_comps$label | (mut=="WT" & TCGA_Study%in%as.character(neo_comps$TCGA_Study)))%>%mutate(TCGA_Study=factor(TCGA_Study,levels=cancerord))%>%mutate(mut=factor(mut,levels=ord))
  
  
  p=ggplot(plotMat, aes(x=mut, y=as.numeric(as.matrix(pheno))))+geom_violin()+geom_boxplot(width=0.1,color="grey",outlier.size=0.2)
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+facet_grid(.~TCGA_Study, space="free",scale="free",drop=T)
  p=p+xlab("")+ylab(short_name)
  #p=p+scale_fill_manual(values=col_value)+scale_color_manual(values=col_value)
  p=p+geom_text(data = neo_comps,aes(x=group2,y=labpos,label=p.adj.signif),colour="#56B4E9")
  p=p+th
  print(p)
  fn = paste0("out/",label,"_",pheno,"_distribution_mutated_BRCAMMR_samples.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
  
  return(neo_comps)
}

genelist=NULL
genelist["HR"]=list(c("BRCA1","BRCA2","PALB2"))
genelist["Polymerase"]=list(c("POLE","POLQ"))
genelist["Sensor"]=list(c("ATM","ATR","CHEK2"))
genelist["MMR"]=list(c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2"))

comp=list(c("sMMR","WT"),c("sHR","WT"),c("gHR","WT"),c("gMMR","WT"))
ord=c("sMMR","sHR","gHR","sPolymerase","gPolymerase","sSensor","gSensor","WT")
cancerord=c("UCEC","COADREAD","CESC","BRCA","BLCA")

#Neoantigene level
col_value=c("sMMR"="#56B4E9","gMMR"="#56B4E9","sSensor"="#56B4E9","gSensor"="#56B4E9","sPolymerase"="#56B4E9","gPolymerase"="#56B4E9","gHR"="#56B4E9","sHR"="#56B4E9","WT"="#dbd9d9")
cancer=unique(clin$type)#c("BRCA","OV","UCEC","COADREAD","BLCA")#c(unique(as.character(contributionNew$ca)))
pheno="SNV_Neoantigens"
snvneo=immuneFeature("Neodist_",cancer,pheno,short_name="log2(SNV Neoantigens)",labpos=15,genelist=genelist,log=TRUE,ord=ord,cancerord=cancerord,th=th,col_value=col_value,w=7,h=3,comparisions = comp)


cancer=unique(clin$type)
pheno="Indel_Neoantigens"
ord=c("sMMR","gMMR","sHR","gHR","sPolymerase","gPolymerase","sSensor","gSensor","WT")
cancerord=unique(clin$type)
col_value=c("sMMR"="#56B4E9","gMMR"="#56B4E9","sPolymerase"="#56B4E9","gPolymerase"="#56B4E9","sSensor"="#56B4E9","gSensor"="#56B4E9","gHR"="#56B4E9","sHR"="#56B4E9","WT"="#dbd9d9")
indelneo=immuneFeature("Neodist_",cancer,pheno,short_name="log2(Indel Neoantigens)",labpos=15,genelist=genelist,ord=ord,cancerord=cancerord,col_value=col_value,log=TRUE,th=th,w=4.5,h=3)

merge=rbind(cbind(Type="SNVneo",snvneo),cbind(Type="Indelneo",indelneo))
write.csv(merge,"Figure3AB_indelneo_test.csv",quote=F)


