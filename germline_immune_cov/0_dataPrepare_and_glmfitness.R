setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/")
source("../global_aes_out.R")
source("../dependency_files_tq.R")


fn = "../../Huang_lab_data/TCGA_PanCanAtlas_2018/Immune_Thorsson_Immunity2018/1-s2.0-S1074761318301213-mmc2.xlsx"
immune_profile = data.frame(readxl::read_xlsx(fn))
immune_profile$TCGA.Study=gsub("^COAD$","COADREAD",gsub("^READ$","COADREAD",immune_profile$TCGA.Study,perl = T),perl = T)

#only keep Eosinophils...41 and Neutrophils...48
#plot(immune_profile$'Eosinophils...41',immune_profile$'Eosinophils...61')
#plot(immune_profile$'Neutrophils...48',immune_profile$'Neutrophils...60')
immune_profile=immune_profile[,-which(colnames(immune_profile)%in%c("Eosinophils...61","Neutrophils...60"))]
colnames(immune_profile)=gsub("Eosinophils...41","Eosinophils",colnames(immune_profile))
colnames(immune_profile)=gsub("Neutrophils...48","Neutrophils",colnames(immune_profile))

colnames(immune_profile) = gsub("\\.","_",colnames(immune_profile))
colnames(immune_profile)[1] = "bcr_patient_barcode"

#distribution of Homologous_Recombination_Defects
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(Homologous_Recombination_Defects))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("HRD") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_HRD.pdf")
#ggsave(fn,height=3,width=3)

#distribution of IFN_gamma_Response
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(IFN_gamma_Response))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("IFN") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_IFN.pdf")
#ggsave(fn,height=3,width=3)


#distribution of TGF_beta_Response
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(TGF_beta_Response))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("TGF_beta_Response") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_TGF.pdf")
#ggsave(fn,height=3,width=3)

#distribution of Macrophage_Regulation
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(Macrophage_Regulation))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("Macrophage_Regulation") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_Macrophage.pdf")
#ggsave(fn,height=3,width=3)

#distribution of Proliferation
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(Proliferation))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("Proliferation") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_Proliferation.pdf")
#ggsave(fn,height=3,width=3)

#distribution of Lymphocyte_Infiltration_Signature_Score
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(Lymphocyte_Infiltration_Signature_Score))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("Lymphocyte_Infiltration") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_Lymphocyte_Infiltration.pdf")
#ggsave(fn,height=3,width=3)


#distribution of Wound_Healing
#p=ggplot(immune_profile, aes(as.numeric(as.matrix(Wound_Healing))))+   geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("Wound_Healing") +
#  ylab("Density") +
#  theme(legend.title=element_blank())#+scale_x_sqrt()
#fn=paste0("./out/distribution_of_Wound_Healing.pdf")
#ggsave(fn,height=3,width=3)

#according to table 1: /Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/Table1. immune_feature.xlsx
####################Genomic biomarkers of immunotherapy#######################
#somatic mutation burden
somatic_f = "../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
somatic<-somatic%>%mutate(bcr_patient_barcode=substr(Tumor_Sample_Barcode,1,12))
somaticExome=somatic[somatic$Variant_Classification!="Intron",]

#exome capture regions
exome=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/gaf_20111020Plusbroad_wex_1.1_hg19.bed",h=F)
exomelength=sum(apply(exome,1,function(x)as.numeric(as.matrix(x[3]))-as.numeric(as.matrix(x[2]))))/1000000

TMB=as.data.frame(table(somaticExome$bcr_patient_barcode)/exomelength)
TMB_SNV=as.data.frame(table(somaticExome$bcr_patient_barcode[-which(somatic$Variant_Classification%in%c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"))])/exomelength)
TMB_INDEL=as.data.frame(table(somaticExome$bcr_patient_barcode[which(somatic$Variant_Classification%in%c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"))])/exomelength)

#update 20191117
#tmbmat=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/Cell-of-Origin_Cell2018_mutation-load_updated.txt",h=T,sep="\t")
#tmbmat$TMB=tmbmat$Silent.per.Mb+tmbmat$Non.silent.per.Mb
#distribution of mutation burden
#p=ggplot(burden, aes(Freq)) +
#  geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("Somatic Mutation Burden") +
#  ylab("Density") +
#  theme(legend.title=element_blank())+scale_x_sqrt()
#fn=paste0("./out/distribution_of_TMB.pdf")
#ggsave(fn,height=3,width=3)


#MSIsensor score
msi<-read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/genomeInstability/Ding_Cell_2018_MSIsensor_Score.txt",sep="\t",h=T)
msicate=as.data.frame(cbind(Participant.Barcode=as.character(msi$Participant.Barcode),genes=as.character(msi$Samples.with.non.silent.mutations.in.MSI.genes),MSICate=ifelse(msi$MSIsensor.score>=4,1,0)))
#distribution of MSIsensor
#p=ggplot(msi, aes(MSIsensor.score)) +
#  geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("MSIsensor Score") +
#  ylab("Density") +
#  theme(legend.title=element_blank())+scale_x_sqrt()
#fn=paste0("./out/distribution_of_MSIsensor.pdf")
#ggsave(fn,height=3,width=3)

#using an MSIsensor score ≥ 10 to define MSI high (MSI-H)
#msistatus=ifelse(msi$MSIsensor.score>=10,1,0)
#names(msistatus)=msi$Participant.Barcode

####################DNA-repair mutational and expression signature#######################
#PARPi7 expression
PARPi7<-read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/PARPi7_expression.txt",sep="\t",h=T)
#distribution of MSIsensor
#p=ggplot(PARPi7, aes(PARPi7)) +
#  geom_density() +
#  theme_bw(base_size = 10) +
#  theme(legend.position ="right") +
#  xlab("PARPi7") +
#  ylab("Density") +
#  theme(legend.title=element_blank())+scale_x_sqrt()
#fn=paste0("./out/distribution_of_PARPi7.pdf")
#ggsave(fn,height=3,width=3)

#wnt signaling signature
raw=fread("/Users/qingtao/Box Sync/GermlineSomatic/analysis/RNF43/TCGA_KEGG_ssGSEA_score.txt-combined.gct",stringsAsFactors = F,h=T,skip = 2,data.table=F)

score=t(raw[,grep("^TCGA|id",colnames(raw),perl=T)])
colnames(score)=score[1,]
score=score[-1,]
score=as.data.frame(score)
score$bcr_patient_barcode=rownames(score)



immune_profile=immune_profile%>%mutate(TMB=TMB$Freq[pmatch(immune_profile$bcr_patient_barcode,TMB$Var1)],TMB_INDEL=TMB_INDEL$Freq[pmatch(immune_profile$bcr_patient_barcode,TMB_INDEL$Var1)],TMB_SNV=TMB_SNV$Freq[pmatch(immune_profile$bcr_patient_barcode,TMB_SNV$Var1)],MSISensor=msi$MSIsensor.score[pmatch(immune_profile$bcr_patient_barcode,msi$Participant.Barcode)],MSICate=msicate$MSICate[pmatch(immune_profile$bcr_patient_barcode,msicate$Participant.Barcode)])%>%mutate(PARPi7=PARPi7$PARPi7[pmatch(immune_profile$bcr_patient_barcode,PARPi7$patient_barcode)])%>%mutate(Wnt=score$WNT_SIGNALING_PATHWAY[pmatch(immune_profile$bcr_patient_barcode,score$bcr_patient_barcode)])

#update 20191117
#immune_profile=immune_profile%>%mutate(TMB=tmbmat$TMB[pmatch(immune_profile$bcr_patient_barcode,tmbmat$Patient_ID)],MSISensor=msi$MSIsensor.score[pmatch(immune_profile$bcr_patient_barcode,msi$Participant.Barcode)])%>%mutate(PARPi7=PARPi7$PARPi7[pmatch(immune_profile$bcr_patient_barcode,PARPi7$patient_barcode)])

####################Immune infiltrate composition#######################

#xCell
xCell<-read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/xCell_Aran_GenomeBio2017/xCell_TCGA_RSEM.txt",sep="\t",h=F)
xCell=t(xCell)
colnames(xCell)=paste0("xCell_",gsub(" ","_",xCell[1,]))
xCell=xCell[-1,]
colnames(xCell)[1]="bcr_patient_barcode"
xCell[,1]=substr(gsub("\\.","-",xCell[,1]),1,12)
xCell=as.data.frame(xCell)

#write.table(colnames(xCell),"xCell",quote=F,row.names = F,col.names = F)

overlap=intersect(immune_profile$bcr_patient_barcode,xCell$bcr_patient_barcode)
immune_profile=cbind(immune_profile,xCell[pmatch(immune_profile$bcr_patient_barcode,xCell$bcr_patient_barcode),-1])

#chech whether the data frames are correctly merged
#j=9500
#immune_profile[j,"xCell_Macrophages_M2"]
#xCell[which(xCell$bcr_patient_barcode==immune_profile$bcr_patient_barcode[j]),"xCell_Macrophages_M2"]

#for(i in grep("xCell",colnames(immune_profile),value=T)){
#  cybersortName=gsub("xCell_","",i)
#  if(any(colnames(immune_profile)==cybersortName)){
    #png(paste0("./out/",cybersortName,".png"),height=800,width=800,res=120)
    #plot(as.numeric(as.matrix(immune_profile[,i])),as.numeric(as.matrix(immune_profile[,cybersortName])),xlab="xCell",ylab="CyberSort",main=cybersortName)
    #dev.off()
#    plotMat=as.data.frame(cbind(Xcell=as.numeric(as.matrix(immune_profile[,i])),cybersort=as.numeric(as.matrix(immune_profile[,cybersortName]))))
#    plotMat$Xcell=as.numeric(as.matrix(plotMat$Xcell))
#    plotMat$cybersort=as.numeric(as.matrix(plotMat$cybersort))
#    subplotMat=plotMat[!(is.na(plotMat$Xcell) | is.na(plotMat$cybersort)),]
#    cor=round(cor( subplotMat$Xcell, subplotMat$cybersort),digits=2)
#    library(ggrepel)
#    p = ggplot(plotMat,aes(y=cybersort, x =Xcell))
#    p = p + geom_point(stroke=0,alpha = 0.6,size=1)
#    p = p + geom_smooth(data=plotMat,aes(y=cybersort, x =Xcell),colour="blue",method = "glm",size=0.2)
#    p = p + ylab("Cybersort") + xlab("Xcell") + theme_bw() + ggtitle(paste0(cybersortName,"(cor=",cor,")"))
#    p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) 
#    p = p + theme(legend.position = "none",plot.title = element_text(size=10, face="bold"),legend.key=element_blank(),axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))+ guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA)))
    #p = p + expand_limits(x = 0,y=0) + scale_y_sqrt(limits = c(0,axis_mat1))+ scale_x_sqrt(limits = c(0,axis_mat1))
#    fn = paste0("./out/",cybersortName,".pdf")
#    ggsave(fn,h=3,w=3,useDingbat=F)
#  }
#}

##################Immune Expression#########################
igMat<-fread("./out/Immunity_Thorsson_2018_Pardoll_2012_ImmuneCheckpoint_Exp.tsv",data.table = F)
colnames(igMat)[1]="bcr_patient_barcode"
#igMat=igMat[,-1]

immune_profile=cbind(immune_profile,igMat[pmatch(immune_profile$bcr_patient_barcode,igMat$bcr_patient_barcode),-1])

#immunotherapy expression signature
#CYT Score expression,GZMA, PRF1 "https://science.sciencemag.org/content/suppl/2019/05/01/364.6439.485.DC1; https://www.sciencedirect.com/science/article/pii/S0092867414016390" The cytolytic activity of the immune infiltrate (CYT score) was defined based on the transcript levels of two key cytolytic effectors, namely GZMA (granzyme A) and PRF1 (perforin), which are upregulated upon CD8+ T cell activation (25). For this, RNASeqV2 data were obtained from the Broad’s Institute Firehose (01/28/2016). Following [https://gdac.broadinstitute.org], we took, per tumor, the geometric mean of the expression in transcripts (‘scaled estimate’ values in the expression data, representing TPM as calculated by RSEM) in the two genes specified above"

CYTScore= apply(igMat[,c("GZMA","PRF1")],1,mean)
names(CYTScore)=igMat[,1]

immune_profile=cbind(immune_profile,CYTScore=CYTScore[pmatch(immune_profile$bcr_patient_barcode,names(CYTScore))])

#CTL signature 	CD8A, CD8B, GZMA, GZMB and PRF1 	https://www.nature.com/articles/s41591-018-0136-1	we used the average expression level of CD8A, CD8B, GZMA, GZMB and PRF1 to estimate the cytotoxic T lymphocyte (CTL) level in a tumor
CTLScore= apply(igMat[,c("IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA")],1,mean)
names(CTLScore)=igMat[,1]
immune_profile=cbind(immune_profile,CTLScore=CTLScore[pmatch(immune_profile$bcr_patient_barcode,names(CTLScore))])


####TIL fraction#####
til=data.frame(readxl::read_xlsx("/Users/qingtao/Box Sync/GermlineSomatic/reference_files/TCGA_DeepLearningPathTILs_Saltz_CellRep2018/mmc2.xlsx"))

immune_profile=cbind(immune_profile,TIL=til$til_percentage[pmatch(immune_profile$bcr_patient_barcode,til$ParticipantBarcode)])


##################Immune Signature#############################
mut_signature_f = "../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/signature_profile_sample.txt.gz"
mut_signature = read.table(sep="\t",header=T,file=gzfile(mut_signature_f), stringsAsFactors=FALSE)
mut_signature_tcga = mut_signature[mut_signature$Country=="United States",]
mut_signature_tcga$cancer = gsub("-US","",mut_signature_tcga$project_code)
mut_signature_tcga$Signature = gsub("\\.","-",mut_signature_tcga$Signature)
mut_signature_tcga$bcr_patient_barcode = gsub("\\.","-",mut_signature_tcga$Tumor_Sample_Barcode)
mut_signature_tcga_brief = mut_signature_tcga[,c("cancer" , "Signature", "Contribution", "bcr_patient_barcode")]


#reverse reshapre 
mut_signature_tcga_all=dcast(data = mut_signature_tcga_brief,formula =bcr_patient_barcode~Signature,fun.aggregate = sum,value.var = "Contribution")

#merge signature
immune_profile=cbind(immune_profile,mut_signature_tcga_all[pmatch(immune_profile$bcr_patient_barcode,mut_signature_tcga_all$bcr_patient_barcode),-1])

#only include 10389 samples
pathVarImmune_profile=immune_profile[which(immune_profile$bcr_patient_barcode%in%clin$bcr_patient_barcode),]#remain 10260

#include covariates, age,sex,race
covariates=fread("../../Huang_lab_data/TCGA_PanCanAtlas_2018/covariates/TCGA_ancestry_PC.txt",data.table = F)

merge=cbind(pathVarImmune_profile,covariates[pmatch(pathVarImmune_profile$bcr_patient_barcode,covariates$bcr_patient_barcode),c("age_at_initial_pathologic_diagnosis","washu_assigned_ethnicity","PC1","PC2")])
colnames(merge)=gsub("age_at_initial_pathologic_diagnosis","Age",colnames(merge))
colnames(merge)=gsub("washu_assigned_ethnicity","Race",colnames(merge))

#reshape to check score distribution
library(reshape2)
mdata <- melt(merge, id=c("bcr_patient_barcode","TCGA_Study","Immune_Subtype","TCGA_Subtype","Age","Race","PC1","PC2"))

tn = "./out/pca10260_immuneprofile_covariates.txt"
write.table(mdata, quote=F, sep="\t", file = tn, row.names = F)


#overlap=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_n9967.txt",h=F)[,1]

#mdataoverlap=mdata[mdata$bcr_patient_barcode%in%overlap,]
#tn = "./out/pca9967_overlaped_immuneprofile_covariates.txt"
#write.table(mdataoverlap, quote=F, sep="\t", file = tn, row.names = F)

#overlaphyper=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9631.txt",h=F)[,1]

#mdataoverlaprmhyper=mdata[mdata$bcr_patient_barcode%in%overlaphyper,]
#tn = "./out/pca9631_overlaped_rmhyper_immuneprofile_covariates.txt"
#write.table(mdataoverlaprmhyper, quote=F, sep="\t", file = tn, row.names = F)
####################END###################

library(ggplot2)
#figure 1a distribution of deleterious mutation
p <- ggplot(mdata, aes(x=variable, y=as.numeric(as.matrix(value)),color=TCGA_Study)) + stat_summary(fun.data="mean_sdl",colour="grey",fun.args = list(mult=1), geom="crossbar", width=0.5)+ geom_dotplot(binaxis='y',binwidth=0.1, stackdir='center',alpha=0.1,dotsize=0.6)
p <- p + ylab("Immune Score")+ xlab("")+ggtitle("") + theme_bw() + theme(legend.position="none",axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=10,angle=90,hjust=1), axis.text.y = element_text(colour="black", size=10))

fn="immune_score_distribution.png"
ggsave(fn,height=10,width=15)


#check the correlation between covariats and immune variables 
cormat=NULL
for(covariat in c("Age","Race","Sex","MutLoad")){
  for(s in colnames(merge)[5:64]){
    tmp<-cbind(as.numeric(as.matrix(merge[,covariat])),as.numeric(as.matrix(merge[,s])))
    tmp<-tmp[-c(which(is.na(tmp[,1])),which(is.na(tmp[,2]))),]
    cc=cor(tmp[,1],tmp[,2])
    p=cor.test(tmp[,1],tmp[,2])$p.value
    cormat=rbind(cormat,c(covariat,s,cc,p))
  }
}

colnames(cormat)=c("covariates","immuneScore","correlation","pvalue")

tn = "./out/covariates_and_variable_correlation.txt"
write.table(cormat, quote=F, sep="\t", file = tn, row.names = F)


