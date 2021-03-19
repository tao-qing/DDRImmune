setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")
library(data.table)
mainContribution<-fread("./out/PLS-PM_germline_somatic_contribution_MMR_MSI_genes_20200204.tsv")
mainContribution1=mainContribution[grep("Neoantigen",mainContribution$Phenotype),]

subContribution<-read.table("./out/PLS-PM_genelevel_contribution_MMR_MSI_genes_20200204.tsv",h=T,stringsAsFactors = F)
subContribution1=subContribution[grep("Neoantigen",subContribution$phenotype),]


#####Figure A  germline somatic contribution#####
contribution=NULL
for(ge in c("BRCA1","BRCA2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")){
  for(ca in unique(mainContribution1$Cancer)){
    for(ph in unique(mainContribution1$Phenotype)){
      gcont=as.numeric(subContribution1[subContribution1$genes%in%paste0("G:",ge) & subContribution1$cancer==ca & subContribution1$phenotype==ph,"contribution"]) * as.numeric(mainContribution1[mainContribution1$Cancer==ca & mainContribution1$Phenotype==ph,"germline_Contribution"])
      scont=as.numeric(subContribution1[subContribution1$genes%in%paste0("S:",ge) & subContribution1$cancer==ca & subContribution1$phenotype==ph,"contribution"]) * as.numeric(mainContribution1[ mainContribution1$Cancer==ca & mainContribution1$Phenotype==ph,"somatic_Contribution"])
      if(length(gcont)==0){
        gcont=0
      }
      if(length(scont)==0){
        scont=0
      }
      contribution=rbind(contribution,cbind(ge,ca,ph,gcont,scont))
    }
  }
}

contribution=as.data.frame(contribution)
contribution$ge=as.character(contribution$ge)
contribution$ge[which(contribution$ge%in%c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2"))]="MMR"

contributionNew=melt(contribution,id=c("ge","ca","ph"))
contributionNew$Cate=apply(contributionNew,1,function(x)if(x[4]=="gcont"){paste0("g",x[1])}else{paste0("s",x[1])})
contributionNew$value=as.numeric(as.matrix(contributionNew$value))

#sort by mean
conttrivalue=tapply(contributionNew$value,contributionNew$Cate,mean)
contributionNew$Cate=factor(contributionNew$Cate,levels=names(conttrivalue)[order(conttrivalue,decreasing = T)])

library(ggplot2)
library(ggsignif)
th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black", size=12, angle = 60,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

p=ggplot(contributionNew, aes(x=Cate, y=as.numeric(as.matrix(value))))+geom_violin()+geom_boxplot(width=0.1,color="grey",outlier.size = 0.2)+ylim(0,0.4)
#p = p + facet_wrap(ph~.)
p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
p=p+ylab("Contribution to neoantigens")+xlab("")
p=p+scale_y_sqrt()
p=p+stat_signif(comparisons = list(c("sBRCA2","gBRCA1"),c("sBRCA2","sBRCA1"),c("sBRCA2","sMMR")),step_increase = 0.1,map_signif_level=F,test = "wilcox.test")
p=p+th
p
fn="./out/Figure5_germlineBRCA12_contribution2neoantigen.pdf"
ggsave(fn,useDingbat=F,width=2.5,height=3.8)


#Figure B  germline somatic contribution only mutant samples
contributionNewDel0=contributionNew[-which(contributionNew$value==0),]#remove gene do not have any contribution
conttrivaluedel0=tapply(contributionNewDel0$value,contributionNewDel0$Cate,mean)
contributionNewDel0$Cate=factor(contributionNewDel0$Cate,levels=names(conttrivaluedel0)[order(conttrivaluedel0,decreasing = T)])


p=ggplot(contributionNewDel0, aes(x=Cate, y=as.numeric(as.matrix(value))))+geom_violin()+geom_boxplot(width=0.1,color="grey",outlier.size = 0.2)+ylim(0,0.4)
#p = p + facet_wrap(ph~.)
p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
p=p+ylab("Contribution to neoantigens")+xlab("")
p=p+scale_y_sqrt()
p=p+stat_signif(comparisons = list(c("gBRCA2","sBRCA2"),c("sBRCA2","sBRCA1"),c("sBRCA2","sMMR"),c("gBRCA2","sMMR")),step_increase = 0.1,map_signif_level=F,test = "wilcox.test")
p=p+th
p
fn="./out/Figure5_germlineBRCA12_contribution2neoantigen_withdata.pdf"
ggsave(fn,useDingbat=F,width=2.5,height=3.8)


##########Neoantigen level in germline and somatic mutated sample######### 
#load data 
source("../dependency_files_tq.R")
rmhypermutator=TRUE
#immune score
immuneProfile=read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]

#cancer samples
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]

#if(rmhypermutator){
#  hyper<-read.table("../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/Driver_BaileyCell2018/HypermutatorSamples.txt",h=F,sep="\t",stringsAsFactors=FALSE)[,1]
  immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]
  pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%samples),]
  label="_rmhypermutator_"
  source("../load_somatic_rmhypermut.R")
  somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$bcr_patient_barcode%in%samples,]
#}else{
#  label="_"
#}


#plot(density(immuneProfile$value[immuneProfile$variable=="Indel_Neoantigens"]))

#short name
xlabel=scale_y_discrete(position = "left",labels=c("Aneuploidy_Score" = "Aneuploidy", "BCR_Evenness" = "BCR_E","BCR_Richness" = "BCR_R","BCR_Shannon"="BCR_S","CTA_Score"="CTA","Homologous_Recombination_Defects"="HRD","TIL_Regional_Fraction"="TIL","Leukocyte_Fraction"="Leukocyte","Macrophage_Regulation"="Macrophage","Lymphocyte_Infiltration_Signature_Score"="Lymphocyte","IFN_gamma_Response"="IFN","TGF_beta_Response"="TGF_beta","Th1_Cells"="Th1","Th2_Cells"="Th2","Th17_Cells"="Th17"))

th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black", size=12, angle = 60,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))


library(ggsignif)
immuneFeature<-function(label,cancer,pheno,short_name,labpos,col_value,genelist,log=TRUE,ord,th,w,h,comparisions){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]
  mut_status=sapply(sample,function(x){
    g=pathVarP$HUGO_Symbol[pathVarP$bcr_patient_barcode==x]
    s=somatic_likelyfunctional_driver$Hugo_Symbol[somatic_likelyfunctional_driver$bcr_patient_barcode==x]
    gindex=which(unlist(lapply(genelist,function(x)any(x%in%g))))
    sindex=which(unlist(lapply(genelist,function(x)any(x%in%s))))
    label=paste(paste(ifelse(length(gindex)==0,"",paste0("g",names(genelist)[gindex])),collapse = "-"),paste(ifelse(length(sindex)==0,"",paste0("s",names(genelist)[sindex])),collapse = "-"),collapse = "-")
    return(label)
  })
  
  mutMat=as.data.frame(cbind(sample,pheno=subimmnue$value[pmatch(sample,subimmnue$bcr_patient_barcode)],cancer=cancers,mut=mut_status)) 
  
  mutMat$mut=gsub(" ","",mutMat$mut)
  mutMat$mut=ifelse(mutMat$mut=="","NA",mutMat$mut)
  
  mutMat$label=paste0(mutMat$cancer,mutMat$mut)
  mutMat=mutMat[!mutMat$label%in%names(which(table(mutMat$label)<=3)),]
  if(log){
    mutMat$pheno=log2(as.numeric(as.matrix(mutMat$pheno))+0.1)
  }else{
    mutMat$pheno=as.numeric(as.matrix(mutMat$pheno))
  }
  
  #within facet ordering
  mutMat$mut=factor(as.character(mutMat$mut),levels=c(ord))
  
  #mutMat= mutMat%>%ungroup() %>% arrange(cancer,pheno) %>%  mutate(.r = row_number())
  if(length(which(is.na(mutMat$mut)))>0){
    mutMat=mutMat[-which(is.na(mutMat$mut)),] 
  }
  
  mutMat$pheno=as.numeric(as.matrix(mutMat$pheno))
  mutMat=mutMat[-which(is.na(mutMat$pheno)),]
  
  mutype=as.character(unique(mutMat$mut))
  mutype=mutype[-which(mutype=="NA")]

  ann_text=NULL
  pvalueMat=NULL
  for(ca in as.character(unique(mutMat$cancer))){
    subMat=mutMat[mutMat$cancer==ca,]
    for(m in mutype){
      if(!any(subMat$mut==m)){
        #plab="" 
        next
      }else{
        tmpPvalue=wilcox.test(subMat$pheno[subMat$mut==m],subMat$pheno[subMat$mut=="NA"])$p.value
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
      }
      ann_text=rbind(ann_text,t(cbind(c(cancer=ca,mut=m,pheno=pheno,lab=plab))))
      pvalueMat=rbind(pvalueMat,t(cbind(c(cancer=ca,mut=m,pheno=pheno,lab=as.numeric(tmpPvalue)))))
    }
  }
  ann_text=as.data.frame(ann_text)
  pvalueMat=as.data.frame(pvalueMat)
  pvalueMat$lab=as.numeric(as.matrix(pvalueMat$lab))
  
  ann_text$mut=as.character(ann_text$mut)
  ann_text$cancer=as.character(ann_text$cancer)
  
  mutMat$cancer=factor(mutMat$cancer,levels=cancer)
  p=ggplot(mutMat, aes(x=mut, y=as.numeric(as.matrix(pheno)),color=I(mut),fill=I(mut)))+geom_violin()+geom_boxplot(width=0.1,color="grey",outlier.size=0.2)
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+facet_grid(~cancer, space="free",scale="free",drop=T)
  p=p+xlab("")+ylab(short_name)
  p=p+scale_fill_manual(values=col_value)+scale_color_manual(values=col_value)
  p=p+geom_text(data = ann_text,aes(x=mut,y=labpos,label=ann_text$lab),colour="#56B4E9")
  p=p+th
  print(p)
  fn = paste0("out/",label,"_",pheno,"_distribution_mutated_BRCAMMR_samples.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}

genelist=NULL
genelist["BRCAgene"]=list(c("BRCA1","BRCA2"))
genelist["MMR"]=list(c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2"))

comp=list(c("sMMR","NA"),c("sBRCAgene","NA"),c("gBRCAgene","NA"),c("gMMR","NA"))
ord=c("sMMR","sBRCAgene","gBRCAgene","NA")
#Neoantigene level
col_value=c("sMMR"="#56B4E9","gMMR"="#56B4E9","gBRCAgene"="#56B4E9","sBRCAgene"="#56B4E9","NA"="#dbd9d9")
cancer=c("BRCA","OV","UCEC","COADREAD","BLCA")#c(unique(as.character(contributionNew$ca)))
pheno="SNV_Neoantigens"
immuneFeature("Neodist_",cancer,pheno,short_name="log2(SNV Neoantigens)",labpos=15,genelist=genelist,log=TRUE,ord=ord,th=th,col_value=col_value,w=7,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COADREAD","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="Indel_Neoantigens"
ord=c("sMMR","gMMR","sBRCAgene","gBRCAgene","NA")

immuneFeature("Neodist_",cancer,pheno,short_name="log2(Indel Neoantigens)",labpos=15,genelist=genelist,ord=ord,col_value=col_value,log=TRUE,th=th,w=6,h=3.6)

#PD1,PDL1 level
cancer=c("BRCA","OV","UCEC")#c(unique(as.character(contributionNew$ca)))
pheno="PDCD1"
ord=c("sMMR","gMMR","sBRCAgene","gBRCAgene","NA")
immuneFeature("Neodist_",cancer,pheno,short_name="PD1",labpos=13,genelist=genelist,log=FALSE,ord=ord,col_value=col_value,th=th,w=4,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COAD","BLCA","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="CD274"
immuneFeature("Neodist_",cancer,pheno,short_name="PD-L1",labpos=110,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=5.5,h=3.6,comparisions =comp)


cancer=c("BRCA","OV","UCEC","COAD","BLCA","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="TIL_Regional_Fraction"
immuneFeature("Neodist_",cancer,pheno,short_name="TIL_Regional_Fraction",labpos=110,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=5.5,h=3.6,comparisions =comp)


cancer=c("UCEC")#c(unique(as.character(contributionNew$ca)))
pheno="Leukocyte_Fraction"
immuneFeature("Neodist_",cancer,pheno,short_name="Leukocyte_Fraction",labpos=1,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=2.2,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COAD","BLCA","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="BCR_Shannon"
immuneFeature("Neodist_",cancer,pheno,short_name="BCR_Shannon",labpos=110,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=5.5,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COAD","BLCA","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="TCR_Shannon"
immuneFeature("Neodist_",cancer,pheno,short_name="TCR_Shannon",labpos=110,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=5.5,h=3.6,comparisions = comp)


cancer=c("BRCA","UCEC")#c(unique(as.character(contributionNew$ca)))
pheno="CYTScore"
immuneFeature("Neodist_",cancer,pheno,short_name="CYTScore",labpos=15,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=3,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COAD","BLCA","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="CTLScore"
immuneFeature("Neodist_",cancer,pheno,short_name="CTLScore",labpos=110,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=5.5,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COADREAD","BLCA","PAAD")#c(unique(as.character(contributionNew$ca)))
pheno="TIL"
immuneFeature("Neodist_",cancer,pheno,short_name="TIL",labpos=50,genelist=genelist,log=FALSE,th,ord=ord,col_value=col_value,w=5.5,h=3.6,comparisions =comp)



##################20200205#######################
########individual gene distribution#############
genelist=NULL
genelist["BRCA1"]="BRCA1"
genelist["BRCA2"]="BRCA2"
genelist["MLH1"]="MLH1"
genelist["MSH2"]="MSH2"
genelist["MSH3"]="MSH3"
genelist["MSH6"]="MSH6"
genelist["PMS1"]="PMS1"
genelist["PMS2"]="PMS2"

col_value=NULL
col_value["BRCA1"]="#56B4E9"
col_value["BRCA2"]="#56B4E9"
col_value["MLH1"]="#56B4E9"
col_value["MSH2"]="#56B4E9"
col_value["MSH3"]="#56B4E9"
col_value["MSH6"]="#56B4E9"
col_value["PMS1"]="#56B4E9"
col_value["PMS2"]="#56B4E9"
col_value["NA"]="#dbd9d9"

ord=c(paste0("s",c("BRCA1","BRCA2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")),"NA")
#Neoantigene level
#col_value= as.character(paste0('"',ord,'"','="#dbd9d9"'))

cancer=c("BRCA","OV","UCEC","COADREAD","BLCA","LUAD","LUSC")#c(unique(as.character(contributionNew$ca)))
pheno="SNV_Neoantigens"
immuneFeature("IndividualgeneNeodist_",cancer,pheno,short_name="log2(SNV Neoantigens)",labpos=15,genelist=genelist,log=TRUE,ord=ord,th=th,col_value=col_value,w=7,h=3.6,comparisions = comp)


cancer=c("BRCA","OV","UCEC","COADREAD","BLCA","LUAD","LUSC")#c(unique(as.character(contributionNew$ca)))
pheno="PDCD1"
immuneFeature("IndividualgenePDCD1_",cancer,pheno,short_name="PDCD1",labpos=15,genelist=genelist,log=TRUE,ord=ord,th=th,col_value=col_value,w=7,h=3.6,comparisions = comp)

cancer=c("BRCA","OV","UCEC","COADREAD","BLCA","LUAD","LUSC")#c(unique(as.character(contributionNew$ca)))
pheno="CYTScore"
immuneFeature("IndividualgeneCYTScore_",cancer,pheno,short_name="CYTScore",labpos=15,genelist=genelist,log=TRUE,ord=ord,th=th,col_value=col_value,w=7,h=3.6,comparisions = comp)









