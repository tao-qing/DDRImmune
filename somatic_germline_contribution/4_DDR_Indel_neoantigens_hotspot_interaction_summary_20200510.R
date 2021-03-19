setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")
source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../stat_functions.R")
source("../load_somatic.R")

samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]
subclin=clin[clin$bcr_patient_barcode%in%samples,]
cancerSize=table(subclin$type)
###############################################################
####################1. data prepare###########################
##############################################################
#remove hyper mutator: TRUE of FALSE
rmhypermutator=TRUE
#immune score
immuneProfile<-read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile<-immuneProfile[-which(is.na(immuneProfile$value)),]


#msi high samples (MSISensor>=4)
msi=immuneProfile[immuneProfile$variable=="MSISensor",]
msispl=msi$bcr_patient_barcode[msi$value>=4]

#immune response features for consideration
response=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)

if(rmhypermutator){
  hyper<-read.table("../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/Driver_BaileyCell2018/HypermutatorSamples.txt",h=T,sep="\t",stringsAsFactors=FALSE)[,1]
  immuneProfile=immuneProfile[-which(immuneProfile$bcr_patient_barcode%in%hyper),]
  pathVarP=pathVarP[-which(pathVarP$bcr_patient_barcode%in%hyper),]
  label="_rmhypermutator_"
}else{
  label="_"
}


mmrgene=c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")
gbrcaspl=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%c("BRCA1","BRCA2","PALB2")])
sbrcaspl=unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%c("BRCA1","BRCA2","PALB2")])

gpolespl=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%c("POLE","POLQ")])
spolespl=unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%c("POLE","POLQ")])

gatmspl=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%c("ATM","ATR","CHEK2")])
satmspl=unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%c("ATM","ATR","CHEK2")])

gbrcaspl=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%c("BRCA1","BRCA2","PALB2")])
sbrcaspl=unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%c("BRCA1","BRCA2","PALB2")])

gmmrspl=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%mmrgene])
smmrspl=unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%mmrgene])


library(data.table)
neoantigen<-fread("~/Box Sync/Others/controledData/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.indel_neoantigens_10062017.tsv",data.table = F,stringsAsFactors = F)
neoantigen$cancer_type=gsub("COAD|READ","COADREAD",neoantigen$cancer_type)
tmp0=neoantigen[-which(duplicated(neoantigen$indel)),]
tmp=as.data.frame(table(neoantigen$indel,neoantigen$cancer_type))
tmp$Var1=as.character(tmp$Var1)
tmp$Var2=as.character(tmp$Var2)
tmp=cbind(tmp,gene_name=tmp0[unlist(sapply(tmp$Var1,function(x)which(tmp0$indel==x))),"gene_name"])


#colnames(tmpSort)=gsub("Freq","Count",colnames(tmp))
#tmpSort=tmpSort[order(tmpSort$Count,decreasing = T),]
#tmpSort$Freq=tmpSort$Count/cancerSize[tmpSort$Var2]
#write.table(tmpSort,"./out/hotspot_somatic_Indel_Neoantigens.csv")


#hotspot variants
th= theme_bw()+ theme(legend.position = c(0.8, 0.95),legend.key.size =unit(.2, "cm"),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black", size=8,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=10,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

#plotmat1=tmp[tmp$Freq>=5,]
#plotmat1$mut=factor(plotmat1$mut,levels=as.character(plotmat1$mut[order(as.numeric(as.matrix(plotmat1$Freq)),decreasing = F)]))

#p=ggplot(plotmat1) + geom_bar(aes(mut,as.numeric(as.matrix(Freq))), stat='identity', alpha=0.5)+coord_flip()
#p=p+ ylab("# of carriers")+ xlab("Predicted Indel Neoantigens")+ggtitle("")+ labs(fill = "Types")
#p=p+th
#ggsave("./out/number_predictedIndelNeoantigenes_TCGA.pdf",useDingbat=F,width=4,height=6)

#write.csv(tmp,"./out/hotspot_somatic_Indel_Neoantigens.csv",quote=F)

#whether some specific neoantigens enriched in germline and somatic mutant samples
#focus on BRCA1/2, MMR genes.
#gene level, hotspots level
sample=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",stringsAsFactors = F)[,1]
#sample=sample[-which(sample%in%hyper)]
mat=clin[clin$bcr_patient_barcode%in%sample,c("bcr_patient_barcode","type")]
mat$gHR=ifelse(mat$bcr_patient_barcode%in%gbrcaspl,1,0)
mat$sHR=ifelse(mat$bcr_patient_barcode%in%sbrcaspl,1,0)
mat$gDP=ifelse(mat$bcr_patient_barcode%in%gpolespl,1,0)
mat$sDP=ifelse(mat$bcr_patient_barcode%in%spolespl,1,0)
mat$gDS=ifelse(mat$bcr_patient_barcode%in%gatmspl,1,0)
mat$sDS=ifelse(mat$bcr_patient_barcode%in%satmspl,1,0)
mat$gMMR=ifelse(mat$bcr_patient_barcode%in%gmmrspl,1,0)
mat$sMMR=ifelse(mat$bcr_patient_barcode%in%smmrspl,1,0)

subneoantigen=neoantigen[neoantigen$indel%in%tmp$Var1[which(tmp$Freq>=4)],]
subneoantigen=as.data.frame(table(subneoantigen$sample,subneoantigen$indel))
subneoantigen$Var1=as.character(subneoantigen$Var1)
subneoantigen$Var2=as.character(subneoantigen$Var2)
subneoantigen$Freq=as.numeric(subneoantigen$Freq)

mutstat=NULL
for(ge in c("gHR","sHR","gMMR","sMMR","gDP","sDP","gDS","sDS")){
  for(mut in tmp$Var1[which(tmp$Freq>=4)]){
    submat=mat
    mutspl=unique(subneoantigen$Var1[which(subneoantigen$Var2==mut & subneoantigen$Freq!=0)])
    if(!any(submat$bcr_patient_barcode%in%mutspl)){next}
    submat$Type=ifelse(submat$bcr_patient_barcode%in%mutspl,1,0)
    mat_tmp=table(submat$Type,submat[,ge])
    mutnoneo=mat_tmp["0","1"]
    wtnoneo=mat_tmp["0","0"]
    mutneo=mat_tmp["1","1"]
    wtneo=mat_tmp["1","0"]
    p=fisher.test(submat$Type,submat[,ge])$p.value
    or=as.numeric(fisher.test(submat$Type,submat[,ge])$estimate)
    mutstat=rbind(mutstat,cbind(gene=ge,Indel=mut,mutationWithNeo=mutneo,wildwithNeo=wtneo,mutationNon=mutnoneo,wildNon=wtnoneo,OR=or,P=p))
  }
}

mutstat=as.data.frame(mutstat)
mutstat$mutGene=as.character(sapply(mutstat$Indel,function(x)tmp0$gene_name[which(tmp0$indel==x)]))
mutstat$P=as.numeric(as.matrix(mutstat$P))
mutstat$FDR=p.adjust(mutstat$P,method="fdr")
write.table(mutstat,"./out/pancancer_mutation_Indel_neoantigen_association.txt",row.names = F,sep="\t")


mutstatcancer=NULL
for(ca in unique(clin$type)){
  for(ge in c("gBRCAgene","sBRCAgene","gMMR","sMMR","gPOLE_POLQ","sPOLE_POLQ","gATM_ATR","sATM_ATR")){
    for(mut in tmp$Var1[which(tmp$Freq>5)]){
      submat=mat[mat$type==ca,]
      mutspl=unique(subneoantigen$Var1[which(subneoantigen$Var2==mut & subneoantigen$Freq!=0)])
      submat$Type=ifelse(submat$bcr_patient_barcode%in%mutspl,1,0)
      if(sum(submat$Type)!=0 & sum(submat[,ge])!=0){
        mat_tmp=table(submat$Type,submat[,ge])
        mutnoneo=mat_tmp["0","1"]
        wtnoneo=mat_tmp["0","0"]
        mutneo=mat_tmp["1","1"]
        wtneo=mat_tmp["1","0"]
        p=fisher.test(submat$Type,submat[,ge])$p.value
        or=as.numeric(fisher.test(submat$Type,submat[,ge])$estimate)
        mutstatcancer=rbind(mutstatcancer,cbind(cancer=ca,gene=ge,Indel=mut,mutationWithNeo=mutneo,wildwithNeo=wtneo,mutationNon=mutnoneo,wildNon=wtnoneo,OR=or,P=p))
      }
    }
  }
}
mutstatcancer=as.data.frame(mutstatcancer)
mutstatcancer$mutGene=as.character(sapply(mutstatcancer$Indel,function(x)tmp0$gene_name[which(tmp0$indel==x)]))
mutstatcancer$P=as.numeric(as.matrix(mutstatcancer$P))

write.table(mutstatcancer,"./out/mutation_Indel_neoantigen_association_byCancer.txt",row.names = F,sep="\t")



