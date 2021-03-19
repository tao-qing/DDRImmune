setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/DDRHypermutator/")

source("../dependency_files_tq.R")
source("../load_somatic.R")

library(ggplot2)
library(ggsignif)
library(ggrepel)


allspl=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]
hyper=read.table("../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/Driver_BaileyCell2018/HypermutatorSamples.txt",h=T,sep="\t",stringsAsFactors=FALSE)[,1]


###############################################################
####################1. data prepare###########################
##############################################################
#immune score
immuneProfile<-read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile<-immuneProfile[-which(is.na(immuneProfile$value)),]

#immune response features for consideration
response=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)

immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%allspl),]
pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%allspl),]



#MSISensor
cancer=c("STAD","UCEC","COADREAD","BRCA")#c("STAD","UCEC","COADREAD","LUSC","LUAD","SKCM","BRCA","BLCA","CESC","HNSC")

pheno="Leukocyte_Fraction"
genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR","CHEK2")
#genes=c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")
#col_value=c("both"="#eb0e0e","germline"="#071eed","somatic"="#f57c02","NA"="#dbd9d9")

immuneFeature<-function(label,cancer,pheno,short_name,gene,log,lpos,hpos,col_value,th,w,h){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable%in%pheno),]
  #sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  if(log==TRUE){
    subimmnue$value=log2(subimmnue$value+0.1)
  }
  mutMat=as.data.frame(clin[clin$type%in%cancer,c("bcr_patient_barcode","type")])
  
  for(ge in genes){
    smutspl=unique(c(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%ge]))
    gmutspl=unique(c(pathVar$bcr_patient_barcode[pathVar$HUGO_Symbol%in%ge]))
    allmutspl=unique(c(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%genes],pathVar$bcr_patient_barcode[pathVar$HUGO_Symbol%in%genes]))
    eval(parse(text=paste0("mutMat$s",ge,"<-ifelse(mutMat$bcr_patient_barcode%in%smutspl,\"MUT\",ifelse(mutMat$bcr_patient_barcode%in%allmutspl,NA,\"WT\"))")))
    eval(parse(text=paste0("mutMat$g",ge,"<-ifelse(mutMat$bcr_patient_barcode%in%gmutspl,\"MUT\",ifelse(mutMat$bcr_patient_barcode%in%allmutspl,NA,\"WT\"))")))
  }
  
 phenoMat=NULL
 for(ge in c(paste0("g",genes),paste0("s",genes))){
    for(ca in cancer){
        mutspl=mutMat$bcr_patient_barcode[mutMat$type==ca & mutMat[,ge]=="MUT"]
        wtspl=mutMat$bcr_patient_barcode[mutMat$type==ca & mutMat[,ge]=="WT"]
        mutvalue=subimmnue$value[which(subimmnue$bcr_patient_barcode%in%mutspl)]
        wtvalue=subimmnue$value[which(subimmnue$bcr_patient_barcode%in%wtspl)]
        if(length(mutvalue)<4 || length(wtvalue)<4){
          next
        }
        mutmean=mean(mutvalue)
        wtmean=mean(wtvalue)
        pvalue=t.test(mutvalue,wtvalue)$p.value
        phenoMat=rbind(phenoMat,cbind(Cancer=ca,Gene=ge,MUT=mutmean,WT=wtmean,P=pvalue))
    }  
 }
 
 phenoMat=as.data.frame(phenoMat)
 phenoMat$P=as.numeric(as.matrix(phenoMat$P))
 phenoMat$MUT=as.numeric(as.matrix(phenoMat$MUT))
 phenoMat$WT=as.numeric(as.matrix(phenoMat$WT))
 phenoMat$FDR=p.adjust(phenoMat$P,method="fdr")
 phenoMat$Sig=ifelse(phenoMat$FDR<0.05,"Significant (FDR<0.05)",ifelse(phenoMat$FDR<0.15,"Suggestive (FDR<0.15)","None"))
 
 return_mat=phenoMat
 #th= theme_bw()+ theme(legend.key.size =unit(.2, "cm"),legend.title = element_text(size=8) ,panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=14,hjust = 0.95), axis.text.y = element_text(colour="black", size=14,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=14,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.5, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))
 
 phenoMatSub=phenoMat[phenoMat$FDR<0.15,]
 phenoMatSub$label=paste0(phenoMatSub$Cancer,":",phenoMatSub$Gene)
   
 p = ggplot(phenoMatSub, aes(y=WT, x=MUT)) #
 p = p + geom_point(aes(color=factor(phenoMatSub$Sig)),alpha=0.5,size=-log10(phenoMatSub$P))+xlim(lpos,hpos)+ylim(lpos,hpos)#+scale_shape_manual(values=c(1, 5))#,shape =Type$ ,alpha = I(0.6)
 p = p + geom_text_repel(aes(label=ifelse(phenoMatSub$FDR<0.05,label,"")),size=3,segment.alpha =0.5,fontface="bold.italic",segment.colour ="grey",min.segment.length = 0) 
 #p = p+facet_grid(.~Gene,drop=T)#, space="free",scale="free"
 p = p + ylab(paste0(short_name," in WT")) + xlab(paste0(short_name," in MUT")) + theme_bw()
 p = p + geom_abline(intercept = 0, slope=1, alpha=0.2)
 
 p = p + guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA)))+ theme(plot.title = element_text(size=14, face="bold"),axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14)) + labs(size = "-log10(P)") +scale_size_continuous(limits=c(0,50),breaks = c(1,10,20,30,40,50))
 #p = p + scale_colour_manual("", values = c("Significant (FDR<0.05)" = "red",  "Suggestive (FDR<0.15)" = "blue","None" = "grey"))
 p = p + scale_color_manual("FDR<0.05", values = c("Significant (FDR<0.05)" = "red","Suggestive (FDR<0.15)"="blue","None" = "grey"))
# p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
# p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
 p=p+th
 p
 fn = paste0("out/",label,"_",pheno,"_DDRgeneMutvsWT_20200211.pdf")
 ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
 
 return(return_mat)
}

th= theme_bw()+ theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=12,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))


pheno="CD274"
#genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR")
cd274mat=immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=FALSE,lpos=3,hpos=7,short_name="PD-L1",gene=genes,th=th,w=3,h=3)

pheno="PDCD1"
#genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR")
pd1mat=immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=FALSE,lpos=4,hpos=8,short_name="PD1",gene=genes,th=th,w=3,h=3)


#pheno=c("xCell_CD8+_T-cells","xCell_CD4+_T-cells","xCell_B-cells","xCell_Th1_cells","xCell_Th2_cells","xCell_DC")
#immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=FALSE,lpos=4,hpos=8,short_name="ImmuneCellSignature",gene=genes,th=th,w=4,h=3)

pheno="Lymphocyte_Infiltration_Signature_Score"
#genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR")
tilmat=immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=FALSE,lpos=-1,hpos=2,short_name="TIL",gene=genes,th=th,w=3,h=3)


th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=12,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

pheno="CYTScore"
#genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR")
cytmat=immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=FALSE,lpos=6,hpos=9,short_name="CYTScore",gene=genes,th=th,w=5.1,h=3)


mergemat=rbind(pd1mat,cd274mat,cytmat,tilmat)
write.csv(mergemat)

length(mergemat$Sig[mergemat$Sig=="Significant (FDR<0.05)"])

#cancer=unique(clin$type)#[clin$type%in%cancer])
pheno="SNV_Neoantigens"
genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR")
immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=TRUE,lpos=5,hpos=15,short_name="SNV Neoantigens",gene=genes,th=th,w=3,h=3)

#cancer=unique(clin$type[clin$type%in%cancer])
pheno="Indel_Neoantigens"
genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR")
immuneFeature("Figure_",cancer=cancer,pheno=pheno,log=TRUE,lpos=3,hpos=8.5,short_name="Indel Neoantigens",gene=genes,th=th,w=3,h=3)

  