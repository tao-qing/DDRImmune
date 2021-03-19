setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/")

#draw a scatter plot to compare the coefficient from germline and somatic level. 

###################################Data Prepare###################################
source("../global_aes_out.R")
source("../dependency_files_tq.R")
library(grid)
library("scales")
library("ggplot2")
library("ggrepel")
library(magrittr)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
#wrapping too long facet title
swr = function(string, nwrap=10) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

#load data
#remove hyper mutator: TRUE of FALSE
rmhypermutator=TRUE
#immune score
immuneProfile=read.table("./out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]

#cancer samples
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]

#genomic and immune signatures
  immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]
  pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%samples),]
  datatype="overlaped_n9738"




###################################Figures###################################
xlabel=scale_x_discrete(position = "right",labels=c("Aneuploidy_Score" = "Aneuploidy", "BCR_Evenness" = "BCR_E","BCR_Richness" = "BCR_R","BCR_Shannon"="BCR_S","CTA_Score"="CTA","Homologous_Recombination_Defects"="HRD","TIL_Regional_Fraction"="TIL","Leukocyte_Fraction"="Leukocyte","Macrophage_Regulation"="Macrophage","Lymphocyte_Infiltration_Signature_Score"="Lymphocyte","IFN_gamma_Response"="IFN","TGF_beta_Response"="TGF_beta","Th1_Cells"="Th1","Th2_Cells"="Th2","Th17_Cells"="Th17"))
dotshape=scale_shape_manual(labels=c("Significant (FDR<0.05)","Suggestive (FDR<0.15)","None"),values=c(17,15,16))


########################Figure 2############################
pheno=c("DNA-repair expression signature","DNA-repair mutational signature")
allIterm=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms",sep="\t")
iterm=allIterm$V2[which(allIterm$V1%in%pheno)]

#all core DDR genes
ddr=read.csv("../../Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/DDR_Pathways.csv",h=T)
ddrgene=unique(as.character(unlist(ddr)))[-1]
load("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/DDR_NanoString_Pathways.RData")

#th= theme_bw()+ theme(legend.position =legend,legend.key.size =unit(.2, "cm"),legend.title = element_text(size=8) ,panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
#th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black", size=14, angle = 60,hjust=0.90), axis.text.y = element_text(colour="black", size=14,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=18,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))
gscontributiongene<-function(pheno,label,llim,hlim,w,h,th=th,legend="none"){
  ####load data####
  #germline 
  type="pathVarP_"
  #datatype="_"#_rmhypermutator_
  fn=paste0("./out/",type,"ImmuneAssoc_genelevel",datatype,"byCancerType_cov_20200409.txt")
  gassoc = read.table(sep="\t",header=T,file=fn, stringsAsFactors=FALSE)
  gassoc = gassoc[gassoc$gene_path_count>=4,]
  if(length(which(is.na(gassoc$p.value)))!=0){
    gassoc = gassoc[-which(is.na(gassoc$p.value)),]
  }
  
  
 #select variables 
  gsubAssoc=gassoc[which(gassoc$clinicPhenotype%in%pheno),]
  gsubAssoc=gsubAssoc[gsubAssoc$gene%in%ddrgene,]
  gsubAssoc$Type="germline"
  #FDR
  gsubAssoc$FDR=p.adjust(gsubAssoc$p.value,method="fdr")
  gsubAssoc$Association = "None" 
  gsubAssoc$Association[gsubAssoc$FDR<0.15 ] = "Suggestive (FDR<0.15)"
  gsubAssoc$Association[gsubAssoc$FDR<0.05 ] = "Significant (FDR<0.05)"
  gsubAssoc$Association<-factor(gsubAssoc$Association,levels=c("Significant (FDR<0.05)","Suggestive (FDR<0.15)","None")) 

  #somatic
  type="somaticDriver_"
  #datatype="_"#_rmhypermutator_
  fn=paste0("./out/",type,"ImmuneAssoc_genelevel",datatype,"byCancerType_cov_20200409.txt")
  sassoc = read.table(sep="\t",header=T,file=fn, stringsAsFactors=FALSE)
  sassoc = sassoc[sassoc$gene_path_count>=4,]
  if(length(which(is.na(sassoc$p.value)))!=0){
    sassoc = sassoc[-which(is.na(sassoc$p.value)),]
  }
  #select variables 
  ssubAssoc=sassoc[which(sassoc$clinicPhenotype%in%pheno),]
  ssubAssoc=ssubAssoc[ssubAssoc$gene%in%ddrgene,]
  ssubAssoc$Type="somatic"
  #FDR
  ssubAssoc$FDR=p.adjust(ssubAssoc$p.value,method="fdr")
  ssubAssoc$Association = "None" 
  ssubAssoc$Association[ssubAssoc$FDR<0.15 ] = "Suggestive (FDR<0.15)"
  ssubAssoc$Association[ssubAssoc$FDR<0.05 ] = "Significant (FDR<0.05)"
  ssubAssoc$Association<-factor(ssubAssoc$Association,levels=c("Significant (FDR<0.05)","Suggestive (FDR<0.15)","None")) 
  
  gsubAssoc1=gsubAssoc
  ssubAssoc1=ssubAssoc
##########Construct a new data frame to show the scatter plot######
  gsubAssoc1$id=paste0(gsubAssoc1$cancer,"|",gsubAssoc1$gene,"|",gsubAssoc1$clinicPhenotype)
  ssubAssoc1$id=paste0(ssubAssoc1$cancer,"|",ssubAssoc1$gene,"|",ssubAssoc1$clinicPhenotype)
  id=unique(c(gsubAssoc1$id,ssubAssoc1$id))
  newMat=as.data.frame(cbind(id=id,cancer=as.character(sapply(id,function(x)strsplit(x,split="\\|")[[1]][1])),gene=as.character(sapply(id,function(x)strsplit(x,split="\\|")[[1]][2])),pheno=as.character(sapply(id,function(x)strsplit(x,split="\\|")[[1]][3]))))

  newMat$gcontribution=gsubAssoc1$coefficient[pmatch(newMat$id,gsubAssoc1$id)]
  newMat$scontribution=ssubAssoc1$coefficient[pmatch(newMat$id,ssubAssoc1$id)]
  newMat$gcontribution=ifelse(is.na(newMat$gcontribution),0,newMat$gcontribution)
  newMat$scontribution=ifelse(is.na(newMat$scontribution),0,newMat$scontribution)
  
  newMat$gSig=as.character(gsubAssoc1$Association)[pmatch(newMat$id,gsubAssoc1$id)]
  newMat$sSig=as.character(ssubAssoc1$Association)[pmatch(newMat$id,ssubAssoc1$id)]
  
  #P value
  newMat$gP=gsubAssoc1$p.value[pmatch(newMat$id,gsubAssoc1$id)]
  newMat$sP=ssubAssoc1$p.value[pmatch(newMat$id,ssubAssoc1$id)]
  
  newMat$gFDR=gsubAssoc1$FDR[pmatch(newMat$id,gsubAssoc1$id)]
  newMat$sFDR=ssubAssoc1$FDR[pmatch(newMat$id,ssubAssoc1$id)]
  
  newMat$gse=gsubAssoc1$StdError[pmatch(newMat$id,gsubAssoc1$id)]
  newMat$sse=ssubAssoc1$StdError[pmatch(newMat$id,ssubAssoc1$id)]
  
  newMat$label=paste0(newMat$cancer,":",newMat$gene)
  
  ####add categorie####
  if(label=="HRD"){
    BRCAgene=c("BRCA1","BRCA2")
    HRgene=c("ATM")
    #MMR=c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")
    MMR=pathlist$DDR$CORE_Mismatch_Repair
    #OtherDDR=ddrgene[-which(ddrgene%in%c("PALB2",BRCAgene,MMR,HRgene))]#
    OtherHR=c("PALB2","BRIP1","FANCA","FANCM","RAD51C","RAD51D","RAD54L","HDAC2","CHEK2","NBN")
    
    #tmp9=pathVarP[pathVarP$HUGO_Symbol%in%OtherHR,]
    newMat$Cate=ifelse(newMat$gene=="BRCA1","BRCA1",ifelse(newMat$gene=="BRCA2","BRCA2",ifelse(newMat$gene%in%"ATM","ATM",ifelse(newMat$gene%in%OtherHR,"Other HR genes","NA"))))
    #newMat1$Cate=ifelse(newMat1$gene=="BRCA1","BRCA1",ifelse(newMat1$gene=="BRCA2","BRCA2",ifelse(newMat1$gene%in%"ATM","ATM",ifelse(newMat1$gene=="PALB2","PALB2",ifelse(newMat1$gene%in%"FANCM","FANCM","NA")))))
    newMat=newMat[-which(newMat$Cate=="NA"),]
    newMat$Cate=factor(newMat$Cate,levels=c("BRCA1","BRCA2","ATM","Other HR genes"))
    for(g1 in c(BRCAgene,HRgene)){
      newMat$label=gsub(paste0(":",g1),"", newMat$label)
    }
    
  }else{
    BRCAgene=c("BRCA1","BRCA2")
    #HRgene=c("ATM","ATR")
    MMR=c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")
    #MMR=pathlist$DDR$CORE_Mismatch_Repair
    #NER=pathlist$DDR$CORE_Nucleotide_Excision_Repair
    #TS=pathlist$DDR$CORE_Translesion_Synthesis
    OtherDDR=ddrgene[-which(ddrgene%in%c(BRCAgene,MMR))]
    newMat$Cate=ifelse(newMat$gene=="BRCA1","BRCA1",ifelse(newMat$gene=="BRCA2","BRCA2",ifelse(newMat$gene%in%MMR,"Mismatch repair",ifelse(newMat$gene%in%OtherDDR,"Other DDR","Other DDR"))))
    newMat$Cate=factor(newMat$Cate,levels=c("BRCA1","BRCA2","Mismatch repair","Other DDR"))
    
    # newMat1=newMat2
    newMat$label[which(newMat$Cate=="Other DDR" & newMat$scontribution<0.5 & newMat$gcontribution<=0.5)]=""
    newMat$label[which(newMat$Cate=="Mismatch repair" & newMat$scontribution<0.5 & newMat$gcontribution<=0.5)]=""
    
    for(g1 in c(BRCAgene)){
      newMat$label=gsub(paste0(":",g1),"", newMat$label)
    }
  }
  
  newMat0=newMat[newMat$pheno==pheno,]
  
  anysig=as.character(unique(newMat0$id[newMat0$gSig%in%c("Significant (FDR<0.05)","Suggestive (FDR<0.15)") | newMat0$sSig%in%c("Significant (FDR<0.05)","Suggestive (FDR<0.15)") ]))
  newMat1=newMat0[newMat0$id%in%anysig,]
 
  newMat1$gP=ifelse(newMat1$gP<1e-50,1e-50,newMat1$gP)
  newMat1$sP=ifelse(newMat1$sP<1e-50,1e-50,newMat1$sP)
  
  newMat1$P=ifelse(is.na(newMat1$gP) & is.na(newMat1$sP),1,ifelse(is.na(newMat1$gP),newMat1$sP,ifelse(is.na(newMat1$sP),newMat1$gP,ifelse(newMat1$gP>newMat1$sP,newMat1$sP,newMat1$gP))))
  
  newMat1$gFDR=ifelse(newMat1$gFDR<1e-50,1e-50,newMat1$gFDR)
  newMat1$sFDR=ifelse(newMat1$sFDR<1e-50,1e-50,newMat1$sFDR)
  
  newMat1$FDR=ifelse(is.na(newMat1$gFDR) & is.na(newMat1$sFDR),1,ifelse(is.na(newMat1$gFDR),newMat1$sFDR,ifelse(is.na(newMat1$sFDR),newMat1$gFDR,ifelse(newMat1$gFDR>newMat1$sFDR,newMat1$sFDR,newMat1$gFDR))))
  
  
  #add stand error 
  newMat1$x_l=newMat1$gcontribution-newMat1$sse
  newMat1$x_h=newMat1$gcontribution+newMat1$sse
  newMat1$y_l=newMat1$scontribution-newMat1$gse
  newMat1$y_h=newMat1$scontribution+newMat1$gse
    
  #label show on scatter plot 
  #newMat1$label=paste0(newMat1$cancer,":",newMat1$gene)
  
  #add color 
  newMat1$color="none"
  gsig=newMat1$id[which( newMat1$gSig=="Significant (FDR<0.05)")]
  ssig=newMat1$id[which( newMat1$sSig=="Significant (FDR<0.05)")]
 
  overlap=intersect(gsig,ssig)
  gonly=gsig[!gsig%in%overlap]
  sonly=ssig[!ssig%in%overlap]
  newMat1$color[newMat1$id%in%gonly]="germline"
  newMat1$color[newMat1$id%in%sonly]="somatic"
  newMat1$color[newMat1$id%in%overlap]="both"
  
  newMat1$color=factor(newMat1$color,levels=c("both","germline","somatic","none"))
  #newMat1$label[newMat1$gSig!="Significant (FDR<0.05)" & newMat1$sSig!="Significant (FDR<0.05)"]=""
  newMat1$label[newMat1$color=="none"]=""
  
  
  #for germline and somatic, dot size only shows the most significant one
#  ids=unique(as.character(newMat1$id))
#  newMat2=newMat1
#  for(id in ids){
#    index=which(newMat1$id==id)
#    submet=newMat1[index,]
#    if(submet$P[1]>submet$P[2]){
#      newMat2[index[1],"P"]=1
#    }else{
#      newMat2[index[2],"P"]=1
#    }
#  }

  p = ggplot(newMat1, aes(y=scontribution, x=gcontribution)) 
  p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
  p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
  p = p + geom_abline(intercept = 0, slope=1, alpha=0.2)
  #p = p + geom_point(aes(colour=newMat1$Sig,size=-log10(newMat1$P),shape =Type) )+scale_shape_manual(values=c(1,5))+xlim(llim,hlim)+ylim(llim,hlim)
  p = p + geom_point(aes(color=newMat1$color,size=-log10(newMat1$FDR)),alpha=I(0.6))+scale_shape_manual(values=c(1,5))+xlim(llim,hlim)+ylim(llim,hlim)
  #p = p + geom_errorbar(aes(ymin = y_l,ymax = y_h),color="lightblue")
  #p = p + geom_errorbarh(aes(xmin = x_l ,xmax = x_h),color="lightblue")
  p = p + geom_text_repel(aes(label=newMat1$label),size=5,segment.alpha =0.5,segment.colour ="grey",min.segment.length = 0)
  #fill=ifelse(newMat1$label=="","",as.character(newMat1$cancer))
  #p = p + getPCACancerColor()
  p = p+facet_grid(.~Cate,drop=T)#, space="free",scale="free"
  p = p + ylab(paste0("Somatic coefficients (",label,")")) + xlab(paste0("Germline coefficients (",label,")")) + theme_bw() 
 # p = p +theme_bw()+  guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA)))+ theme(plot.title = element_text(size=14, face="bold"),axis.title = element_text(size=18), axis.text.x = element_text(colour="black", size=18), axis.text.y = element_text(colour="black", size=18))
  p = p + labs(size = "-log10(FDR)") +scale_size_continuous(limits=c(0,50),breaks = c(1,20,30,40,50))
 #p = p + scale_colour_manual("", values = c("Significant (FDR<0.05)" = "red",  "Suggestive (FDR<0.15)" = "blue","None" = "black"))
  p = p + scale_color_manual("FDR<0.05", values = c("both" = "red",  "germline" = "blue","somatic"="orange","none" = "grey"))
  p=p+th
  p
  fn = paste0("out/Figure2_",label,"_germline_somatic_contribution.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
  return(newMat)
}
w=16
h=4

th= theme_bw()+theme(legend.position = "none",axis.text.x = element_text(colour="black", size=20,hjust = 0.95), axis.text.y = element_text(colour="black", size=20,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=18,face="bold"))+ theme(strip.placement = "outside",strip.text.x = element_text(size=18),strip.text.y = element_text(angle = 0,size=20,face="italic"),panel.spacing = unit(0.5, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#hrd=gscontributiongene("Homologous_Recombination_Defects","HRD",llim=-0.5,hlim=1.0,w=w,h=h)
tmb=gscontributiongene("TMB","TMB",w=w,h=h,llim=-2.0,hlim=5.0,th=th,legend="none")
snv=gscontributiongene("SNV_Neoantigens","SNV Neo",w=w,h=h,llim=-1,hlim=3.5,th=th,legend="none")
indel=gscontributiongene("Indel_Neoantigens","Indel Neo",w=w,h=h,llim=-1,hlim=2.5,th=th,legend="none")

merge=rbind(tmb,snv,indel)
write.csv(merge,"./out/figure2_allAssoc_suppleTable2.csv",quote=F)

