setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")

source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../stat_functions.R")
source("../load_somatic.R")
library(plspm)

#cancer samples
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt",h=F,stringsAsFactors = F)[,1]


###############################################################
####################1. data prepare###########################
##############################################################
#remove hyper mutator: TRUE of FALSE
rmhypermutator=TRUE
#immune score
immuneProfile<-read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile<-immuneProfile[-which(is.na(immuneProfile$value)),]

#immune response features for consideration
response=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)

#if(rmhypermutator){
#  hyper<-read.table("../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/Driver_BaileyCell2018/HypermutatorSamples.txt",h=T,sep="\t",stringsAsFactors=FALSE)[,1]
  immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]
  pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%samples),]
#  label="_rmhypermutator_"
#}else{
#  label="_"
#}

cancers=unique(immuneProfile$TCGA_Study)

genes_TMB=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS2")
genes_list_TMB=genes_TMB
names(genes_list_TMB)=genes_TMB

genes_HRD=c("BRCA1","BRCA2","PALB2","FANCM","ATM")
genes_list_HRD=genes_HRD
names(genes_list_HRD)=genes_HRD
###############################################################
################2. PLS-PM ANALYSIS####################
##############################################################
#pheno="MSISensor"#"Nonsilent_Mutation_Rate","MSISensor"

#function for extracting the results of PLS-PM
rus_pls_sum<-function(rus_pls,factor1="germline",factor2="somatic"){
  feature=rus_pls$inner_model$feature#significant test
  effects=rus_pls$effects#latent variable association
  loading=rus_pls$crossloadings#individuals gene contribution
  
  mainMat=NULL#germline->somatic->feature
  subMat=NULL#gene->germline, gene->somatic
  if(any(grepl(factor1,colnames(rus_pls$path_coefs))) & !any(grepl(factor2,colnames(rus_pls$path_coefs)))){
    mainMat=c(feature[factor1,"Pr(>|t|)"],0,effects[effects$relationships==paste0(factor1," -> feature"),"total"],0)
    names(mainMat)=c(paste0(c(factor1,factor2),"_P"),paste0(c(factor1,factor2),"_Contribution"))
    iterm=loading$name
    subMat=c(loading[grep("G:",iterm),factor1],0)
    names(subMat)=c(as.character(iterm[grep("G:",iterm)]),NA)
    
  }else if(!any(grepl(factor1,colnames(rus_pls$path_coefs))) & any(grepl(factor2,colnames(rus_pls$path_coefs)))){
  
    mainMat=c(0,feature[factor2,"Pr(>|t|)"],0,effects[effects$relationships==paste0(factor2," -> feature"),"total"])
    names(mainMat)=c(paste0(c(factor1,factor2),"_P"),paste0(c(factor1,factor2),"_Contribution"))
    iterm=loading$name
    subMat=c(0,loading[grep("S:",iterm),factor2])
    names(subMat)=c(NA,as.character(iterm[grep("S:",iterm)]))
  }else{
    mainMat=c(feature[factor1,"Pr(>|t|)"],feature[factor2,"Pr(>|t|)"],effects[effects$relationships==paste0(factor1," -> feature"),"total"],effects[effects$relationships==paste0(factor2," -> feature"),"total"])
    names(mainMat)=c(paste0(c(factor1,factor2),"_P"),paste0(c(factor1,factor2),"_Contribution"))
    iterm=loading$name
    subMat=c(loading[grep("G:",iterm),factor1],loading[grep("S:",iterm),factor2])
    names(subMat)=c(as.character(iterm[grep("G:",iterm)]),as.character(iterm[grep("S:",iterm)]))
  }
  
  return(list(main=mainMat,loading=subMat))
}

mainContribution=NULL
subContribution=NULL
for(pheno in c("TMB","MSISensor","SNV_Neoantigens","Indel_Neoantigens","Homologous_Recombination_Defects")){#"TMB","MSISensor","SNV_Neoantigens",
    for(cancer in cancers){
      samples=clin$bcr_patient_barcode[clin$type==cancer]
      mutMat=samples
      latent=NULL
      gval=NULL
      sval=NULL
      counts=NULL
      
      subimmune=immuneProfile[immuneProfile$variable==pheno,]
      if(length(which(is.na(subimmune$value)))>0){
        subimmune=subimmune[-which(is.na(subimmune$value)),]
      }
      
      if(pheno=="Homologous_Recombination_Defects"){
        genes_list=genes_list_HRD
      }else{
        genes_list=genes_list_TMB
      }
      #germline mutation
      for(g in names(genes_list)){
        gsample=pathVarP$bcr_patient_barcode[which(pathVarP$HUGO_Symbol==g)]
        
        gstatus=ifelse(samples%in%gsample,1,0)
        gcarrier=intersect(samples[samples%in%gsample],subimmune$bcr_patient_barcode)
        
        if(length(gcarrier)>3){
          eval(parse(text=(paste0("mutMat=cbind(mutMat,\"G:",g,"\"=gstatus)"))))
          latent=c(latent,"germline")
          gval=c(gval,paste0("G:",g))
        }
      }
      
      #somatic mutation
      for(s in names(genes_list)){
        ssample=somatic_likelyfunctional_driver$bcr_patient_barcode[which(somatic_likelyfunctional_driver$Hugo_Symbol==s)]
        sstatus=ifelse(samples%in%ssample,1,0)
        scarrier=intersect(samples[samples%in%ssample],subimmune$bcr_patient_barcode)
        
        if(length(scarrier)>3){
          eval(parse(text=(paste0("mutMat=cbind(mutMat,\"S:",s,"\"=sstatus)"))))
          latent=c(latent,"somatic")
          sval=c(sval,paste0("S:",s))
        }
      }
      
      #somatic BRAF_V600E
      if(cancer=="COAD" | cancer=="READ"){
        ssample=somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol=="BRAF" & somatic_likelyfunctional_driver$HGVSp_Short == "p.V600E"]
        scarrier=intersect(samples[samples%in%ssample],subimmune$bcr_patient_barcode)
        sstatus=ifelse(samples%in%ssample,1,0)
        if(length(scarrier)>3){
          mutMat=cbind(mutMat,"S:BRAF_V600E"=sstatus)
          latent=c(latent,"somatic")
          sval=c(sval,"S:BRAF_V600E")
        }
      }
      
      #add feature: mutation burden, MSIsensor 
      for(f in  c("TMB","MSISensor","MSIStatus","SNV_Neoantigens","Indel_Neoantigens","Homologous_Recombination_Defects")){
        latent=c(latent,f)
        subimmune=immuneProfile[immuneProfile$variable==f,]
        eval(parse(text=paste0("mutMat=cbind(mutMat,",f,"=subimmune$value[pmatch(samples,subimmune$bcr_patient_barcode)])")))
      }
      
      #add sample and tissue type
      mutMat=cbind(mutMat,cancer=clin$type[pmatch(samples,clin$bcr_patient_barcode)])
      colnames(mutMat)[1]="sample"
      mutMat=as.data.frame(mutMat)
      
      mutMat$TMB=log2(as.numeric(as.matrix(mutMat$TMB))+1)
      mutMat$SNV_Neoantigens=log2(as.numeric(as.matrix(mutMat$SNV_Neoantigens))+1)
      mutMat$Indel_Neoantigens=log2(as.numeric(as.matrix(mutMat$Indel_Neoantigens))+1)
 
        subMat=mutMat[which(mutMat$cancer==cancer),]
        if(length(table(subMat[,pheno]))==0){
          next
        }
        
        if(is.null(gval) & is.null(sval)){
          next
        }else if(is.null(gval)){
          somatic = c(0, 0)
          feature = c(1, 0)
          rus_path = rbind(somatic, feature)
          # add optional column names
          colnames(rus_path) = rownames(rus_path)
          
          rus_blocks = list(sval,pheno)
          rus_modes = c("B","A")
          
          subMat[,pheno]=as.numeric(as.matrix(subMat[,pheno]))
          if(length(which(is.na(subMat[,pheno])))>0){
            subMat=subMat[-which(is.na(subMat[,pheno])),]
          }
          rus_pls = plspm(subMat, rus_path, rus_blocks, modes = rus_modes)
    
          tmp=rus_pls_sum(rus_pls,factor1="germline",factor2="somatic")
          mainContribution=rbind(mainContribution,c(pheno,cancer,tmp$main))
          subContribution=rbind(subContribution,cbind(pheno,cancer,tmp$loading))
          
        }else if(is.null(sval)){
          germline = c(0, 0)
          feature = c(1, 0)
          rus_path = rbind(germline, feature)
          # add optional column names
          colnames(rus_path) = rownames(rus_path)
          
          subMat[,pheno]=as.numeric(as.matrix(subMat[,pheno]))
          if(length(which(is.na(subMat[,pheno])))>0){
            subMat=subMat[-which(is.na(subMat[,pheno])),]
          }
          
          rus_blocks = list(gval,pheno)
          rus_modes = c("B","A")
          rus_pls = plspm(subMat, rus_path, rus_blocks, modes = rus_modes)
          
          tmp=rus_pls_sum(rus_pls,factor1="germline",factor2="somatic")
          mainContribution=rbind(mainContribution,c(pheno,cancer,tmp$main))
          subContribution=rbind(subContribution,cbind(pheno,cancer,tmp$loading))
        }else{
          germline = c(0, 0, 0)
          somatic = c(0, 0, 0)
          feature = c(1, 1, 0)
          rus_path = rbind(germline, somatic, feature)
          # add optional column names
          colnames(rus_path) = rownames(rus_path)
          
          rus_blocks = list(gval,sval,pheno)
          rus_modes = c("B","B","A")
          
          subMat[,pheno]=as.numeric(as.matrix(subMat[,pheno]))
          if(length(which(is.na(subMat[,pheno])))>0){
          subMat=subMat[-which(is.na(subMat[,pheno])),]
          }
          
          rus_pls = plspm(subMat, rus_path, rus_blocks, modes = rus_modes)
          tmp=rus_pls_sum(rus_pls,factor1="germline",factor2="somatic")
          mainContribution=rbind(mainContribution,c(pheno,cancer,tmp$main))
          subContribution=rbind(subContribution,cbind(pheno,cancer,tmp$loading))
        }
        print(paste0(pheno,":",cancer))
    }
}



#germline,somatic->feature
mainContribution=as.data.frame(mainContribution)
colnames(mainContribution)=c("Phenotype","Cancer","germline_P","somatic_P","germline_Contribution","somatic_Contribution")
write.table(mainContribution,"./out/PLS-PM_germline_somatic_contribution_MMR_MSI_genes_20200204.tsv",quote = F,row.names =F,sep="\t")


#gene->germline, somatic
subContribution=cbind(rownames(subContribution),subContribution)
colnames(subContribution)=c("genes","phenotype","cancer","contribution")
write.table(subContribution,"./out/PLS-PM_genelevel_contribution_MMR_MSI_genes_20200204.tsv",quote = F,row.names =F,sep="\t")

###############################################################
#####################2. Visulization##########################
##############################################################
mainContribution<-fread("./out/PLS-PM_germline_somatic_contribution_MMR_MSI_genes_20191117.tsv")
subContribution<-fread("./out/PLS-PM_genelevel_contribution_MMR_MSI_genes_20191117.tsv")

subContribution=as.data.frame(subContribution)

finalMat=NULL
for(pheno in c("TMB","MSISensor","SNV_Neoantigens","Indel_Neoantigens","Homologous_Recombination_Defects")){
       contribution=as.data.frame(mainContribution[mainContribution$Phenotype==pheno,])
       contribution$germline_P=as.numeric(as.matrix(contribution$germline_P))
       contribution$germline_P=ifelse(contribution$germline_P==0,1,contribution$germline_P)
      
       contribution$somatic_P=as.numeric(as.matrix(contribution$somatic_P))
       contribution$somatic_P=ifelse(contribution$somatic_P==0,1,contribution$somatic_P)
      
       contribution$germline_Contribution=as.numeric(as.matrix(contribution$germline_Contribution))
       contribution$somatic_Contribution=as.numeric(as.matrix(contribution$somatic_Contribution))
  
      #germline FDR 
       contribution$G_FDR=p.adjust(contribution$germline_P)
       contribution$G_Sig=ifelse(contribution$G_FDR<0.15,"Suggestive (FDR<0.15)","None")
       contribution$G_Sig=ifelse(contribution$G_FDR<0.05,"Significant (FDR<0.05)",contribution$G_Sig)
        #somatic FDR 
       contribution$S_FDR=p.adjust(contribution$somatic_P)
       contribution$S_Sig=ifelse(contribution$S_FDR<0.15,"Suggestive (FDR<0.15)","None")
       contribution$S_Sig=ifelse(contribution$S_FDR<0.05,"Significant (FDR<0.05)",contribution$S_Sig)
       
      g=contribution[,c("Cancer","germline_P","germline_Contribution","somatic_Contribution","G_FDR","G_Sig")]
      s=contribution[,c("Cancer","somatic_P","germline_Contribution","somatic_Contribution","S_FDR","S_Sig")]
      colnames(g)=colnames(s)=c("Cancer","P","germline_Contribution","somatic_Contribution","FDR","Sig")
      
      plotMat=cbind(rbind(g,s),c(rep(c("Germline","Somatic"),each=dim(contribution)[1])))
      colnames(plotMat)[7]="Type"
      #plotMat$Sig=gsub("Germline ","",plotMat$Sig)
      #plotMat$Sig=gsub("Somatic ","",plotMat$Sig)
      plotMat$Sig=factor(plotMat$Sig,levels=c("Significant (FDR<0.05)","Suggestive (FDR<0.15)","None"))
     #plotMat$Cancer=c(as.character(plotMat$Cancer[1:dim(contribution)[1]]),rep("",dim(contribution)[1]))
      plotMat=cbind(plotMat,"Pheno"=pheno)
      finalMat=rbind(finalMat,plotMat)
}

#tmp=finalMat[finalMat$Cancer%in%c("COAD","UCEC")&finalMat$Pheno%in%c("SNV_Neoantigens","Indel_Neoantigens","TMB"),]
#tmp=finalMat[finalMat$Cancer%in%c("COAD","UCEC","LUSC","STAD")&finalMat$Pheno%in%c("MSISensor"),]
#View(tmp)


th= theme_bw()+ theme(legend.key.size =unit(.2, "cm"),legend.title = element_text(size=8) ,panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=14,hjust = 0.95), axis.text.y = element_text(colour="black", size=14,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=14,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.5, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

finalMat$Pheno=gsub("Homologous_Recombination_Defects","HRD",finalMat$Pheno)
finalMat$Pheno=factor(finalMat$Pheno,levels=c("HRD","TMB","SNV_Neoantigens","Indel_Neoantigens","MSISensor"))

finalMat$Cancer[finalMat$Type=="Somatic"]=""#keep only one label
p = ggplot(finalMat, aes(y=somatic_Contribution, x=germline_Contribution)) 
p = p + geom_point(aes(colour=factor(finalMat$Sig),size=-log10(finalMat$P),shape =Type) ,alpha = I(0.6))+scale_shape_manual(values=c(1, 5))#+xlim(-0.22,0.6)+ylim(-0.22,0.6)
p = p + geom_text_repel(aes(label=Cancer),size=3,segment.alpha =0.5,segment.colour ="grey")
p = p+facet_grid(.~Pheno,drop=T)#, space="free",scale="free"
p = p + ylab("Somatic coefficients") + xlab("Germline coefficients") + theme_bw() 
p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) 
p = p + guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA)))+ theme(plot.title = element_text(size=14, face="bold"),axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14)) + labs(size = "-log10(P)") +scale_size_continuous(limits=c(0,50),breaks = c(1,10,20,30,40,50))
p = p + scale_colour_manual("", values = c("Significant (FDR<0.05)" = "red",  "Suggestive (FDR<0.15)" = "blue","None" = "grey"))
p=p+th
p
fn = paste0("out/immune_germline_somatic_contribution_20190920.pdf")
ggsave(fn,w = 12, h = 3, useDingbat=F,limitsize = FALSE)


#
#draw germline somatic contribution by cancer types

library(diagram)

mainContribution$Phenotype=gsub("Homologous_Recombination_Defects","HRD",mainContribution$Phenotype)
subContribution$phenotype=gsub("Homologous_Recombination_Defects","HRD",subContribution$phenotype)


for(pheno in  c("TMB","MSISensor","SNV_Neoantigens","Indel_Neoantigens","HRD")){
#pheno="Homologous_Recombination_Defects"
mainContribution1=as.data.frame(mainContribution[mainContribution$Phenotype==pheno,])
subContribution1=as.data.frame(subContribution[subContribution[,2]==pheno,])


cancers=c("BRCA","OV","PAAD","UCEC","READ","COAD","STAD","BLCA","LUAD","LUSC","HNSC")
for(cancer in cancers){
    mainContribution2=mainContribution1[mainContribution1$Cancer==cancer,]
    mainContribution2$germline_Contribution=as.numeric(as.matrix(mainContribution2$germline_Contribution))
    mainContribution2$somatic_Contribution=as.numeric(as.matrix(mainContribution2$somatic_Contribution))
    subContribution2=subContribution1[subContribution1$cancer==cancer,]
    subContribution2$contribution=as.numeric(as.matrix(subContribution2$contribution))
    
    if(dim(mainContribution2)[1]==0){
      next
    }
    
    g=grep("G:",subContribution2$genes,value=T)
    s=grep("S:",subContribution2$genes,value=T)
    ng=length(g)#number of genes involved
    ns=length(s)
    pdf(paste0("./out/",cancer,"_",pheno,"_plsms_contributionfigure_20190920.pdf"),h=5,w=5.5)
    par(mar = c(1, 1, 2, 1))
    
    openplotmat(cex=2)
    title(cancer, adj=0)
    #feature nodes coordinates
    #feature
    elpos_f=c(0.85,0.5)
    #latent variable
    if(mainContribution2$germline_Contribution!=0 & mainContribution2$somatic_Contribution!=0){
      elpos_l=rbind(c(0.5,0.75),c(0.5,0.25))
      if(ng==1){
        elpos_g=as.matrix(cbind(0.1,0.95))
      }else{
        elpos_g <-cbind(0.1,seq(from=0.54,to=0.96,by=(0.96-0.55)/(ng-1)))
        elpos_g <-elpos_g[order(elpos_g[,2],decreasing = T),]
      }
      
      if(ns==1){
        elpos_s=as.matrix(cbind(0.1,0.25))
      }else{
        elpos_s <-cbind(0.1,seq(from=0.04,to=0.46,by=(0.46-0.04)/(ns-1)))
        elpos_s <-elpos_s[order(elpos_s[,2],decreasing = T),]
      }
      #draw connection: gene-> latent variable
      
      #the position of gene coefficient
      for ( i in 1:ng){
        arrpos_g<-straightarrow(from = elpos_g[i, ], to = elpos_l[1, ], lty = 1, lcol = "grey")
        text(arrpos_g[1] - 0, arrpos_g[2],cex=1, round(as.numeric(subContribution2[subContribution2$genes==g[i],"contribution"]),digits=3))
      }
      
      for ( i in 1:ns){
        arrpos_s<-straightarrow(from = elpos_s[i, ], to = elpos_l[2, ], lty = 1, lcol = "grey")
        text(arrpos_s[1] + 0, arrpos_s[2],cex=1, round(as.numeric(subContribution2[subContribution2$genes==s[i],"contribution"]),digits=3))
      }
      
      #draw connection: latent variable->feature
      #the position of latent variable coefficient
      arrpos_g_f=straightarrow(from = elpos_l[1,], to = elpos_f, lty = 1, lcol = "grey")
      text(arrpos_g_f[1] + 0, arrpos_g_f[2],cex=1, round(as.numeric(mainContribution2[,"germline_Contribution"]),digits=3))
      arrpos_s_f=straightarrow(from = elpos_l[2,], to = elpos_f, lty = 1, lcol = "grey")
      text(arrpos_s_f[1] + 0, arrpos_s_f[2],cex=1, round(as.numeric(mainContribution2[,"somatic_Contribution"]),digits=3))
      
      #draw labels
      textrect (elpos_f,0.05,0.02, lab = pheno, cex = 1, box.col = "green",lcol ="green")
      textrect (elpos_l[1,],0.05,0.02, lab = "germline", cex = 1.5, box.col = "grey",lcol ="grey")
      textrect (elpos_l[2,],0.05,0.02, lab = "somatic", cex = 1.5, box.col = "grey",lcol ="grey")
      
      for ( i in 1:ng)
        textrect (elpos_g[i,],0.05,0.02, lab = gsub("G:","",g[i]), cex = 1.5,box.col = "lightblue",lcol ="lightblue",font=3)
      for ( i in 1:ns)
        textrect (elpos_s[i,],0.05,0.02, lab = gsub("S:","",s[i]), cex = 1.5,box.col = "orange",lcol ="orange",font=3)

      
    }else if(mainContribution2$germline_Contribution==0){
      elpos_l=c(0.5,0.5)
      if(ns==1){
        elpos=c(0.2,0.5)
        arrpos_s<-straightarrow(from = elpos, to = elpos_l, lty = 1, lcol = "grey")
        text(arrpos_s[1] + 0, arrpos_s[2],cex=1, round(as.numeric(subContribution2[subContribution2$genes==g[i],"contribution"]),digits=3))
      }else{
        elpos <-cbind(0.2,seq(from=0.15,to=0.85,by=(0.85-0.15)/ns))
        for ( i in 1:ns){
          arrpos_s<-straightarrow(from = elpos[i, ], to = elpos_l, lty = 1, lcol = "grey")
          text(arrpos_s[1] + 0, arrpos_s[2],cex=1, round(as.numeric(subContribution2[subContribution2$genes==s[i],"contribution"]),digits=3))
        }
        #draw connection: latent variable->feature
        arrpos_s_f=straightarrow(from = elpos_l, to = elpos_f, lty = 1, lcol = "grey")
        text(arrpos_s_f[1] - 0, arrpos_s_f[2],cex=1, round(as.numeric(mainContribution2[,"somatic_Contribution"]),digits=3))
        #variable label
        textrect (elpos_f,0.05,0.02, lab = pheno, cex = 1, box.col = "green",lcol ="green")
        textrect (elpos_l,0.05,0.02, lab = "somatic", cex = 1.5, box.col = "grey",lcol ="grey")
        for ( i in 1:ns)
          textrect (elpos[i,],0.05,0.02, lab = gsub("S:","",s[i]), cex = 1.5,box.col = "orange",lcol ="orange",font=3)
      }
    }
    dev.off()
}
}


##################3. Figure 4D-I Distribution of germline somatic carrier########################
#figure parameters
th= theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='none',axis.text.x = element_text(colour="black",angle = 90, size=12), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=12,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))


library(ggsignif)

#HRD
cancer=c("BRCA","OV","PAAD","STAD")
pheno="Homologous_Recombination_Defects"
genes=genes_HRD#c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","TGFBR2")
col_value=c("both"="#56B4E9","germline"="#56B4E9","somatic"="#56B4E9","NA"="#dbd9d9")
#label="Figure5_";short_name="HRD";gene=genes;
immuneFeature<-function(label,cancer,pheno,labpos,short_name,gene,col_value,th,w,h){
  #label="Figure5_";cancer=cancer;pheno=pheno;labpos=500;short_name="TMB";gene=genes;col_value=col_value;th=th;w=6;h=3
    subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]
  
  mut_status=sapply(sample,function(x){
    g=pathVarP$HUGO_Symbol[pathVarP$bcr_patient_barcode==x]
    s=somatic_likelyfunctional_driver$Hugo_Symbol[somatic_likelyfunctional_driver$bcr_patient_barcode==x]
    if(any(genes%in%g) & any(genes%in%s)){
      return("both")
    }else if(any(genes%in%g) & !any(genes%in%s)){
      return("germline")
    }else if(!any(genes%in%g) & any(genes%in%s)){
      return("somatic")
    }else{
      return("NA")
    }
  })
  
  mutMat=as.data.frame(cbind(sample,pheno=subimmnue$value[pmatch(sample,subimmnue$bcr_patient_barcode)],cancer=cancers,mut=mut_status))
  
  if(any(is.na(mutMat$pheno))){
    mutMat=mutMat[-which(is.na(mutMat$pheno)),]
  }
  
  mutMat$mut=factor(mutMat$mut,levels=c("both","germline","somatic","NA"))
  mutMat$pheno=as.numeric(as.matrix(mutMat$pheno))
  
  mutype=as.character(unique(mutMat$mut))
  mutype=mutype[-which(mutype=="NA")]
  
  ann_text=NULL
  pvalueMat=NULL
  for(ca in unique(mutMat$cancer)){
    subMat=mutMat[mutMat$cancer==ca,]
    for(m in mutype){
      if(!any(subMat$mut==m)){next}
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
      ann_text=rbind(ann_text,t(cbind(c(cancer=ca,mut=m,pheno=pheno,lab=plab))))
      pvalueMat=rbind(pvalueMat,t(cbind(c(cancer=ca,mut=m,pheno=pheno,lab=as.numeric(tmpPvalue)))))
    }
  }
  ann_text=as.data.frame(ann_text)
  pvalueMat=as.data.frame(pvalueMat)
  pvalueMat$lab=as.numeric(as.matrix(pvalueMat$lab))
  tmp=pvalueMat[pvalueMat$cancer%in%c("BRCA","UCEC") & pvalueMat$pheno%in%c("SNV_Neoantigens","Indel_Neoantigens"),]
  
  p=ggplot(mutMat, aes(x=mut, y=as.numeric(as.matrix(pheno)),color=I(mut),fill=I(mut)))+geom_violin()+geom_boxplot(width=0.1,color="grey",outlier.size = 0.2)
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+facet_grid(.~cancer, space="free",scale="free",drop=T)
  p=p+scale_fill_manual(values=col_value)+scale_color_manual(values=col_value)
  p=p+xlab("")+ylab(paste0(short_name))+scale_y_log10(labels = function(x) format(x, scientific = FALSE))
  p=p+geom_text(data = ann_text,label=ann_text$lab)
  p=p+th
  print(p)
  fn = paste0("out/",label,"_",pheno,"_MMRgenes_20190830.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}

immuneFeature("testFigure5_",cancer=cancer,pheno=pheno,labpos=100,short_name="HRD",gene=genes,col_value=col_value,th=th,w=6,h=3)


#TMB
cancer=c("BRCA","OV","PAAD","STAD","UCEC")
pheno="TMB"
genes=genes_TMB#c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","TGFBR2")
col_value=c("both"="#56B4E9","germline"="#56B4E9","somatic"="#56B4E9","NA"="#dbd9d9")
immuneFeature("testFigure5_",cancer=cancer,pheno=pheno,labpos=500,short_name="TMB",gene=genes,col_value=col_value,th=th,w=6,h=3)



#SNV_Neoantigens
cancer=c("BRCA","UCEC")
pheno="SNV_Neoantigens"
genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","TGFBR2")
col_value=c("both"="#56B4E9","germline"="#56B4E9","somatic"="#56B4E9","NA"="#dbd9d9")
immuneFeature("Figure5_",cancer=cancer,pheno=pheno,labpos=12000,short_name="SNV Neoantigens",gene=genes,col_value=col_value,th=th,w=3,h=3)


#Indel_Neoantigens
cancer=c("BRCA","UCEC")
pheno="Indel_Neoantigens"
genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","TGFBR2")
col_value=c("both"="#56B4E9","germline"="#56B4E9","somatic"="#56B4E9","NA"="#dbd9d9")

immuneFeature("Figure5_",cancer=cancer,pheno=pheno,labpos=10000,short_name="Indel Neoantigens",gene=genes,col_value=col_value,th=th,w=3,h=3)




#MSISensor
cancer=c("STAD","UCEC","COAD","LUSC","SKCM")
pheno="MSISensor"
genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","TGFBR2")
col_value=c("both"="#56B4E9","germline"="#56B4E9","somatic"="#56B4E9","NA"="#dbd9d9")
immuneFeature<-function(label,cancer,pheno,short_name,gene,col_value,th,w,h){
  subimmnue=immuneProfile[(immuneProfile$TCGA_Study%in%cancer)&(immuneProfile$variable==pheno),]
  sample=clin$bcr_patient_barcode[clin$type%in%cancer]
  cancers=clin$type[clin$type%in%cancer]
  
  mut_status=sapply(sample,function(x){
    g=pathVarP$HUGO_Symbol[pathVarP$bcr_patient_barcode==x]
    s=somatic_likelyfunctional_driver$Hugo_Symbol[somatic_likelyfunctional_driver$bcr_patient_barcode==x]
    if(any(genes%in%g) & any(genes%in%s)){
      return("both")
    }else if(any(genes%in%g) & !any(genes%in%s)){
      return("germline")
    }else if(!any(genes%in%g) & any(genes%in%s)){
      return("somatic")
    }else{
      return("NA")
    }
  })
  
  mutMat=as.data.frame(cbind(sample,pheno=subimmnue$value[pmatch(sample,subimmnue$bcr_patient_barcode)],cancer=cancers,mut=mut_status))
  
  mutMat$mut=factor(mutMat$mut,levels=c("both","germline","somatic","NA"))
  
  
  mutMat=mutMat[-which(is.na(mutMat$pheno)),]
  mutMat$pheno=as.numeric(as.matrix(mutMat$pheno))
  p=ggplot(mutMat, aes(x=mut, y=log2(as.numeric(as.matrix(pheno+0.01))),color=I(mut),fill=I(mut)))+geom_violin()+geom_boxplot(width=0.1,color="grey",outlier.size = 0.2)
  p=p+geom_jitter(height = 0, width = 0.1,size=0.4,alpha=0.2,color="black")
  p=p+facet_grid(.~cancer, space="free",scale="free",drop=T)
  p=p+scale_fill_manual(values=col_value)+scale_color_manual(values=col_value)
  p=p+xlab("")+ylab(paste0(short_name))
  p=p+th
  #print(p)
  fn = paste0("out/",label,"_",pheno,"_MMRgenes_20190830.pdf")
  ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}

immuneFeature("Figure5_",cancer=cancer,pheno=pheno,short_name="log2(MSISensor+0.01)",gene=genes,col_value=col_value,th=th,w=6,h=3)

