setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/DDRNeo_DDRImmune_Assoc/")

library(ggplot2)
library(ggsignif)
library(ggrepel)
library(dplyr)

source("../dependecy_file_immune.R")

cancer = unique(clin$type)
genes = list()
genes=NULL
genes["HR"]=list(c("BRCA1","BRCA2","PALB2"))
genes["Polymerase"]=list(c("POLE","POLQ"))
genes["Sensor"]=list(c("ATM","ATR","CHEK2"))
genes["MMR"]=list(c("MLH1","MSH2","MSH3","MSH6","PMS2"))

#germline 
g<-fread("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/germline_ImmuneAssoc_pathwayleveloverlaped_n9738byCancerType_cov_20200603.txt",data.table=F)%>%filter(gene%in%names(genes) & clinicPhenotype%in%immunopheno & gene_path_count >=4 )%>%rename('pvalue'='p-value')%>%mutate(FDR=p.adjust(pvalue,method="fdr"))%>%mutate(Sig=ifelse(FDR<0.05,"Significant (FDR<0.05)","None"))%>%select(c("cancer","gene","gene_path_count","coefficient","pvalue","FDR","Sig","clinicPhenotype"))%>%rename("Cancer"="cancer","Gene"="gene","#samples with mutations"="gene_path_count","Coefficient"="coefficient","P"="pvalue","Pheno"="clinicPhenotype")%>%mutate(Gene=paste0("g",Gene))

s<-fread("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/somaticDriver_ImmuneAssoc_pathwayleveloverlaped_n9738byCancerType_cov_20200603.txt",data.table=F)%>%filter(gene%in%names(genes) & clinicPhenotype%in%immunopheno & gene_path_count >=4 )%>%rename('pvalue'='p-value')%>%mutate(FDR=p.adjust(pvalue,method="fdr"))%>%mutate(Sig=ifelse(FDR<0.05,"Significant (FDR<0.05)","None"))%>%select(c("cancer","gene","gene_path_count","coefficient","pvalue","FDR","Sig","clinicPhenotype"))%>%rename("Cancer"="cancer","Gene"="gene","#samples with mutations"="gene_path_count","Coefficient"="coefficient","P"="pvalue","Pheno"="clinicPhenotype")%>%mutate(Gene=paste0("s",Gene))

m=rbind(g,s)

#enrichment level in mutant vs wildtype
phenoMat=NULL 
for(ge in names(genes)){
  for(ca in unique(clin$type)){
    allspl=intersect(clin$bcr_patient_barcode[clin$type==ca],samples)
    for(ph in immunopheno){
      gmutspl=intersect(unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%genes[[ge]]]),allspl)
      gmutmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%gmutspl])
      
      gwtspl=allspl[!allspl%in%gmutspl]
      gwtmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%gwtspl])
      
      smutspl=intersect(unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol%in%genes[[ge]]]),allspl)
      smutmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%smutspl])
      
      swtspl=allspl[!allspl%in%smutspl]
      swtmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%swtspl])
      
      gratio=(length(gmutspl)/length(allspl))*100
      sratio=(length(smutspl)/length(allspl))*100
      
      results=rbind(Gene=cbind(paste0("g",ge),Cancer=ca,Pheno=ph,Ratio=gratio,MUT=gmutmean,WT=gwtmean),cbind(Gene=paste0("s",ge),Cancer=ca,Pheno=ph,Ratio=sratio,MUT=smutmean,WT=swtmean))
      phenoMat=rbind(phenoMat,results)
    }
  }
  print(ge)
}


phenoMat=phenoMat%>%as.data.frame()%>%rename(Gene=V1)%>%mutate(label=paste0(Gene,"|",Cancer,"|",Pheno))

m1=m%>%mutate(label=paste0(Gene,"|",Cancer,"|",Pheno))%>%left_join(.,phenoMat[,c("label","Ratio","MUT","WT")],by="label")%>%select(c("Cancer","Gene","Ratio","Pheno","MUT","WT","Coefficient","P","FDR","Sig"))

m2=m1%>%mutate(Pheno=gsub("Lymphocyte_Infiltration_Signature_Score","TILs",gsub("PDCD1","PD1",gsub("CD274","PD-L1",Pheno))))
fwrite(m2,"DDR_immune_signature_association_Pathwaylevel_TableS5.csv")

m1=m1%>%mutate(label=paste0(Cancer,":",Gene),MUT=as.numeric(as.matrix(MUT)),WT=as.numeric(as.matrix(WT)),FDR=as.numeric(as.matrix(FDR)))


immuneFeature<-function(mat,label,pheno,short_name,lpos,hpos,col_value,th,w,h){
 p = ggplot(mat, aes(y=WT, x=MUT)) #
 p = p + geom_point(aes(color=factor(mat$Sig)),alpha=0.5,size=-log10(mat$FDR))+xlim(lpos,hpos)+ylim(lpos,hpos)#+scale_shape_manual(values=c(1, 5))#,shape =Type$ ,alpha = I(0.6)
 p = p + geom_text_repel(aes(label=ifelse(mat$FDR<0.05,mat$label,"")),size=4,segment.alpha =0.5,fontface="bold.italic",segment.colour ="black",min.segment.length = 0) 
 #p = p+facet_grid(.~Gene,drop=T)#, space="free",scale="free"
 p = p + ylab(paste0(short_name," in WT")) + xlab(paste0(short_name," in MUT")) + theme_bw()
 p = p + geom_abline(intercept = 0, slope=1, alpha=0.2)
 #p = p + geom_smooth(data=phenoMatSub,aes(x=MUT, y =WT),colour="red",method = "glm",alpha=0.1,size=1.5)
 p = p + guides(color=guide_legend(override.aes=list(fill=NA)),linetype=guide_legend(override.aes=list(fill=NA)))+ theme(plot.title = element_text(size=14, face="bold"),axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14)) + labs(size = "-log10(P)") +scale_size_continuous(limits=c(0,50),breaks = c(1,10,20,30,40,50))
 #p = p + scale_colour_manual("", values = c("Significant (FDR<0.05)" = "red",  "Suggestive (FDR<0.15)" = "blue","None" = "grey"))
 p = p + scale_color_manual("FDR<0.05", values = c("Significant (FDR<0.05)" = "red","None" = "grey"))
 #p = p +labs(title=paste0("\nslope=",slope), size = 16, colour = "red")#+xlab(shortName1)+ylab(shortName2)
# p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
# p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
 p=p+th
 p
 fn = paste0("out/",label,"_",pheno,"_DDRgeneMutvsWT_20201005.pdf")
 ggsave(fn,w = w, h = h, useDingbat=F,limitsize = FALSE)
}


th= theme_bw()+ theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=15,hjust=0.90), axis.text.y = element_text(colour="black", size=15,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=14,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))


m3=m1[grep("^s",m1$Gene,perl=T),]

pheno="PDCD1"
immuneFeature(mat=m3[m3$Pheno==pheno,],"Figure_",pheno=pheno,lpos=4,hpos=7.5,short_name="PD1",th=th,w=3,h=3)


pheno="Lymphocyte_Infiltration_Signature_Score"
immuneFeature(mat=m3[m3$Pheno==pheno,],"Figure_",pheno=pheno,lpos=-1,hpos=1,short_name="TILs",th=th,w=3,h=3)

pheno="CYTScore"
immuneFeature(mat=m3[m3$Pheno==pheno,],"Figure_",pheno=pheno,lpos=6,hpos=8.5,short_name="CYTScore",th=th,w=3,h=3)


pheno="CD274"
immuneFeature(mat=m3[m3$Pheno==pheno,],"Figure_",pheno=pheno,lpos=2.7,hpos=8,short_name="PD-L1",th=th,w=3.0,h=3)
