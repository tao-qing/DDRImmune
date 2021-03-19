setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/DDRNeo_DDRImmune_Assoc/")

library(data.table)
library(dplyr)

source("../dependecy_file_immune.R")

#enrichment level in mutant vs wildtype
phenoMat=NULL
for(ge in ddrgene){
  for(ca in unique(clin$type)){
    allspl=intersect(clin$bcr_patient_barcode[clin$type==ca],samples)
    for(ph in immunopheno){
      gmutspl=intersect(unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol==ge]),allspl)
      gmutmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%gmutspl])
      
      gwtspl=allspl[!allspl%in%gmutspl]
      gwtmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%gwtspl])
      
      smutspl=intersect(unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol==ge]),allspl)
      smutmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%smutspl])
      
      swtspl=allspl[!allspl%in%smutspl]
      swtmean=mean(immuneProfile$value[immuneProfile$variable==ph & immuneProfile$bcr_patient_barcode%in%swtspl])
      
      gratio=(length(gmutspl)/length(allspl))*100
      sratio=(length(smutspl)/length(allspl))*100
      
      results=rbind(Gene=cbind(paste0("g",ge),Cancer=ca,Pheno=ph,Ratio=gratio,MUT=gmutmean,WT=gwtmean),cbind(Gene=paste0("s",ge),Cancer=ca,Pheno=ph,Ratio=sratio,MUT=smutmean,WT=swtmean))
      phenoMat=rbind(phenoMat,results)
    }
  }
}

phenoMat=as.data.frame(phenoMat)
fwrite(phenoMat,"./out/mean_immune_level_mutvswt.txt")

phenoMat<-fread("./out/mean_immune_level_mutvswt.txt",data.table=F,stringsAsFactors = F)
#germline 
g<-fread("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/pathVarP_ImmuneAssoc_geneleveloverlaped_n9738byCancerType_cov_20200409.txt",data.table=F)%>%filter(gene%in%ddrgene)%>%filter(clinicPhenotype%in%immunopheno)%>%rename('pvalue'='p-value')%>%mutate(FDR=p.adjust(pvalue,method="fdr"))%>%mutate(Sig=ifelse(FDR<0.05,"Significant (FDR<0.05)","None"))%>%select(c("cancer","gene","gene_path_count","coefficient","pvalue","FDR","Sig","clinicPhenotype"))%>%rename("Cancer"="cancer","Gene"="gene","#samples with mutations"="gene_path_count","Coefficient"="coefficient","P"="pvalue","Pheno"="clinicPhenotype")%>%mutate(Gene=paste0("g",Gene))

s<-fread("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/somaticDriver_ImmuneAssoc_geneleveloverlaped_n9738byCancerType_cov_20200409.txt",data.table=F)%>%filter(gene%in%ddrgene)%>%filter(clinicPhenotype%in%immunopheno)%>%rename('pvalue'='p-value')%>%mutate(FDR=p.adjust(pvalue,method="fdr"))%>%mutate(Sig=ifelse(FDR<0.05,"Significant (FDR<0.05)","None"))%>%select(c("cancer","gene","gene_path_count","coefficient","pvalue","FDR","Sig","clinicPhenotype"))%>%rename("Cancer"="cancer","Gene"="gene","#samples with mutations"="gene_path_count","Coefficient"="coefficient","P"="pvalue","Pheno"="clinicPhenotype")%>%mutate(Gene=paste0("s",Gene))


m=rbind(g,s)


phenoMat=phenoMat%>%rename(Gene=V1)%>%mutate(label=paste0(Gene,"|",Cancer,"|",Pheno))

m1=m%>%mutate(label=paste0(Gene,"|",Cancer,"|",Pheno))%>%left_join(.,phenoMat[,c("label","Ratio","MUT","WT")],by="label")%>%select(c("Cancer","Gene","Ratio","Pheno","MUT","WT","Coefficient","P","FDR","Sig"))%>%mutate(Pheno=gsub("Lymphocyte_Infiltration_Signature_Score","TILs",gsub("PDCD1","PD1",gsub("CD274","PD-L1",Pheno))))



fwrite(m1,"DDR_immune_signature_association_genelevel_TableS4.csv")








