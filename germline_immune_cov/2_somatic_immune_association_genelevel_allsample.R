#germline variant clinical score association
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/")
source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../stat_functions.R")
source("../load_somatic.R")

#remove hyper mutator: TRUE of FALSE
rmhypermutator=FALSE

### MAIN ###
#immune response features for consideration
response=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)

samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_n10080.txt",h=F,stringsAsFactors = F)[,1]

label="overlaped_n10080"

#immune score
immuneProfile=read.table("./out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]
immuneProfile=immuneProfile[which(immuneProfile$variable%in%response[,2]),]
immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]

#somatic
somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[which(somatic_likelyfunctional_driver$bcr_patient_barcode%in%samples),]

# read input files
gene_sample = data.frame(table(somatic_likelyfunctional_driver$Hugo_Symbol,somatic_likelyfunctional_driver$bcr_patient_barcode))
colnames(gene_sample) = c("Gene","bcr_patient_barcode","Freq")
gene_sample=gene_sample[gene_sample$Freq!=0,]

##### individual cancer type analysis #####
cancers = unique(clin$type)
genes = unique(somatic_likelyfunctional_driver$Hugo_Symbol)
#save.image(paste0("2_somaticDriver_immune_association_genelevel",label,".RData"))


#this is a time consuming step, so I move the analysis to a super computer cluster for parallel computating

#setwd("/gpfs/ysm/project/tq37/GermlineSomatic/analysis/germline_immune_cov")
setwd("/home/tq37/project/GermlineSomatic/analysis/germline_immune_cov")
load("2_somaticDriver_immune_association_geneleveloverlaped_n10080.RData")

library(data.table)
core=20
require(foreach)
require(doParallel)
registerDoParallel(core)

#response=read.table("./out/allAssoc_iterms_plot 2 neo normal",sep="\t",h=F)

finalMatrix <- foreach(gene=genes, .combine=rbind) %dopar% {
  tt=NULL
  gene_sample_g = gene_sample[gene_sample$Gene==gene,]
  var_exp_g = merge(immuneProfile,gene_sample_g,by="bcr_patient_barcode",all.x=T)
  var_exp_g$Freq[is.na(var_exp_g$Freq)] = 0
  var_exp_g$Freq[var_exp_g$Freq != 0 ] = 1
  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$TCGA_Study %in% cancer,]
    for (clinicPheno in unique(var_exp_g_c$variable)){
      var_exp_g_c_s = var_exp_g_c[var_exp_g_c$variable == clinicPheno,]
      gene_path_count=sum(var_exp_g_c_s$Freq)
      if((dim(var_exp_g_c_s)[1]>10) && sum(var_exp_g_c_s$Freq) > 2){
        pheno=as.character(as.matrix(response[response[,2]==clinicPheno,]))
        df=var_exp_g_c_s
        if(all(is.na(df$value)) | all(is.na(var_exp_g_c_s$value[var_exp_g_c_s$Freq!=0])) | all(is.na(var_exp_g_c_s$value[var_exp_g_c_s$Freq==0]))){
          next
        }
        
        # run GLM
        w = wilcox.test(var_exp_g_c_s$value[var_exp_g_c_s$Freq==0],var_exp_g_c_s$value[var_exp_g_c_s$Freq!=0])
        wP = w$p.value
        wWstat = w$statistic
        #cancer_gene_stat = run_glm(var_exp_g_c_s,yi="value",xi="Freq",covi=c("Age","Sex","Race"),ytype="Continuous")
        
        if(pheno[4]=="log2"){
          df$value=log2(df$value+as.numeric(pheno[5]))
        }
        
        
        if(clinicPheno=="TMB"){
          #pheno[3]="Poisson"
          if(length(which(df$value==0))>0){
            df=df[df$value!=0,]
          }
          df$value=log2(as.numeric(as.matrix(df$value)))
          #df$value=log2(as.numeric(as.matrix(df$value))+0.001)
        }
        
        #if(){
        #cancer_gene_stat=run_glm(df,yi="value",xi="Freq",covi=c("Age","PC1","PC2"),ytype=pheno[3],gene=gene,cancer=cancer)
        #}else{
          cancer_gene_stat=run_glm(df,yi="value",xi="Freq",covi=c("Age","PC1","PC2"),ytype=pheno[3],gene=gene,cancer=cancer)
        #}
        #,"TCGA_Subtype"
        if(pheno[3]=="Poisson" | pheno[3]=="Binary"){
          cancer_gene_stat=cancer_gene_stat[-which(colnames(cancer_gene_stat)%in%c("F","z.value","xi_lvl1"))]
        }else{
          cancer_gene_stat=cancer_gene_stat[-which(colnames(cancer_gene_stat)%in%c("F","t.value","xi_lvl1"))]
        }
        colnames(cancer_gene_stat) = c("yi","ytype","xi","df","Deviance","Resid. Df","Resid. Dev","Pr(>F)","Estimate","Std..Error","Pr...t..","covars");
        # compile results
        #full_cancer_gene_stat = cbind(cancer,gene,clinicPheno,gene_path_count,wP,wWstat,cancer_gene_stat)
        full_cancer_gene_stat = cbind(cancer,gene,clinicPheno,gene_path_count,wP,wWstat,cancer_gene_stat)
        tt = rbind(tt, full_cancer_gene_stat)
      }
    }
  }
  print(grep(gene,genes))
  return(tt)
}


#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
#tt=finalMatrix
colnames(finalMatrix) = c("cancer","gene","clinicPhenotype","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","p-value(Anova)","coefficient","StdError","p-value","covariants");

#colnames(tt) = c("cancer","gene","clinicPhenotype","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","F_statistic","p-value(Anova)","coefficient","StdError","p-value","covariants");
#tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
#tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")

finalMatrix=finalMatrix[order(finalMatrix$`p-value`, decreasing=FALSE),]
tn = paste0("./out/somaticDriver_ImmuneAssoc_genelevel",label,"byCancerType_cov_20200605.txt")
#tn = paste0("./out/somaticDriver_ImmuneAssoc_genelevel",label,"byCancerType_cov_20200409_neo_normal.txt")

write.table(finalMatrix, quote=F, sep="\t", file = tn, row.names = F)



