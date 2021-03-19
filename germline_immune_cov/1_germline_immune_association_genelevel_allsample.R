#germline variant clinical score association analysis

###################################Prepare Data####################################
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/")
source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../stat_functions.R")

#remove hyper mutator: TRUE of FALSE
rmhypermutator=FALSE

### MAIN ###
#immune response features for consideration
response=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms_plot 2",sep="\t",h=F)
samples=read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_n10080.txt",h=F,stringsAsFactors = F)[,1]

label="overlaped_n10080"

#immune score with covariates
immuneProfile=read.table("./out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
immuneProfile=immuneProfile[-which(is.na(immuneProfile$value)),]
immuneProfile=immuneProfile[which(immuneProfile$variable%in%response[,2]),]
immuneProfile=immuneProfile[which(immuneProfile$bcr_patient_barcode%in%samples),]

pathVarP=pathVarP[pathVarP$bcr_patient_barcode%in%samples,]

germlineGene=unique(pathVarP$HUGO_Symbol)

#if(rmhypermutator){
#  hyper<-read.table("../../Huang_lab_data/TCGA_PanCanAtlas_2018/somatic/Driver_BaileyCell2018/HypermutatorSamples.txt",h=T,sep="\t",stringsAsFactors=FALSE)[,1]
#  immuneProfile=immuneProfile[-which(immuneProfile$bcr_patient_barcode%in%hyper),]
#  pathVarP=pathVarP[-which(pathVarP$bcr_patient_barcode%in%hyper),]
#  label="_rmhypermutator_"
#}else{
#  label="_"
#}

# read input files
gene_sample = data.frame(table(pathVarP$HUGO_Symbol,pathVarP$bcr_patient_barcode))
colnames(gene_sample) = c("Gene","bcr_patient_barcode","Freq")
gene_sample=gene_sample[gene_sample$Freq!=0,]

##### individual cancer type analysis #####
cancers = unique(pathVarP$cancer)
genes = unique(pathVarP$HUGO_Symbol)
# limit runs to cancers with at least 5 likely patho/pathogenic variants
#save.image(paste0("2_germline_immune_association_genelevel",label,".RData"))


#module load R/3.5.0-foss-2016b-avx2
setwd("/home/tq37/project/GermlineSomatic/analysis/germline_immune_cov")
load("2_germline_immune_association_geneleveloverlaped_n10080.RData")

#response=read.table("./out/allAssoc_iterms_plot 2 neo normal",sep="\t",h=F)

library(data.table)
core=10
require(foreach)
require(doParallel)
registerDoParallel(core)

finalMatrix <- foreach(gene=genes, .combine=rbind) %dopar% {

tt=NULL
  gene_sample_g = gene_sample[gene_sample$Gene==gene,]
  
  #gene_sample_g_other=unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol%in%germlineGene[!germlineGene%in%gene]])
  
  var_exp_g = merge(immuneProfile,gene_sample_g,by="bcr_patient_barcode",all.x=T)
  var_exp_g$Freq[is.na(var_exp_g$Freq)] = 0
  var_exp_g$Freq[var_exp_g$Freq != 0 ] = 1
  
  #immuneProfile=immuneProfile[-which(immuneProfile$bcr_patient_barcode%in%gene_sample_g_other),]
  
  for (cancer in cancers){
    
    var_exp_g_c = var_exp_g[var_exp_g$TCGA_Study %in% cancer,]

      for (clinicPheno in unique(var_exp_g_c$variable)){
        var_exp_g_c_s = var_exp_g_c[var_exp_g_c$variable == clinicPheno,]
        gene_path_count=sum(var_exp_g_c_s$Freq)
        if((dim(var_exp_g_c_s)[1]>10) && sum(var_exp_g_c_s$Freq) > 2){

          pheno=as.character(as.matrix(response[response[,2]==clinicPheno,]))
          df=var_exp_g_c_s
          if(all(is.na(df$value)) | all(is.na(var_exp_g_c_s$value[var_exp_g_c_s$Freq!=0])) | all(is.na(var_exp_g_c_s$value[var_exp_g_c_s$Freq==0]))){
            full_cancer_gene_stat=NULL
          }else{
            # run GLM
            w = wilcox.test(var_exp_g_c_s$value[var_exp_g_c_s$Freq==0],var_exp_g_c_s$value[var_exp_g_c_s$Freq!=0])
            wP = w$p.value
            wWstat = w$statistic
            
            #df=cbind(var_exp_g_c_s,covariats[pmatch(var_exp_g_c_s$bcr_patient_barcode,covariats$bcr_patient_barcode),])
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
          
            cancer_gene_stat=run_glm(df,yi="value",xi="Freq",covi=c("Age","PC1","PC2"),ytype=pheno[3],gene=gene,cancer=cancer)#"TCGA_Subtype",
            if(pheno[3]=="Poisson" | pheno[3]=="Binary"){
              cancer_gene_stat=cancer_gene_stat[-which(colnames(cancer_gene_stat)%in%c("F","z.value","xi_lvl1"))]
            }else{
              cancer_gene_stat=cancer_gene_stat[-which(colnames(cancer_gene_stat)%in%c("F","t.value","xi_lvl1"))]
            }
            colnames(cancer_gene_stat) = c("yi","ytype","xi","df","Deviance","Resid. Df","Resid. Dev","Pr(>F)","Estimate","Std..Error","Pr...t..","covars");
   
            # compile results
            full_cancer_gene_stat = cbind(cancer,gene,clinicPheno,gene_path_count,wP,wWstat,cancer_gene_stat)
          }
            tt = rbind(tt, full_cancer_gene_stat)
        }
      }
  }
  
  print(grep(gene,genes))
  return(tt)
}

#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
colnames(finalMatrix) = c("cancer","gene","clinicPhenotype","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","p-value(Anova)","coefficient","StdError","p-value","covariants");


finalMatrix=finalMatrix[order(finalMatrix$`p-value`, decreasing=FALSE),]
tn = paste0("out/pathVarP_ImmuneAssoc_genelevel",label,"byCancerType_cov_20200605.txt")

#tn = paste0("out/pathVarP_ImmuneAssoc_genelevel",label,"byCancerType_cov_20200409_neo_normal.txt")

write.table(finalMatrix, quote=F, sep="\t", file = tn, row.names = F)


