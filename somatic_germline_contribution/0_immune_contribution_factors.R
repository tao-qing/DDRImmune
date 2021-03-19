setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")

#####################gene level contributors######################
sigGeneAssco=function(assoc,pheno,iterm,FDR,mutCounts){
  subAssoc=assoc[which(assoc$clinicPhenotype%in%iterm),]
  subAssoc=subAssoc[subAssoc$gene_path_count>=mutCounts,]
  #FDR
  #samplesize=table(subAssoc$cancer)
  #subAssoc$FDR = sapply(1:dim(subAssoc)[1],function(x){p.adjust(subAssoc$p.value[x],method="fdr",n=samplesize[subAssoc$cancer[x]])})
  subAssoc$FDR=p.adjust(subAssoc$p.value,method="fdr")
  subAssoc=subAssoc[subAssoc$FDR<FDR,]
  if(length(which(is.na(subAssoc$cancer)))){
    subAssoc=subAssoc[-which(is.na(subAssoc$cancer)),]
  }
  return(subAssoc)
}
pheno="Immune"

sPathImmune=read.table("../germline_immune_cov/out/somaticDriver_ImmuneAssoc_genelevel_rmhypermutator_byCancerType_cov_20200204.txt",h=T,sep="\t",stringsAsFactors = F)
gPathImmune=read.table("../germline_immune_cov/out/pathVarP_ImmuneAssoc_genelevel_rmhypermutator_byCancerType_cov_20200204.txt",h=T,sep="\t",stringsAsFactors = F)


  iterm=as.character(read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/allAssoc_iterms",sep="\t")[,2])
  
  gPathImmuneSig=sigGeneAssco(gPathImmune,pheno="Immune",iterm,FDR=0.05,mutCounts=4)
  sPathImmuneSig=sigGeneAssco(sPathImmune,pheno="Immune",iterm,FDR=0.05,mutCounts=4)
  
  immunePhenotype=unique(c(gPathImmuneSig$clinicPhenotype,sPathImmuneSig$clinicPhenotype))
  cancers=unique(c(gPathImmuneSig$cancer,sPathImmuneSig$cancer))
  
  cMat=NULL
  for(c in cancers){
    for(i in immunePhenotype){
      germlineFactor=paste(unique(gPathImmuneSig$gene[intersect(which(gPathImmuneSig$cancer==c),which(gPathImmuneSig$clinicPhenotype==i))]),collapse = " + ")
      somaticFactor=paste(unique(sPathImmuneSig$gene[intersect(which(sPathImmuneSig$cancer==c),which(sPathImmuneSig$clinicPhenotype==i))]),collapse = " + ")
      if(length(germlineFactor)==0 || germlineFactor==""){
        germlineFactor=NA
      }
      if(length(somaticFactor)==0 || somaticFactor==""){
        somaticFactor=NA
      }
      contributorMat=cbind(c,i,germlineFactor,somaticFactor)
      
      cMat=rbind(cMat,contributorMat)
    }
  }
  
  cMat=cMat[-intersect(which(is.na(cMat[,3])),which(is.na(cMat[,4]))),]
  
  colnames(cMat)=c("Cancer","ImmunePhenotype","Germline","Somatic")
  cn =paste0("./out/",pheno,"_geneLevel_germline_somatoc_contributors_20200204.tsv")
  write.table(cMat, quote=F, sep="\t", file = cn, row.names = F)

  
  
  
    


#####################pathway level contributors######################
#load association between immune and pathway level alteration 

#select true positive assosication, (1) FDR < 0.01 within cancer type, (2) # of mutated samples >=4
sigPathAssco=function(assoc,pheno,iterm,path,FDR,mutCounts){
  subAssoc=assoc[which(assoc$clinicPhenotype%in%iterm),]
  subAssoc=subAssoc[subAssoc$gene_path_count>=mutCounts,]
  
  subAssoc=subAssoc[subAssoc$pathway%in%path,]
  #FDR
  samplesize=table(subAssoc$cancer)
  subAssoc$FDR = sapply(1:dim(subAssoc)[1],function(x){p.adjust(subAssoc$p.value[x],method="fdr",n=samplesize[subAssoc$cancer[x]])})
  subAssoc=subAssoc[subAssoc$FDR<FDR,]
  return(subAssoc)
}


sPathImmune=read.table("../germline_immune_cov/out/somaticDriver_ImmuneAssoc_pathwaylevel_rmhypermutator_byCancerType_cov_20190708.txt",h=T,sep="\t",stringsAsFactors = F)
gPathImmune=read.table("../germline_immune_cov/out/pathVarP_ImmuneAssoc_pathwaylevel_rmhypermutator_byCancerType_cov_20190708.txt",h=T,sep="\t",stringsAsFactors = F)

for(pheno in c("Immune")){#"Cybersort","xCell"
    pathway_info=read.table(sep="\t",header=F,file="../germline_immune_cov/out/pathway_info.tsv", stringsAsFactors=FALSE)
    path=pathway_info[which(pathway_info$V2%in%c("DDR","NanoString")),1]
    
    fn=paste0("../germline_immune_cov/out/",pheno,"_Iterms")
    iterm= as.character(read.table(sep="\t",header=F,file=fn, stringsAsFactors=FALSE)[,1])
 
    gPathImmuneSig=sigPathAssco(gPathImmune,pheno="Immune",iterm,path,FDR=0.05,mutCounts=4)
    sPathImmuneSig=sigPathAssco(sPathImmune,pheno="Immune",iterm,path,FDR=0.05,mutCounts=4)
    
    immunePhenotype=unique(c(gPathImmuneSig$clinicPhenotype,sPathImmuneSig$clinicPhenotype))
    cancers=unique(c(gPathImmuneSig$cancer,sPathImmuneSig$cancer))
    pathways=unique(c(gPathImmuneSig$pathway,sPathImmuneSig$pathway))
    
    cMat=NULL
    for(c in cancers){
        for(i in immunePhenotype){
            germlineFactor=paste(unique(gPathImmuneSig$pathway[intersect(which(gPathImmuneSig$cancer==c),which(gPathImmuneSig$clinicPhenotype==i))]),collapse = " + ")
            somaticFactor=paste(unique(sPathImmuneSig$pathway[intersect(which(sPathImmuneSig$cancer==c),which(sPathImmuneSig$clinicPhenotype==i))]),collapse = " + ")
            if(length(germlineFactor)==0 || germlineFactor==""){
              germlineFactor=NA
            }
            if(length(somaticFactor)==0 || somaticFactor==""){
              somaticFactor=NA
            }
            contributorMat=cbind(c,i,germlineFactor,somaticFactor)
            
            cMat=rbind(cMat,contributorMat)
      }
    }
    
    cMat=cMat[-intersect(which(is.na(cMat[,3])),which(is.na(cMat[,4]))),]
    
    colnames(cMat)=c("Cancer","ImmunePhenotype","Germline","Somatic")
    cn =paste0("./out/",pheno,"_DDR_NanoString_pathwalLevel_germline_somatoc_contributors.tsv")
    write.table(cMat, quote=F, sep="\t", file = cn, row.names = F)
}  



#####################BRCA1/2 and total germline variants contributors######################

#select true positive assosication, (1) FDR < 0.01 within cancer type, (2) # of mutated samples >=4
sigPathAssco=function(assoc,pheno,iterm,path,FDR,mutCounts){
  subAssoc=assoc[which(assoc$clinicPhenotype%in%iterm),]
  subAssoc=subAssoc[subAssoc$gene_path_count>=mutCounts,]
  
  subAssoc=subAssoc[subAssoc$pathway%in%path,]
  #FDR
  samplesize=table(subAssoc$cancer)
  subAssoc$FDR = sapply(1:dim(subAssoc)[1],function(x){p.adjust(subAssoc$p.value[x],method="fdr",n=samplesize[subAssoc$cancer[x]])})
  subAssoc=subAssoc[subAssoc$FDR<FDR,]
  return(subAssoc)
}


sPathImmune=read.table("../germline_immune_cov/out/somaticDriver_ImmuneAssoc_pathwaylevel_rmhypermutator_byCancerType_cov_20190708.txt",h=T,sep="\t",stringsAsFactors = F)
gPathImmune=read.table("../germline_immune_cov/out/pathVarP_ImmuneAssoc_pathwaylevel_rmhypermutator_byCancerType_cov_20190708.txt",h=T,sep="\t",stringsAsFactors = F)

for(pheno in c("Immune")){#"Cybersort","xCell"
  pathway_info=read.table(sep="\t",header=F,file="../germline_immune_cov/out/pathway_info.tsv", stringsAsFactors=FALSE)
  path=pathway_info[which(pathway_info$V2%in%c("DDR","NanoString")),1]
  
  fn=paste0("../germline_immune_cov/out/",pheno,"_Iterms")
  iterm= as.character(read.table(sep="\t",header=F,file=fn, stringsAsFactors=FALSE)[,1])
  
  gPathImmuneSig=sigPathAssco(gPathImmune,pheno="Immune",iterm,path,FDR=0.05,mutCounts=4)
  sPathImmuneSig=sigPathAssco(sPathImmune,pheno="Immune",iterm,path,FDR=0.05,mutCounts=4)
  
  immunePhenotype=unique(c(gPathImmuneSig$clinicPhenotype,sPathImmuneSig$clinicPhenotype))
  cancers=unique(c(gPathImmuneSig$cancer,sPathImmuneSig$cancer))
  pathways=unique(c(gPathImmuneSig$pathway,sPathImmuneSig$pathway))
  
  cMat=NULL
  for(c in cancers){
    for(i in immunePhenotype){
      germlineFactor=paste(unique(gPathImmuneSig$pathway[intersect(which(gPathImmuneSig$cancer==c),which(gPathImmuneSig$clinicPhenotype==i))]),collapse = " + ")
      somaticFactor=paste(unique(sPathImmuneSig$pathway[intersect(which(sPathImmuneSig$cancer==c),which(sPathImmuneSig$clinicPhenotype==i))]),collapse = " + ")
      if(length(germlineFactor)==0 || germlineFactor==""){
        germlineFactor=NA
      }
      if(length(somaticFactor)==0 || somaticFactor==""){
        somaticFactor=NA
      }
      contributorMat=cbind(c,i,germlineFactor,somaticFactor)
      
      cMat=rbind(cMat,contributorMat)
    }
  }
  
  cMat=cMat[-intersect(which(is.na(cMat[,3])),which(is.na(cMat[,4]))),]
  
  colnames(cMat)=c("Cancer","ImmunePhenotype","Germline","Somatic")
  cn =paste0("./out/",pheno,"_DDR_NanoString_pathwalLevel_germline_somatoc_contributors.tsv")
  write.table(cMat, quote=F, sep="\t", file = cn, row.names = F)
}  






