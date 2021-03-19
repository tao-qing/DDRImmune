#Figure 4 heatmap
#heatmap for DDR and neoantigen interaction
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_contribution/")

#######################
####Part1 Load data####
#######################
source("../global_aes_out.R")
source("../dependency_files_tq.R")


########################################
####annotate with protein alteration####
####time consuming, precalculated#######
########################################
anno=read.table("./out/somatic_mutation_annotation.txt",sep="\t",h=T)
##############Indel annotation with amino acid change######################
tnpan = "./out/pancancer_mutation_Indel_neoantigen_association.txt"
ttpan = read.table(sep="\t",header=T,file=tnpan, stringsAsFactors=FALSE)
neonumb=ttpan$mutationWithNeo+ttpan$wildwithNeo
mutnumb=ttpan$mutationWithNeo+ttpan$mutationNon

ttpan=ttpan[neonumb>=4 & mutnumb>=4 & ttpan$mutationWithNeo>1,]

tmp=sapply(ttpan$Indel,function(x){
  chr=strsplit(x,split="_")[[1]][1] 
  start=strsplit(x,split="_")[[1]][2] 
  end=strsplit(x,split="_")[[1]][3] 
  ref=strsplit(x,split="_")[[1]][4] 
  alt=strsplit(x,split="_")[[1]][5]
  if(nchar(ref)==1){
    newref="-"
    newalt=substr(alt,2,nchar(alt))
    mut=paste0(c(chr,start,end,newref,newalt),collapse = "_")
  }
  
  if(nchar(alt)==1){
    newstart=as.numeric(as.matrix(start))+1
    newref=substr(ref,2,nchar(ref))
    newalt="-"
    mut=paste0(c(chr,newstart,end,newref,newalt),collapse = "_")
  }
  return(mut)
})
ttpan=cbind(ttpan,anno[unlist(sapply(tmp,function(x)grep(x,anno$mutation))),])
write.table(ttpan,"./out/mutation_Indel_neoantigen_association_annotated_new.txt",row.names = F,quote = F,sep = "\t")


tn = "./out/mutation_Indel_neoantigen_association_byCancer.txt"
tt = read.table(sep="\t",header=T,file=tn, stringsAsFactors=FALSE)
neonumb=tt$mutationWithNeo+tt$wildwithNeo
mutnumb=tt$mutationWithNeo+tt$mutationNon
tmp=sapply(tt$Indel,function(x){
  chr=strsplit(x,split="_")[[1]][1] 
  start=strsplit(x,split="_")[[1]][2] 
  end=strsplit(x,split="_")[[1]][3] 
  ref=strsplit(x,split="_")[[1]][4] 
  alt=strsplit(x,split="_")[[1]][5]
  if(nchar(ref)==1){
    newref="-"
    newalt=substr(alt,2,nchar(alt))
    mut=paste0(c(chr,start,end,newref,newalt),collapse = "_")
  }
  
  if(nchar(alt)==1){
    newstart=as.numeric(as.matrix(start))+1
    newref=substr(ref,2,nchar(ref))
    newalt="-"
    mut=paste0(c(chr,newstart,end,newref,newalt),collapse = "_")
  }
  return(mut)
})
tt=cbind(tt,anno[sapply(tmp,function(x)grep(x,anno$mutation)),])
write.table(tt,"./out/mutation_Indel_neoantigen_association_byCancer_annotated.txt",row.names = F,quote = F,sep = "\t")

##################################
####Figure4C pancancer heatmap####
##################################
#pan-cancer level Indel
ttpan = read.table(sep="\t",header=T,file="./out/mutation_Indel_neoantigen_association_annotated_new.txt", stringsAsFactors=FALSE)
ttpan$FDR=p.adjust(ttpan$P,method="fdr")
ttpan1=ttpan[ttpan$P<0.05,]
ttpanG=ttpan1
ttpanG$label=paste0(ttpanG$mutGene," | ",ttpanG$aa)
# pre-plotting
ttpanG$minusLogP = -log10(ttpanG$P) 
ttpanG$plotP = ttpanG$minusLogP
ttpanG$plotP[ttpanG$OR < 1] = -ttpanG$plotP[ttpanG$OR < 1] # opposite effect size (mutual exclusivity)
ttpanG$plotP[ttpanG$plotP > 5] = 5
ttpanG$plotP[ttpanG$plotP < -5] = -5

### plotting ###
ttpanG=ttpanG[ttpanG$wildwithNeo>1,]

###select variant####
#variants with at least one FDR<0.15
neosig=unique(ttpanG$label[ttpanG$FDR<0.15])
ttpanG=ttpanG[ttpanG$label%in%neosig,]
###variant order###
neonumb=ttpanG$mutationWithNeo+ttpanG$wildwithNeo
mutnumb=ttpanG$mutationWithNeo+ttpanG$mutationNon
freq=neonumb/(neonumb+mutnumb)
index=unique(ttpanG$label[order(freq,decreasing = T)])
ttpanG$label=factor(ttpanG$label,levels=index)

ttpanG$gene=gsub("sDS","sSensor",gsub("sDP","sPolymerase",ttpanG$gene))

ttpanG$freq=round(((ttpanG$mutationWithNeo)/(ttpanG$mutationWithNeo+ttpanG$mutationNon))*100,digits=1)

write.table(ttpanG,"./out/DDRsignificant_enriched_Indel_neoantigens_rmhypermutator_pancancerLevel.txt",sep="\t")


ttpanG=read.table("./out/DDRsignificant_enriched_Indel_neoantigens_rmhypermutator_pancancerLevel.txt",h=T,stringsAsFactors = F)
hs=c("DOCK3 | p.P1852Qfs*45","DAZAP1 | p.P257Rfs*78","RNF43 | p.G659Vfs*41")
ttpanG1=ttpanG[ttpanG$label%in%hs,]

p = ggplot(data=ttpanG)
p = p + geom_tile(data=ttpanG,aes(y=gene, x=label, fill= as.numeric(as.matrix(freq))), linetype="blank") + scale_fill_gradientn(name= "Frequency (%)", colours=c("white","orange","red"), na.value=NA, limit=c(0,6))#+coord_flip()
#p = p + geom_text(data=ttpanG,aes(y=gene, x=label, label =as.numeric(as.matrix(freq))), color="black", size=3)
p = p + geom_tile(data=ttpanG[ttpanG$FDR < 0.15,],aes(y=gene, x=label), color="grey",fill=NA, size=1) #+ scale_color_gradien
p = p + geom_tile(data=ttpanG[ttpanG$FDR < 0.05,],aes(y=gene, x=label), color="black",fill=NA, size=1) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
p = p  + theme_bw() + theme_nogrid() +
  theme(legend.position = "top",legend.title=element_text(size=16),legend.text=element_text(size=14),axis.text.x = element_text(colour="black", size=16, angle=60, hjust = 1,face="italic"), axis.text.y = element_text(colour="black", size=16,hjust = 0.95),axis.title=element_text(size=16,face="bold"),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + labs(title="",y="DDR genes",x = "Indel Neoantigens")
p
fn = './out/mutation_Indelneoantigen_association_heatmap.pdf'
ggsave(fn,height=6,width=14.5)

##################################
####Figure3C heatmap by cancer####
##################################
tt = read.table(sep="\t",header=T,file="./out/mutation_Indel_neoantigen_association_byCancer_annotated.txt", stringsAsFactors=FALSE)
tt=tt[tt$P!=1,]
tt$FDR=p.adjust(tt$P)
tt1=tt[tt$P<0.05,]
ttG=tt1
ttG$label=paste0(ttG$cancer," | ",ttG$gene.1," | ",ttG$aa)
# pre-plotting
ttG$minusLogP = -log10(ttG$P) 
ttG$plotP = ttG$minusLogP
ttG$plotP[ttG$OR < 1] = -ttG$plotP[ttG$OR < 1] # opposite effect size (mutual exclusivity)
ttG$plotP[ttG$plotP > 5] = 5
ttG$plotP[ttG$plotP < -5] = -5

### plotting ###
ttG=ttG[ttG$wildwithNeo>1,]
ttG$gene=gsub("BRCAgene","BRCA1/2",ttG$gene)
ttG$gene=gsub("ATM_ATR","ATM/ATR",ttG$gene)
ttG$gene=gsub("POLE_POLQ","POLE/POLQ",ttG$gene)

p = ggplot(data=ttG)
p = p + geom_tile(data=ttG,aes(y=gene, x=label, fill= OR), linetype="blank") + scale_fill_gradientn(name= "Odds Ratio", colours=c("white","orange","red"), na.value=NA, limit=c(3,40))#+coord_flip()
p = p + geom_text(data=ttG,aes(y=gene, x=label, label = mutationWithNeo), color="black", size=3)
p = p + geom_tile(data=ttG[ttG$P < 0.05,],aes(y=gene, x=label), color="grey",fill=NA, size=1) #+ scale_color_gradien
p = p + geom_tile(data=ttG[ttG$FDR < 0.05,],aes(y=gene, x=label), color="black",fill=NA, size=1) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
p = p  + theme_bw() + theme_nogrid() +
  theme(legend.position = "top",axis.text.x = element_text(colour="black", size=10, angle=60, hjust = 1,face="italic"), axis.text.y = element_text(colour="black", size=14,face="italic",hjust = 0.95),axis.title=element_text(size=14,face="bold"),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + labs(title="",y="DDR genes",x = "Indel Neoantigens")
p
fn = './out/mutation_Indelneoantigen_association_heatmap_byCancertype.pdf'
ggsave(fn,height=4.5,width=5)


#significant neoantigens
neo=ttpanG$Indel
write.table(neo,"./out/DDRsignificant_enriched_Indel_neoantigens_rmhypermutator.txt")









