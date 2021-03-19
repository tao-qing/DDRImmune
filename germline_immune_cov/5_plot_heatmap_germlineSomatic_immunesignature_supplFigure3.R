#Figure 4 heatmap
#heatmap for DDR and neoantigen interaction
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/")

#######################
####Part1 Load data####
#######################
source("../global_aes_out.R")

genes=c("BRCA1","BRCA2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","PALB2","ATM","ATR","POLE","POLQ","CHEK2")
phenos=c("PDCD1","CD274","xCell_CD8+_T-cells","xCell_CD4+_T-cells","xCell_B-cells","xCell_Th1_cells","xCell_Th2_cells","xCell_DC","Lymphocyte_Infiltration_Signature_Score","CYTScore")

somatic=read.table("./out/somaticDriver_ImmuneAssoc_geneleveloverlaped_n9738byCancerType_cov_20200409.txt",h=T,sep="\t")
somatic=somatic[somatic$gene%in%genes,]
somatic=somatic[somatic$clinicPhenotype%in%phenos,]
germline=read.table("./out/pathVarP_ImmuneAssoc_geneleveloverlaped_n9738byCancerType_cov_20200409.txt",h=T,sep="\t")
germline=germline[germline$gene%in%genes,]
germline=germline[germline$clinicPhenotype%in%phenos,]
somatic$clinicPhenotype=gsub("_"," ",gsub("Lymphocyte_Infiltration_Signature_Score","TIL",gsub("Th2_cells","T helper type 2 (Th2) cells",gsub("Th1_cells","T helper type 1 (Th1) cells",gsub("^DC$","Dendritic cell",gsub("xCell_","",somatic$clinicPhenotype))))))
germline$clinicPhenotype=gsub("_"," ",gsub("Lymphocyte_Infiltration_Signature_Score","TIL",gsub("Th2_cells","T helper type 2 (Th2) cells",gsub("Th1_cells","T helper type 1 (Th1) cells",gsub("^DC$","Dendritic cell",gsub("xCell_","",germline$clinicPhenotype))))))

somatic$gene=paste0("s",somatic$gene)
germline$gene=paste0("g",germline$gene)

somatic=somatic[somatic$gene_path_count>3,]
germline=germline[germline$gene_path_count>3,]

merge=rbind(somatic,germline)

##################################
####pancancer heatmap####
##################################
merge$FDR=sapply(merge$p.value,function(x)p.adjust(x,n=28,method="fdr"))
plotMat=merge[merge$p.value<0.05,]
plotMat$label=paste0(plotMat$cancer," | ",plotMat$clinicPhenotype)
# pre-plotting
plotMat$plotP = -log10(plotMat$p.value) 

#sort gene
geneindex=names(table(plotMat$gene))[order(table(plotMat$gene),decreasing = T)]
plotMat$gene=factor(plotMat$gene,levels=geneindex)
#sort label
labelindex=names(table(plotMat$label))[order(table(plotMat$label),decreasing = F)]
plotMat$label=factor(plotMat$label,levels=labelindex)


p = ggplot(data=plotMat)
p = p + geom_tile(data=plotMat,aes(y=label, x=gene, fill= coefficient), linetype="blank") + scale_fill_gradient2(name= "", low="blue",mid="white",high="red",  midpoint = 0,na.value=NA, limit=c(-1.5,2.5))
p = p + geom_text(data=plotMat,aes(y=label, x=gene, label = round(coefficient,digits = 2)), color="black", size=2)
p = p + geom_tile(data=plotMat[plotMat$FDR < 0.15,],aes(y=label, x=gene), color="grey",fill=NA, size=1)
p = p + geom_tile(data=plotMat[plotMat$FDR < 0.05,],aes(y=label, x=gene), color="black",fill=NA, size=1)
p = p  + theme_bw() + theme_nogrid() +
  theme(legend.position = "right",axis.text.x = element_text(colour="black", size=12, angle=90, hjust = 0.95,face="italic"), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.title=element_text(size=14,face="bold"),axis.ticks = element_blank())
p = p + labs(title="",y="Genes",x = "Immune Response Signatures")
p

fn = './out/supplmentaryFigure2_heatmap_mutation_immuneResponseSignature_byCancertype.pdf'
ggsave(fn,height=8,width=8)

