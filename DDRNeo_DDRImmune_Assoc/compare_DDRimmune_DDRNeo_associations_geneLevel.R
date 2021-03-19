setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/DDRNeo_DDRImmune_Assoc/")

library(ggplot2)
library(ggrepel)

th= theme_bw()+ theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=12,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=18,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

genes=c("BRCA1","BRCA2","PALB2","MLH1","MSH2","MSH3","MSH6","PMS1","PMS2","POLE","POLQ","ATM","ATR","CHEK2")
####compare HR with DDRimmune coefficients
label="rmhypermutator"
sddrimmne=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/somaticDriver_ImmuneAssoc_geneleveloverlaped_n9738byCancerType_cov_20200409_neo_normal.txt",h=T,stringsAsFactors = F)
sddrimmne=sddrimmne[sddrimmne$gene%in%genes,]
sddrimmne=sddrimmne[sddrimmne$gene_path_count>=4,]
sddrimmne$gene=paste0("s",sddrimmne$gene)
gddrimmne=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/pathVarP_ImmuneAssoc_geneleveloverlaped_n9738byCancerType_cov_20200409_neo_normal.txt",h=T,stringsAsFactors = F)
gddrimmne=gddrimmne[gddrimmne$gene%in%genes,]
gddrimmne=gddrimmne[gddrimmne$gene_path_count>=4,]
gddrimmne$gene=paste0("g",gddrimmne$gene)
ddrimmne=rbind(sddrimmne,gddrimmne)

#label="allsample"
#ddrimmne=read.table("/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/somaticDriver_ImmuneAssoc_pathwayleveloverlaped_n10080byCancerType_cov_20200603.txt",h=T,stringsAsFactors = F)
ddrimmne$group=paste0(ddrimmne$cancer,":",ddrimmne$gene)

phenoes=c("PDCD1","CD274","CYTScore","Indel_Neoantigens","SNV_Neoantigens","Lymphocyte_Infiltration_Signature_Score")

pheno1="SNV_Neoantigens"
pheno2="Lymphocyte_Infiltration_Signature_Score"

function1<-function(pheno1,pheno2,xllim=-3,xhlim=4,yllim=-2,yhlim=2,shortName1,shortName2){
    ddrimmne_sub1=ddrimmne[ddrimmne$clinicPhenotype==pheno1,]
    ddrimmne_sub2=ddrimmne[ddrimmne$clinicPhenotype==pheno2,]
    overlap=intersect(ddrimmne_sub1$group,ddrimmne_sub2$group)
    #eval(parse(text=paste0("mat=as.data.frame(cbind(",shortName1,"=ddrimmne_sub1$coefficient[pmatch(overlap,ddrimmne_sub1$group)],",shortName2,"=ddrimmne_sub2$coefficient[pmatch(overlap,ddrimmne_sub2$group)],label=overlap))")))
    mat=as.data.frame(cbind(coef1=ddrimmne_sub1$coefficient[pmatch(overlap,ddrimmne_sub1$group)],coef2=ddrimmne_sub2$coefficient[pmatch(overlap,ddrimmne_sub2$group)],label=overlap))
    mat$coef1=as.numeric(as.matrix(mat$coef1))
    mat$coef2=as.numeric(as.matrix(mat$coef2))
    
    cc=round(cor(mat$coef1,mat$coef2),digits=2)
    pvalue=round(cor.test(mat$coef1,mat$coef2)$p.value,digits = 5)
    
    p = ggplot(mat, aes(x=coef1, y=coef2))
    p = p + geom_point(alpha=0.5)#+xlim(xllim,xhlim)+ylim(yllim,yhlim)
    p = p + geom_text_repel(label=mat$label,size=3,segment.alpha =0.5,fontface="bold.italic",segment.colour ="grey",min.segment.length = 0)
    p = p + geom_smooth(data=mat,aes(x=coef1, y =coef2),colour="red",method = "glm",alpha=0.1,size=1.5)
    p = p +labs(title=paste0("\nr=",cc,";p=",pvalue), size = 16, colour = "red")+xlab(shortName1)+ylab(shortName2)
    p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
    p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
    p=p+th
    p
    fn = paste0("./out/TCGA_DDR_",pheno1,"vs",pheno2,".pdf")
    ggsave(fn,w = 5, h = 5, useDingbat=F,limitsize = FALSE)
}

function1("SNV_Neoantigens","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo","DDR:TIL")

function1("Indel_Neoantigens","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo","DDR:TIL")

function1("SNV_Neoantigens","CYTScore",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo","DDR:CYTScore")

function1("Indel_Neoantigens","CYTScore",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo","DDR:CYTScore")

function1("SNV_Neoantigens","PDCD1",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo","DDR:PD1")

function1("Indel_Neoantigens","PDCD1",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo","DDR:PD1")

function1("SNV_Neoantigens","CD274",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo","DDR:PD-L1")

function1("Indel_Neoantigens","CD274",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo","DDR:PD-L1")



function1("TMB","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB","DDR:TIL")

function1("TMB","CYTScore",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB","DDR:CYTScore")

function1("TMB","PDCD1",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo","DDR:TMB")

function1("TMB","CD274",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB","DDR:PD-L1")




#neoantigens and TIL
xllim=-1.8
xhlim=2.2
yllim=-1.5
yhlim=1.8
subddrimmne_indelneo=ddrimmne[ddrimmne$clinicPhenotype=="Indel_Neoantigens",]
subddrimmne_SNVneo=ddrimmne[ddrimmne$clinicPhenotype=="SNV_Neoantigens",]
subddrimmne_TIL=ddrimmne[ddrimmne$clinicPhenotype=="Lymphocyte_Infiltration_Signature_Score",]

overlap=intersect(subddrimmne_indelneo$group,subddrimmne_TIL$group)
mat=as.data.frame(cbind(TCGA_Indel_coef=subddrimmne_indelneo$coefficient[pmatch(overlap,subddrimmne_indelneo$group)],TCGA_TIL_coef=subddrimmne_TIL$coefficient[pmatch(overlap,subddrimmne_TIL$group)],label=overlap))
mat$TCGA_Indel_coef=as.numeric(as.matrix(mat$TCGA_Indel_coef))
mat$TCGA_TIL_coef=as.numeric(as.matrix(mat$TCGA_TIL_coef))

cc=round(cor(mat$TCGA_Indel_coef,mat$TCGA_TIL_coef),digits=2)
pvalue=round(cor.test(mat$TCGA_Indel_coef,mat$TCGA_TIL_coef)$p.value,digits = 2)

p = ggplot(mat, aes(y=TCGA_Indel_coef, x=TCGA_TIL_coef))
p = p + geom_point(alpha=0.5)+xlim(xllim,xhlim)+ylim(yllim,yhlim)
p = p + geom_text_repel(label=mat$label,size=3,segment.alpha =0.5,fontface="bold.italic",segment.colour ="grey",min.segment.length = 0)
p = p + geom_smooth(data=mat,aes(x=TCGA_TIL_coef, y =TCGA_Indel_coef),colour="red",method = "glm",alpha=0.1,size=1.5)
p = p +labs(title=paste0("\nr=",cc,";p=",pvalue), size = 16, colour = "red")
p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
p=p+th
p
fn = paste0("./out/Indel_Neo_TIL_TCGAvsMSKCC.pdf")
ggsave(fn,w = 5, h = 5, useDingbat=F,limitsize = FALSE)





subddrimmne_indelneo=ddrimmne[ddrimmne$clinicPhenotype=="Indel_Neoantigens",]
subddrimmne_SNVneo=ddrimmne[ddrimmne$clinicPhenotype=="SNV_Neoantigens",]
subddrimmne_TIL=ddrimmne[ddrimmne$clinicPhenotype=="Lymphocyte_Infiltration_Signature_Score",]

overlap=intersect(subddrimmne_TIL$group,subddrimmne_SNVneo$group)
mat=as.data.frame(cbind(TCGA_SNV_coef=subddrimmne_SNVneo$coefficient[pmatch(overlap,subddrimmne_SNVneo$group)],TCGA_TIL_coef=subddrimmne_TIL$coefficient[pmatch(overlap,subddrimmne_TIL$group)],label=overlap))
mat$TCGA_SNV_coef=as.numeric(as.matrix(mat$TCGA_SNV_coef))
mat$TCGA_TIL_coef=as.numeric(as.matrix(mat$TCGA_TIL_coef))

cc=round(cor(mat$TCGA_SNV_coef,mat$TCGA_TIL_coef),digits=2)
pvalue=round(cor.test(mat$TCGA_SNV_coef,mat$TCGA_TIL_coef)$p.value,digits = 5)

p = ggplot(mat, aes(y=TCGA_SNV_coef, x=TCGA_TIL_coef))
p = p + geom_point(alpha=0.5)+xlim(xllim,xhlim)+ylim(yllim,yhlim)
p = p + geom_text_repel(label=mat$label,size=3,segment.alpha =0.5,fontface="bold.italic",segment.colour ="grey",min.segment.length = 0)
p = p + geom_smooth(data=mat,aes(x=TCGA_TIL_coef, y =TCGA_SNV_coef),colour="red",method = "glm",alpha=0.1,size=1.5)
p = p +labs(title=paste0("\nr=",cc,";p=",pvalue), size = 16, colour = "red")
p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
p=p+th
p
fn = paste0("./out/SNV_Neo_TIL_TCGAvsMSKCC.pdf")
ggsave(fn,w = 5, h = 5, useDingbat=F,limitsize = FALSE)
