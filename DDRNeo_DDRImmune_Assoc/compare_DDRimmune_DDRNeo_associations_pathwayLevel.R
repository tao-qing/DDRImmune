setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/DDRNeo_DDRImmune_Assoc/")

cancer = unique(clin$type)
genes = list()
genes=NULL
genes["HR"]=list(c("BRCA1","BRCA2","PALB2"))
genes["Polymerase"]=list(c("POLE","POLQ"))
genes["Sensor"]=list(c("ATM","ATR","CHEK2"))
genes["MMR"]=list(c("MLH1","MSH2","MSH3","MSH6","PMS2"))

library(ggplot2)
library(ggrepel)

th= theme_bw()+ theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=20,hjust=0.90), axis.text.y = element_text(colour="black", size=20,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=18,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5,size=18),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))

#cal.pvalue<-function1<-function(pheno1,pheno2,xllim=-3,xhlim=4,yllim=-2,yhlim=2,shortName1,shortName2,cate,type){
#    ddrimmne_sub1=ddrimmne[ddrimmne$clinicPhenotype==pheno1,]
#    ddrimmne_sub2=ddrimmne[ddrimmne$clinicPhenotype==pheno2,]
#    overlap=intersect(ddrimmne_sub1$group,ddrimmne_sub2$group)
    
    
#    mat=as.data.frame(cbind(coef1=ddrimmne_sub1$coefficient[pmatch(overlap,ddrimmne_sub1$group)],coef2=ddrimmne_sub2$coefficient[pmatch(overlap,ddrimmne_sub2$group)],label=overlap,gene=sapply(overlap,function(x)strsplit(x,split="\\:")[[1]][2])))
#    mat$coef1=as.numeric(as.matrix(mat$coef1))
#    mat$coef2=as.numeric(as.matrix(mat$coef2))
    
    #showlabel=intersect(mat$label,ddrimmne_sub2$group[p.adjust(ddrimmne_sub2$p.value,method="fdr")<0.05])
    
#    cc=round(cor(mat$coef1,mat$coef2),digits=2)
#    pvalue=p.adjust(cor.test(mat$coef1,mat$coef2)$p.value,method="fdr",n=4)
#    pvalue=ifelse(pvalue<0.001,formatC(pvalue, format = "e", digits = 2),formatC(signif(pvalue,digits=2), digits=2,format="fg", flag="#"))
#}

function1("SNV_Neoantigens","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo (coef)","DDR:TILs (coef)",cate=cate,type=type)
pheno1="SNV_Neoantigens";pheno2="Lymphocyte_Infiltration_Signature_Score";xllim=-2;xhlim=2.1;yllim=-0.75;yhlim=0.5;shortName1="DDR:SNV_neo (coef)";shortName2="DDR:TILs (coef)";cate=cate;type=type

function1<-function(pheno1,pheno2,xllim=-3,xhlim=4,yllim=-2,yhlim=2,shortName1,shortName2,cate,type){
    ddrimmne_sub1=ddrimmne[ddrimmne$clinicPhenotype==pheno1,]
    ddrimmne_sub2=ddrimmne[ddrimmne$clinicPhenotype==pheno2,]
    overlap=intersect(ddrimmne_sub1$group,ddrimmne_sub2$group)
    

    mat=as.data.frame(cbind(coef1=ddrimmne_sub1$coefficient[pmatch(overlap,ddrimmne_sub1$group)],coef2=ddrimmne_sub2$coefficient[pmatch(overlap,ddrimmne_sub2$group)],label=overlap,gene=sapply(overlap,function(x)strsplit(x,split="\\:")[[1]][2])))
    mat$coef1=as.numeric(as.matrix(mat$coef1))
    mat$coef2=as.numeric(as.matrix(mat$coef2))
    if(pheno1=="Indel_Neoantigens"){
        mat=mat[mat$label!="SARC:sSensor",]
    }
    showlabel=intersect(mat$label,ddrimmne_sub2$group[ddrimmne_sub2$FDR<0.05])
    
    cc=round(cor(mat$coef1,mat$coef2),digits=2)
    #pvalue=p.adjust(cor.test(mat$coef1,mat$coef2)$p.value,method="fdr",n=4)
    pvalue=cor.test(mat$coef1,mat$coef2)$p.value
    pvalue=ifelse(pvalue<0.001,formatC(pvalue, format = "e", digits = 1),formatC(signif(pvalue,digits=1), digits=1,format="fg", flag="#"))
   
    
    p = ggplot(mat, aes(x=coef1, y=coef2,colour=gene))
    p = p + geom_smooth(data=mat,aes(x=coef1, y =coef2),colour="red",method = "glm",alpha=0.1,size=1.5)
    p = p + geom_point(alpha=0.5)#+xlim(xllim,xhlim)+ylim(yllim,yhlim)
    p = p + geom_text_repel(label=ifelse(mat$label%in%showlabel,as.character(mat$label),""),size=5,segment.alpha =0.5,fontface="bold.italic",segment.colour ="black",min.segment.length = 0)
    #if(pvalue<1e-5){
    #    p = p +labs(title=paste0("\nr=",cc,";p<1e-5"), size = 16, colour = "red")+xlab(shortName1)+ylab(shortName2)
   # }else{
        p = p +labs(title=paste0("\nr=",cc,";P=",pvalue), size = 20, colour = "red")+xlab(shortName1)+ylab(shortName2)
   # }
    
    p = p + geom_hline(yintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
    p = p + geom_vline(xintercept=0, linetype="dashed", color = "grey",alpha = I(0.6))
    p=p+th
    p
    fn = paste0("./out/TCGA_DDR_",pheno1,"vs",pheno2,"_",cate,"_",type,"_20201005.pdf")
    ggsave(fn,w = 4, h = 4, useDingbat=F,limitsize = FALSE)
}

pheno1="Indel_Neoantigens"
pheno2="Lymphocyte_Infiltration_Signature_Score"

####compare HR with DDRimmune coefficients
cate="rmhypermutator"
type="somatic"
if(cate=="rmhypermutator"){
    if(type=="somatic"){
       fn="/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/somaticDriver_ImmuneAssoc_pathwayleveloverlaped_n9738byCancerType_cov_20200603.txt"
       label="s"
    }else{
        fn="/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/germline_ImmuneAssoc_pathwayleveloverlaped_n9738byCancerType_cov_20200603.txt"
        label="g"
    }
}else{
    if(type=="somatic"){
        fn="/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/somaticDriver_ImmuneAssoc_pathwayleveloverlaped_n10080byCancerType_cov_20200603.txt"
        label="s"
    }else{
        fn="/Users/qingtao/Box Sync/GermlineSomatic/analysis/germline_immune_cov/out/germline_ImmuneAssoc_pathwayleveloverlaped_n10080byCancerType_cov_20200603.txt"
        label="g"
    }
}


ddrimmne1=read.table(fn,h=T,stringsAsFactors = F)%>%filter(gene%in%names(genes) & clinicPhenotype%in%c("PDCD1","Lymphocyte_Infiltration_Signature_Score","CYTScore","CD274") & gene_path_count>=4)%>%mutate(gene=paste0(label,gene),group=paste0(cancer,":",gene),FDR=p.adjust(p.value,method = "fdr"))

ddrimmne2=read.table(fn,h=T,stringsAsFactors = F)%>%filter(gene%in%names(genes) & clinicPhenotype%in%c("SNV_Neoantigens","Indel_Neoantigens") & gene_path_count>=4)%>%mutate(gene=paste0(label,gene),group=paste0(cancer,":",gene),FDR=1)

ddrimmne=rbind(ddrimmne1,ddrimmne2)

comparisons=list(c("SNV_Neoantigens","Lymphocyte_Infiltration_Signature_Score"),c("Indel_Neoantigens","Lymphocyte_Infiltration_Signature_Score"),c("SNV_Neoantigens","CYTScore"),c("Indel_Neoantigens","CYTScore"),c("SNV_Neoantigens","PDCD1"),c("Indel_Neoantigens","PDCD1"),c("SNV_Neoantigens","CD274"),c("Indel_Neoantigens","CD274"))



function1("SNV_Neoantigens","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo (coef)","DDR:TILs (coef)",cate=cate,type=type)

function1("Indel_Neoantigens","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo (coef)","DDR:TILs (coef)",cate=cate,type=type)

function1("SNV_Neoantigens","CYTScore",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo (coef)","DDR:CYTScore (coef)",cate=cate,type=type)

function1("Indel_Neoantigens","CYTScore",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo (coef)","DDR:CYTScore (coef)",cate=cate,type=type)

function1("SNV_Neoantigens","PDCD1",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:SNV_neo (coef)","DDR:PD1 (coef)",cate=cate,type=type)

function1("Indel_Neoantigens","PDCD1",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:Indel_neo (coef)","DDR:PD1 (coef)",cate=cate,type=type)

function1("SNV_Neoantigens","CD274",xllim=-2,xhlim=4.1,yllim=-0.75,yhlim=2.5,"DDR:SNV_neo (coef)","DDR:PD-L1 (coef)",cate=cate,type=type)

function1("Indel_Neoantigens","CD274",xllim=-5,xhlim=5,yllim=-5,yhlim=5,"DDR:Indel_neo (coef)","DDR:PD-L1 (coef)",cate=cate,type=type)


function1("TMB","Lymphocyte_Infiltration_Signature_Score",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB (coef)","DDR:TIL (coef)",cate=cate,type=type)

function1("TMB","CYTScore",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB (coef)","DDR:CYTScore (coef)",cate=cate,type=type)

function1("TMB","PDCD1",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB (coef)","DDR:PD1 (coef)",cate=cate,type=type)

function1("TMB","CD274",xllim=-2,xhlim=2.1,yllim=-0.75,yhlim=0.5,"DDR:TMB (coef)","DDR:PD-L1 (coef)",cate=cate,type=type)








#figure legends
library(gridExtra)
library(grid)

pheno1="Indel_Neoantigens"
pheno2="CYTScore"
ddrimmne_sub1=ddrimmne[ddrimmne$clinicPhenotype==pheno1,]
ddrimmne_sub2=ddrimmne[ddrimmne$clinicPhenotype==pheno2,]
overlap=intersect(ddrimmne_sub1$group,ddrimmne_sub2$group)
#eval(parse(text=paste0("mat=as.data.frame(cbind(",shortName1,"=ddrimmne_sub1$coefficient[pmatch(overlap,ddrimmne_sub1$group)],",shortName2,"=ddrimmne_sub2$coefficient[pmatch(overlap,ddrimmne_sub2$group)],label=overlap))")))
mat=as.data.frame(cbind(coef1=ddrimmne_sub1$coefficient[pmatch(overlap,ddrimmne_sub1$group)],coef2=ddrimmne_sub2$coefficient[pmatch(overlap,ddrimmne_sub2$group)],label=overlap,gene=sapply(overlap,function(x)strsplit(x,split="\\:")[[1]][2])))
mat$coef1=as.numeric(as.matrix(mat$coef1))
mat$coef2=as.numeric(as.matrix(mat$coef2))

showlabel=intersect(mat$label,ddrimmne_sub2$group[p.adjust(ddrimmne_sub2$p.value,method="fdr")<0.15])

p = ggplot(mat, aes(x=coef1, y=coef2,colour=gene))
p = p + geom_point(alpha=0.5)#+xlim(xllim,xhlim)+ylim(yllim,yhlim)
p=p+theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(colour="black", size=12,hjust=0.90), axis.text.y = element_text(colour="black", size=12,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"),axis.title=element_text(size=18,face="bold"))+ theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=12),strip.text.y = element_text(angle = 0,size=10,face="italic"),panel.spacing = unit(0.1, "lines"),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))


g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

mylegend<-g_legend(p)
legend <- grid.draw(mylegend)

grid.newpage()

pdf("./out/figure_legend.pdf",width=2,height = 2)
p#grid.draw(legend) 
dev.off()



