##### somatic_germline_overlap.R #####

# Find overlap of genes/variants for somatic/germline variants
setwd("/Users/qingtao/Box Sync/GermlineSomatic/analysis/somatic_germline_overlap/")

#epig<-data.frame(readxl::read_xlsx("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/TCGA_DDR_Data_Resources.xlsx",sheet = "DDR epigenetic silencing",col_names =F))
#epig<-epig[,-1]
#epig=t(epig)
#colnames(epig)=epig[1,]
#epig=epig[-1,]

overlapsplrmhyper<-read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/all_overlaped_samples_removehypermutator_n9738.txt")[,1]

#all core DDR genes
ddr=read.csv("../../Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/DDR_Pathways.csv",h=T)
ddrgene=unique(as.character(unlist(ddr)))[-1]

#ddrgene=c("MLH1","MSH2","MSH3","MSH6","PMS1","PMS2")

#ddrepig=as.data.frame(epig[,c("TCGA Sample (tumor type abbr. below)","Gene Symbol",ddrgene)])
#colnames(ddrepig)[c(1,2)]=c("bcr_patient_barcode","cancer")
#ddrepig$bcr_patient_barcode=as.character(ddrepig$bcr_patient_barcode)
#ddrepig$cancer=as.character(ddrepig$cancer)
#apply(ddrepig[,c(2:82)],2,function(x)sum(as.numeric(as.matrix(x[-which(is.na(x))]))))

#DDR genes
#ddr_genes=gene_lists$Nanostring_categories$DNA_repair#read.table("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/TCGA_DDR_Data_Resources/Genes.tsv",h=F)[,1]

### dependencies ###
source("../global_aes_out.R")
source("../dependency_files_tq.R")
source("../load_somatic_column.R")
rmhypermutator=TRUE

pathVarP=pathVarP[which(pathVarP$bcr_patient_barcode%in%overlapsplrmhyper),]
somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[which(somatic_likelyfunctional_driver$bcr_patient_barcode%in%overlapsplrmhyper),]

#germline and somatic  overlap
somatic_gene = unique(somatic_likelyfunctional_driver$Hugo_Symbol)
germline_gene = unique(pathVarP$HUGO_Symbol)

somatic_mut = unique(apply(somatic_likelyfunctional_driver,1,function(x)paste(c(x["Chromosome"],as.numeric(x["Start_Position"]),as.numeric(x["End_Position"]),x["Reference_Allele"],x["Tumor_Seq_Allele2"]),collapse = "|")))
germline_mut =unique(apply(pathVarP,1,function(x)paste(c(x["Chromosome"],as.numeric(x["Start"]),as.numeric(x["Stop"]),x["Reference"],x["Alternate"]),collapse = "|"))) 


library(VennDiagram)
library(grDevices)

par(mar=c(5,5,5,5))
tmp=venn.diagram(
  x = list(germline_gene,somatic_gene,ddrgene),
  category.names = c("" , "",""),
  filename =NULL, #'./out/Figure1A_germline_somatic_overlapgene.tiff',
  output = TRUE ,
  imagetype="tiff" ,
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('purple', 'green',"orange"),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"#,
  #cat.pos = c(-27,135),
  #cat.dist = c(0.055, 0.085),
  #rotation = 1
)

#pdf(file="./out/Figure1A_germline_somatic_overlapgene_rmhypermutator.pdf",height=2,width=2)
pdf(file="./out/Figure1A_germline_somatic_overlapgene.pdf",height=2,width=2)
grid.draw(tmp)
dev.off()

#mut 
library(VennDiagram)
library(grDevices)

par(mar=c(5,5,5,5))
tmp=venn.diagram(
  x = list(somatic_mut,germline_mut),
  category.names = c("" , ""),
  filename =NULL, # './out/Figure1A_germline_somatic_overlapMutations.tiff',
  output = TRUE ,
  imagetype="tiff" ,
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('orange', 'lightblue'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"#,
  #cat.pos = c(-27,135),
  #cat.dist = c(0.055, 0.085),
  #rotation = 1
)

#pdf(file="./out/Figure1A_germline_somatic_overlapMutations_rmhypermutator.pdf",height=2,width=2)
pdf(file="./out/Figure1A_germline_somatic_overlapMutations.pdf",height=2,width=2)
grid.draw(tmp)
dev.off()


# counts of somatic functional mutation by gene
somatic_gene_count = data.frame(table(somatic_likelyfunctional_driver$Hugo_Symbol))
germline_gene_count = data.frame(table(pathVarP$HUGO_Symbol))

colnames(somatic_gene_count) = c("Gene","PredictedFunctionalSomaticMutationCount")
colnames(germline_gene_count) = c("Gene","PathogenicGermlineVariantCount")
gene_count = merge(somatic_gene_count,germline_gene_count,by="Gene",all=T)
gene_count[is.na(gene_count)] = 0

highlight_g = as.character(gene_count$Gene[gene_count$PredictedFunctionalSomaticMutationCount > 400 | gene_count$PathogenicGermlineVariantCount > 10 | (gene_count$PredictedFunctionalSomaticMutationCount > 140 & gene_count$PathogenicGermlineVariantCount > 3)])

highlight_g=highlight_g[-which(highlight_g%in%c("EXT2","POT1","PRDM9","RECQL","COL7A1","GJB2"))]

#core DDR pathway
ddr=as.data.frame(read.csv("../../Huang_lab_data/TCGA_PanCanAtlas_2018/DDR_Knijnenburg_CellReport2018/DDR_Pathways.csv",h=T))
gene_count$GeneClass="Other genes"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Base.Excision.Repair..BER.))[-1]] = "Base Excision Repair"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Nucleotide.Excision.Repair..NER..including.TC.NER.and.GC.NER..))[-1]] = "Nucleotide Excision Repair"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Mismatch.Repair..MMR.))[-1]] = "Mismatch Repair"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Fanconi.Anemia..FA.))[-1]] = "Fanconi Anemia"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Homologous.Recomination..HR.))[-1]] = "Homologous Recomination"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Non.homologous.End.Joining..NHEJ.))[-1]] = "Nonhomologous End Joining"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Direct.Repair..DR.))[-1]] = "Direct Repair"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Translesion.Synthesis..TLS.))[-1]] = "Other genes"
gene_count$GeneClass[gene_count$Gene %in% unique(as.character(ddr$Damage.Sensor.etc.))[-1]] = "Damage Sensor"

gene_count$GeneClass

gene_count$GeneClass=factor(gene_count$GeneClass,levels=c("Homologous Recomination","Mismatch Repair","Nucleotide Excision Repair","Damage Sensor","Fanconi Anemia","Direct Repair","Translesion Synthesis","Other genes"))


colors = c("#ED2891","#C1A72F", "#FAD2D9","#F6B667","#97D1A9", "#B2509E", "#3953A4", "#007EB5")#,"#B2509E","#97D1A9","#ED1C24"
names(colors) =c("Homologous Recomination","Mismatch Repair","Nucleotide Excision Repair","Damage Sensor","Fanconi Anemia","Direct Repair","Translesion Synthesis","Other genes")


p = ggplot(gene_count,aes(x=PredictedFunctionalSomaticMutationCount, y =PathogenicGermlineVariantCount, color = GeneClass))
p = p + geom_point(stroke=0,alpha = 0.2) + theme_bw()
p = p + geom_text_repel(aes(label=ifelse(as.character(Gene) %in% highlight_g,as.character(Gene), NA)),cex=6,min.segment.length = 0)
p = p + theme(legend.position = c(0.74, 0.60),legend.text=element_text(size=16),legend.title=element_text(size=20),axis.title = element_text(size=20), axis.text.x = element_text(colour="black", size=20,vjust=0.5), axis.text.y = element_text(colour="black", size=20))#element_text(colour="black", size=14))
p = p +scale_x_log10() + scale_y_log10()
p = p + expand_limits(x = 0,y=0) + ylim(0,100)+ xlim(0,800)
p = p + xlab("Somatic Variant Count") + ylab("Germline Variant Count")
p = p + scale_color_manual("DDR Pathways",values =colors)
p
fn = "./out/somatic_vs_germline_var_counts_DDR_genes_allsample.pdf"
ggsave(file=fn, width=12, h =4, useDingbats=FALSE)


#######mutation frequency COADREAD-MIS-H, MSI-L####################
#Figure 1E Distribution of DDR mutated genes
#immuneProfile<-read.table("../germline_immune_cov/out/pca10260_immuneprofile_covariates.txt",h=T,sep="\t",stringsAsFactors=FALSE)
#immuneProfileMSI=immuneProfile[immuneProfile$variable=="MSISensor" & immuneProfile$TCGA_Study=="COADREAD",]
#msihspl=immuneProfileMSI$bcr_patient_barcode[which(immuneProfileMSI$value >=4)]
#msshspl=immuneProfileMSI$bcr_patient_barcode[which(immuneProfileMSI$value < 4)]

#clin$typemsi=ifelse(clin$bcr_patient_barcode%in%msihspl,"COADREAD-MSI",ifelse(clin$bcr_patient_barcode%in%msshspl,"COADREAD-MSS",clin$type))
#MSI<-data.frame(readxl::read_xlsx("/Users/qingtao/Box Sync/GermlineSomatic/Huang_lab_data/TCGA_PanCanAtlas_2018/MSI_Isidro_NatureCom2017/41467_2017_BFncomms15180_MOESM259_ESM.xlsx",sheet = "41467_2017_BFncomms15180_MOESM2"))
#MSI=MSI[MSI$Cancer_type%in%c("COAD","READ"),]

#tcgamsi=MSI$Barcode[MSI$MSI_category_nb_from_TCGA_consortium=="msi-h"]
#tcgamss=MSI$Barcode[MSI$MSI_category_nb_from_TCGA_consortium=="mss"]

#clin$typemsi=ifelse(clin$bcr_patient_barcode%in%tcgamsi,"COADREAD-MSI",ifelse(clin$bcr_patient_barcode%in%tcgamss,"COADREAD-MSS",ifelse(clin$type=="COADREAD","COADREAD-Other",clin$type)))

clin$typemsi=clin$type

allspl=table(clin$typemsi[clin$bcr_patient_barcode%in%overlapsplrmhyper])
#allspl=table(clin$typemsi)
freqMat=as.data.frame(cbind(Cancer=names(allspl),SampleSize=allspl))
freqMat$Cancer=as.character(freqMat$Cancer)
freqMat$SampleSize=as.numeric(as.matrix(freqMat$SampleSize))

#gfrequency
tmp0=pathVarP[which(pathVarP$HUGO_Symbol%in%ddrgene),]
if(any(duplicated(tmp0$bcr_patient_barcode))){
  tmp0=tmp0[-which(duplicated(tmp0$bcr_patient_barcode)),]
}
gfreq=table(tmp0$cancer)

#sfrequency
#epigstatus=ddrepig$bcr_patient_barcode[which(apply(ddrepig[,ddrgene],1,function(x)any(as.numeric(as.matrix(x))==1)))]

#ddrepigstatus=ddrepig[,c("bcr_patient_barcode","cancer")]
#ddrepigstatus$silencing=ifelse(ddrepigstatus$bcr_patient_barcode%in%epigstatus,"Yes","No")
somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$bcr_patient_barcode%in%overlapsplrmhyper,]

#somatic_likelyfunctional_driver=somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$bcr_patient_barcode%in%clin$bcr_patient_barcode,]

somatic_likelyfunctional_driver$cancer=clin$typemsi[sapply(somatic_likelyfunctional_driver$bcr_patient_barcode,function(x)which(clin$bcr_patient_barcode==x))] 
tmp1=somatic_likelyfunctional_driver[which(somatic_likelyfunctional_driver$Hugo_Symbol%in%ddrgene),]
#tmp2=rbind(ddrepig[ddrepig$bcr_patient_barcode%in%epigstatus,c("bcr_patient_barcode","cancer")],tmp1[,c("bcr_patient_barcode","cancer")])
#tmp2=tmp2[-which(duplicated(tmp2$bcr_patient_barcode)),]
#sfreq=table(tmp2$cancer)

tmp1=tmp1[-which(duplicated(tmp1$bcr_patient_barcode)),]
sfreq=table(tmp1$cancer)

#g+s frequency
mtmp0=tmp0[,c("bcr_patient_barcode","cancer")]
mtmp1=tmp1[,c("Tumor_Sample_Barcode","cancer")]
colnames(mtmp0)=colnames(mtmp1)=c("Tumor_Sample_Barcode","cancer")
mtmp=rbind(mtmp0,mtmp1)
#mtmp=mtmp[-which(duplicated(mtmp$Tumor_Sample_Barcode)),]
mfreq=table(mtmp$cancer)



freqMat$Germline=as.numeric((gfreq[freqMat$Cancer]/freqMat$SampleSize)*100)
freqMat$Somatic=as.numeric((sfreq[freqMat$Cancer]/freqMat$SampleSize)*100)
freqMat$Merge=as.numeric((mfreq[freqMat$Cancer]/freqMat$SampleSize)*100)


library(reshape2)
plotMat=melt(freqMat,id=c("Cancer","SampleSize"))
plotMat$Label=paste0(plotMat$Cancer," (",plotMat$SampleSize,")")
tmp3=plotMat[plotMat$variable=="Germline",]
or=tmp3$Label[order(tmp3$value,decreasing=T)]
plotMat$Label=factor(plotMat$Label,levels=or)
colnames(plotMat)=gsub("variable","Type",colnames(plotMat))

p= ggplot(plotMat,aes(y=as.numeric(as.matrix(value)),x=Label,fill=Type)) + geom_bar(stat="identity",position='dodge')
#p= p+ geom_text(aes(label=Counts), hjust=1, size=3)
p = p  + theme_bw() 
p = p  + theme(legend.position = c(0.85, 0.85),axis.text.x = element_text(colour="black", size=14, angle=90, vjust = 0.5,hjust = 0.95), axis.text.y = element_text(colour="black", size=14,hjust = 0.95),axis.ticks = element_blank(),plot.title = element_text(hjust = 0,size=16,face="bold"),axis.title=element_text(size=14,face="bold"),panel.border = element_blank(),axis.line= element_line(color='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p = p + labs(title="Samples affected by DDR mutations",y="Percentage",x="Cancer")
p

fn = "out/Figure1E_somaticgermline_frequency_of_DDR_affected_samples_MSI.pdf"
ggsave(file=fn, width=12, h =4, useDingbats=FALSE)


####################Correlated with Overall Response Rate###################
ORR<-data.frame(readxl::read_xlsx("./Clone_ORR.xlsx",sheet = "Sheet2"))
#ORR$types=gsub("COAD_MSI","COADREAD-MSI",gsub("COAD_MSS","COADREAD-MSS",ORR$types))
#ORR$types=gsub("COAD_MSI","COADREAD",gsub("COAD_MSS","COADREAD",ORR$types))
ORR$types=gsub("COAD","COADREAD",ORR$types)


ORR$ORR=ORR$ORR*100
ORR$gFreq=freqMat$Germline[pmatch(ORR$types,freqMat$Cancer)]
ORR$sFreq=freqMat$Somatic[pmatch(ORR$types,freqMat$Cancer)]
ORR$mFreq=freqMat$Merge[pmatch(ORR$types,freqMat$Cancer)]


plot(ORR$sFreq,ORR$ORR)

round(cor.test(ORR$sFreq,ORR$ORR)$p.value,digits = 2)


cc=round(cor(ORR$sFreq,ORR$ORR),digits = 2)
pvalue=round(cor.test(ORR$sFreq,ORR$ORR)$p.value,digits = 2)

p = ggplot(ORR,aes(y=ORR, x =sFreq))
p = p + geom_point(stroke=0,alpha = 0.2)+geom_smooth(method = lm)+geom_text(aes(label=types),cex=3) + theme_bw()#+xlim(0,20)+ylim(0,40)
p = p + geom_abline(intercept = 0, slope=1, alpha=0.2)
#p = p + geom_text_repel(aes(label=types))
p = p + theme(legend.position = c(0.74, 0.78),axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,vjust=0.5), axis.text.y = element_text(colour="black", size=14))
#p = p + scale_x_log10() + scale_y_log10()
#p = p + expand_limits(x = 0,y=0) + ylim(0,1100)
p=p+ggtitle(paste0())
p = p + xlab("Percentage of sample with DDR somatic mutation (%)") + ylab("Overall response rate (%)")
p=p+th
p
fn = "out/somatic_ddr_mutation_vs_overallresponserate.pdf"
ggsave(file=fn, width=5, h =4, useDingbats=FALSE)


plot(ORR$sFreq,ORR$ORR)
cor(ORR$sFreq,ORR$ORR)
cor.test(ORR$sFreq,ORR$ORR)
summary(lm(ORR~sFreq,data=ORR))

plot(ORR$gFreq,ORR$ORR)
cor(ORR$gFreq,ORR$ORR)
cor.test(ORR$gFreq,ORR$ORR)

plot(ORR$mFreq,ORR$ORR)
cor(ORR$mFreq,ORR$ORR)
cor.test(ORR$mFreq,ORR$ORR)

