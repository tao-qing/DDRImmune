##### dependency_files.R #####
# Kuan-lin Huang @ WashU 2017 June
# updated 2018 to clean up paths
# dependent files for analysis in the PCA Germline project

# the pathogenic variant file: please swap to your local secured path don't share on Box
fn = "../../Huang_lab_data/TCGA_PanCanAtlas_2018/PCA_Germline_Cell2018/PCA_pathVar_integrated_filtered_adjusted_wSomaticCounts.tsv"
fn = "~/Box\ Sync/Huang_lab_kuan/control_access_data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv"
pathVar = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)


# gene lists
volg_fn = "../../reference_files/Volgestin2013Science_125genes_class.txt"
volg_class = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=volg_fn)
colnames(volg_class) = c("gene_name","Gene_Classification","Pathway","Process")
volg_TSGs = volg_class$gene_name[volg_class$Gene_Classification=="TSG"]
volg_oncogenes = volg_class$gene_name[volg_class$Gene_Classification=="Oncogene"]

onco_fn = "../../reference_files/GSEA_geneLists/oncogenes.txt"
oncogenes = as.vector(t(read.table(header=F, stringsAsFactors = F, file=onco_fn)))

tsg_fn = "../../reference_files/GSEA_geneLists/tumor_suppressors.txt"
TSGs = as.vector(t(read.table(header=F, stringsAsFactors = F, file=tsg_fn)))

additional_TSGs = c("MAX","SDHA","ATR","BARD1","ERCC1","FANCI","FANCL","FANCM","POLD1","POLE","POLH","RAD50","RAD51","RAD51C","RAD51D","RAD54L")
others = c("ABCB11","ABCC2","AXIN2","CBL","CDKN1B","COL7A1","CYP17A1","CYP1B1","DIS3L2",
           "DKC1","DOCK8","ELANE","FAH","FLCN","GBA","GJB2","HFE","HMBS","LRRK2",
           "MAX","MTAP","MYD88"   
           ,"PRSS1","RHBDF2","RPL5","SDHA","SETBP1","SF3B1","SH2D1A","SLC25A13","SOS1","TMEM127","TRIM37","UROD")

additional_oncogenes = c("AR","STAT3","TERT","MAP2K2","NOTCH3")

all_oncogenes = c(oncogenes,volg_oncogenes,additional_oncogenes)
all_TSGs = c(TSGs,volg_TSGs,additional_TSGs)

all_oncogenes = all_oncogenes[-which(all_oncogenes=="SETBP1")]

gene_fn = "../../reference_files/PCA_feature_gene_list.txt"
glist_f = read.table(header=FALSE, stringsAsFactors = F, file = gene_fn)
featGenes = as.vector(t(glist_f))



# # use LOH_Classification instead, and nominate anything with tumorVAF > 0.6 but not normalVAF as suggestive
# pathVar$LOH_classification = pathVar$LOH_Sig
# pathVar$LOH_classification[pathVar$tumorVAF > 0.6 & pathVar$normalVAF < 0.6 & pathVar$LOH_Sig != "Significant"] = "Suggestive"
# 
# pathVar$AS_LOH[pathVar$LOH_classification %in% c("Suggestive","Significant")] = "Unclassified LOH"
# pathVar$AS_LOH[pathVar$LOH_classification %in% c("Suggestive","Significant") & !is.na(pathVar$CNV_thres) & pathVar$CNV_thres < 0 &
#                         pathVar$tumorVAF > pathVar$normalVAF] = "Copy Number Deletion of WT Allele"
# pathVar$AS_LOH[pathVar$LOH_classification %in% c("Suggestive","Significant") & !is.na(pathVar$CNV_thres) & pathVar$CNV_thres > 0 &
#                         pathVar$tumorVAF > pathVar$normalVAF] = "Copy Number Amplification of Variant Allele"
# 
# fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/data_integration/out/PCA_pathVar_integrated_filtered_adjusted.tsv"
# write.table(pathVar, file=fn, quote=F, sep="\t", col.names=T, row.names=F)

# fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/data_integration/out/PCA_pathVar_integrated_filtered_adjusted_path.tsv"
# write.table(pathVar[pathVar$Overall_Classification %in% c("Pathogenic","Likely Pathogenic"),], file=fn, quote=F, sep="\t", col.names=T, row.names=F)
# fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/data_integration/out/PCA_pathVar_integrated_filtered_adjusted_pVUS.tsv"
# write.table(pathVar[pathVar$Overall_Classification %in% c("Prioritized VUS"),], file=fn, quote=F, sep="\t", col.names=T, row.names=F)

##### subsets #####
pathVarOT = pathVar[!is.na(pathVar$Gene_Classification) & pathVar$Gene_Classification != "None",]

truncations = pathVarOT[pathVarOT$binary_type=="Truncation",]
missenses = pathVarOT[pathVarOT$binary_type=="Missense",]

pathVarFGene = pathVar[pathVar$HUGO_Symbol %in% featGenes,]

pathVarP = pathVar[pathVar$Overall_Classification %in% c("Pathogenic","Likely Pathogenic"),]

pathVarPOT = pathVarP[!is.na(pathVarP$Gene_Classification) & pathVarP$Gene_Classification != "None",]

PCA_count = data.frame(table(pathVarP$HUGO_Symbol))
colnames(PCA_count) = c("Gene","Count")
gene_order = PCA_count$Gene[order(PCA_count$Count,decreasing = T)]
##### clinical files #####
clin_f = "../../Huang_lab_data/TCGA_PanCanAtlas_2018/clinical/PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)

PCs_f = "../../Huang_lab_data/TCGA_PanCanAtlas_2018/PCA_Germline_Cell2018/TCGA_ancestry_PC.txt"
PCs = read.table(header=T, quote = "", sep="\t", fill =T, file = PCs_f, stringsAsFactors=FALSE)
