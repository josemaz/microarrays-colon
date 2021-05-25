library(biomaRt)
library(TCGAbiolinks)
library(SummarizedExperiment)
# library(dplyr)
# library(magrittr)
library(EnhancedVolcano)
# BiocManager::install("biomaRt")

EXP.NAME <- "cesc"
WORKDIR <- "pipeline/TCGA"
dir.create(WORKDIR, recursive = TRUE)
raiz <- getwd()

setwd(WORKDIR)


######### Annotation
bm <- read.table(paste0(raiz,"/Data/biomart-20210507.tsv"), 
                 sep = '\t', header = TRUE)
rownames(bm) <- bm$stableID


######## DOWNLOAD DATA
# query <- GDCquery(project = "TCGA-CESC",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type = "HTSeq - Counts")
# GDCdownload(query)
# data.raw <- GDCprepare(query)
# save(data.raw, file=paste0(EXP.NAME,"-expr.RDATA"), compress="xz")
# load(file=EXP.NAME%+%"-expr.RDATA")
load(file=paste0(EXP.NAME,"-expr.RDATA"))


######### CHECK PAIRED SAMPLES OF CONTROL PATIENTS
designExp <- colData(data.raw)
cases.norm <- 
  unique(designExp$patient[designExp$sample_type=="Solid Tissue Normal"])
rows.qry <- 
  (designExp$patient %in% cases.norm) & (designExp$sample_type=="Primary Tumor")
cases.cancer <- designExp[rows.qry,]
cases.cancer <- unique(cases.cancer$patient)
stopifnot(setequal(cases.norm,cases.cancer))

cases <- designExp[designExp$patient %in% cases.norm,]
CESCrnaseqSE <- data.raw[,data.raw$barcode %in% cases$barcode]
# mat.expr=assay(CESCrnaseqSE)


######### ANALYZE
dataPrep <- TCGAanalyze_Preprocessing(CESCrnaseqSE, cor.cut = 0.6, 
                                      datatype = "HTSeq - Counts")
# BiocManager::install("EDASeq")
dataNorm<- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")
dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")
# # To see outliers
# boxplot(dataPrep, outline = FALSE)
# boxplot(dataNorm, outline = FALSE)
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.15)
filtro <- dataFilt[rowMeans(dataFilt)>10,]
threshold <- round(dim(filtro)[2]/2)
filtro <- filtro[rowSums(filtro == 0) <= threshold, ]
# boxplot(dataNorm, outline = FALSE)
# boxplot(dataFilt, outline = FALSE)
mat.normal <-  filtro[,cases[cases$sample_type == "Solid Tissue Normal",]$barcode]
mat.cancer <-  filtro[,cases[cases$sample_type == "Primary Tumor",]$barcode]
dataDEGs <- TCGAanalyze_DEA(mat1 = mat.normal ,
                            mat2 = mat.cancer ,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            # fdr.cut = 0.0 ,
                            # logFC.cut = 0.0,
                            method = "glmLRT")
dataDEGs$ID <- rownames(dataDEGs)
de <- merge(dataDEGs, bm, by="row.names")
stopifnot(nrow(dataDEGs)==nrow(de))
EnhancedVolcano(de,
                   title = 'TCGA-CESC-SQ-AD',
                   lab = de$gName,
                   x = 'logFC',
                   y = 'FDR',
                   # selectLab = slabs,
                   FCcutoff = 2.5,
                   pCutoff = 0.001,
                   # legendPosition = 'right',
                   # labSize = 3,
                   # legendLabSize = 8,
                   # legendIconSize = 3,
                   # # pointSize = 2,
                   # drawConnectors = TRUE,
                   # widthConnectors = 0.75,
                   # colConnectors = 'black',
                   # boxedLabels = TRUE,
                   ylim = c(0, -log10(dataDEGs[order(dataDEGs$FDR),][1,5]))
                     )
write.table(de, 
            file=paste0(raiz,"/Data/cervix-tumor-de.tsv"), 
            quote=FALSE, sep='\t',row.names=FALSE)


########################################
############ ONLY CANCER

# table(designExp$sample_type)
cases.cancer <- data.raw[,data.raw$sample_type == "Primary Tumor"]
stopifnot(length(unique(cases.cancer$submitter_id)) == 
            length(cases.cancer$submitter_id)) # all patients are unique

dataPrep <- TCGAanalyze_Preprocessing(cases.cancer, cor.cut = 0.6, 
                                      datatype = "HTSeq - Counts",
                                      filename = "cesc-cancer-cor.png")
dataNorm<- TCGAanalyze_Normalization(tabDF = dataPrep,
                                     geneInfo = geneInfoHT,
                                     method = "gcContent")
dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")
filtro <- dataNorm[rowMeans(dataNorm)>10,]
threshold <- round(dim(filtro)[2]/2)
filtro <- filtro[rowSums(filtro == 0) <= threshold, ]
# Quedan mas genes despues de la normalizacion que en el caso de las 3 muestras pareadas
# Annotations
annot <-subset(bm, stableID %in% rownames(filtro))
# stopifnot(nrow(annot) == nrow(filtro)) # Check all genes are annotated
filtro2 <- subset(filtro, rownames(filtro) %in% annot$stableID)
# Writing
write.table(data.frame("ensemble-IDs"=rownames(filtro2),filtro2), 
            file=paste0(raiz,"/Data/cervix-tumor-expr.tsv"), 
            quote=FALSE, sep='\t',row.names=FALSE)
write.table(annot, file=paste0(raiz,"/Data/cervix-tumor-annot.tsv"), 
            quote=FALSE, sep='\t',row.names=FALSE)

