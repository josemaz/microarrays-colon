library(crayon)
library(dplyr)
library(gprofiler2)
library(edgeR)
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
library(oligo)
library(stringr)

dat <- read.celfiles(paste0('Data/',list.celfiles('Data')))
eset <- rma(dat)
annotation(eset)
dim(eset)

library(affycoretools)
eset.main <- getMainProbes(eset)
eset.main <- annotateEset(eset.main, pd.hta.2.0)
gene.symbols <- fData(eset.main)$SYMBOL # if you want data original

M <- exprs(eset.main)
rownames(M) <- gene.symbols
cat("NA gene symbol: ",sum(is.na(rownames(M))))
M <- M[!is.na(rownames(M)),]
head(M)
cat(red("Genes duplicados: ",sum(duplicated(rownames(M))),"\n"))
# # which(duplicated(rownames(M)))
M <- M[!duplicated(rownames(M)),]
cat(red("Numero de genes: ",nrow(M),"\n"))

nfile <- sampleNames(dat)
snames <- sampleNames(eset.main)
tipo <- factor(substr(snames, 1, 2))
tiempo <- factor(substr(snames, 5, 6))
replica <- substr(snames,3,3)
muestra <- paste0(tipo,"-",tiempo)
muestras <- data.frame(nfile, muestra, tipo, tiempo, replica)
muestras <- muestras[order(muestras$muestra),]

deg <- function(m,mdesign, slabs=NULL, tit=""){
  cat(red("Ajuste ...\n"))
  fit <- lmFit(m, mdesign)
  head(fit$coefficients)
  fit <- eBayes(fit)
  topTable(fit,coef=2,genelist = rownames(m))
  dim(fit)
  colnames(fit)
  rownames(fit)[1:10]
  names(fit)
  
  # cat(red("Treat ...\n"))
  # fit2 <- treat(fit,lfc=1.5)
  # topTreat(fit2, coef=2, sort.by="logFC")
  
  cat(red("Toptable ...\n"))
  # volcanoplot(fit,coef=2,highlight=3)
  FCthres = 0.75
  pCut = 0.01
  x <- topTable(fit, coef=2, adjust="fdr", sort.by="P",number=nrow(mtemp),
                genelist = rownames(m))
  if(is.null(slabs)){
    print(head(x))
  }else{
    print(x[x$ID %in% slabs,])
  }
  cat(red("Volcano ...\n"))
  v <- EnhancedVolcano(x,
                  lab = x$ID,
                  # title = "",
                  subtitle = tit,
                  x = 'logFC',
                  y = 'P.Value',
                  # selectLab = slabs,
                  FCcutoff = FCthres,
                  pCutoff = pCut,
                  legendPosition = 'right',
                  labSize = 3,
                  # labCol = 'black',
                  # labFace = 'bold',
                  legendLabSize = 8,
                  legendIconSize = 4,
                  # pointSize = 2,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75,
                  colConnectors = 'black',
                  boxedLabels = TRUE,
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  col=c('black', 'black', 'black', 'red3'),
                  colAlpha = 1,
                  ylim = c(0, -log10(x[order(x$P.Value),][1,5])))
  cat(red("GO ...\n"))
  geneGO <- x %>%
    filter((abs(x$logFC)>FCthres) & (x$P.Value< pCut)) #%>%
    # distinct(ID)
  cat(paste(unlist(geneGO$ID), collapse=' '),"\n")

  gostres <- gost(query = geneGO$ID,
                  organism = "hsapiens", significant = TRUE,
                  user_threshold = 0.05, correction_method = "g_SCS",
                  sources = c("GO:BP","GO:MF","GO:CC"),
                  domain_scope = "annotated")
  if(!is.null(gostres$result)) print(gostres$result)
  return(v)
}


################## 2s:4s
cat(red("Contraste Vh vs Tx for 2s ...\n"))
filtro <- muestras %>% 
  filter(tiempo == "2s") %>%          # can be t
  arrange(factor(tipo, levels = c("Vh","Tx")))
mtemp <- M[,filtro$nfile]
# design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# design <- model.matrix(~ filtro$tipo)
print(filtro)
print(head(mtemp))
print(design)

options(ggrepel.max.overlaps = Inf)
v <- deg(mtemp,design, tit = "11 days")
v
v <- deg(mtemp,design,c('CASTOR2','PCDH17','COL1A2'))
v 

cat(red("Contraste Vh vs Tx for 4s ...\n"))
filtro <- muestras %>% 
  filter(tiempo == "4s") %>%          # can be t
  arrange(factor(tipo, levels = c("Vh","Tx")))
mtemp <- M[,filtro$nfile]
# design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# design <- model.matrix(~ filtro$tipo)
print(filtro)
print(head(mtemp))
print(design)

v <- deg(mtemp,design, tit = "24 days")
v
v <- deg(mtemp,design,c('NR4A1','NR4A3','ATF3','COL1A2'))
v 


################## Tx:Vh
cat(red("Contraste 2s vs 4s for Vh:Vh ...\n"))
filtro <- muestras %>%
  filter(tipo == "Vh") %>%          # can be t
  arrange(factor(tiempo, levels = c("2s","4s")))
mtemp <- M[,filtro$nfile]
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# design <- model.matrix(~ filtro$tiempo)
print(filtro)
print(head(mtemp))
print(design)

v <- deg(mtemp,design)
v


cat(red("Contraste 2s vs 4s for Tx:Tx ...\n"))
filtro <- muestras %>%
  filter(tipo == "Tx") %>%          # can be t
  arrange(factor(tiempo, levels = c("2s","4s")))
mtemp <- M[,filtro$nfile]
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# design <- model.matrix(~ filtro$tiempo)
print(filtro)
print(head(mtemp))
print(design)

v <- deg(mtemp,design)
v









# lp <- -log10(fit$p.value[,1])
# ord <- order(lp, decreasing = TRUE)[1:10]
# x <- fit$coef[ord,1]
# y <- lp[ord]
# text(x, y, names(fit$sigma)[ord], cex=0.8, col="red")
# 
# x <- topTable(fit2, coef=2, adjust="fdr", sort.by="P")
# y <- x[x$adj.P.Val < 0.01 & (x$logFC > 0.5 | x$logFC < -1) & x$AveExpr > 10,]; 
# y; print("Number of genes in this list:"); 
# length(y$ID)
# 
# gesn 

