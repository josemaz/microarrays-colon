library(crayon)
library(dplyr)
library(gprofiler2)
library(edgeR)
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
write.table(M,file="Data/expr.txt")

##################### PYTHON
cat("NA gene symbol: ",sum(is.na(rownames(M))))
M <- M[!is.na(rownames(M)),]
nrow(M); head(M)
cat(red("Genes duplicados: ",sum(duplicated(rownames(M))),"\n"))
# # which(duplicated(rownames(M)))
M <- M[!duplicated(rownames(M)),]
cat(red("Numero de genes: ",nrow(M),"\n"))
##################### PYTHON


nfile <- sampleNames(dat)
snames <- sampleNames(eset.main)
tipo <- factor(substr(snames, 1, 2))
tiempo <- factor(substr(snames, 5, 6))
replica <- substr(snames,3,3)
muestra <- paste0(tipo,"-",tiempo)
muestras <- data.frame(nfile, muestra, tipo, tiempo, replica)
muestras <- muestras[order(muestras$muestra),]

cat(red("Contraste Vh vs Tx for 4s:2s ...\n"))
filtro <- muestras %>% 
  filter(tiempo == "2s") %>%          # can be t
  arrange(factor(tipo, levels = c("Vh","Tx")))
mtemp <- M[,filtro$nfile]
# design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# design <- model.matrix(~ filtro$tipo)


cat(red("Ajuste ...\n"))
fit <- lmFit(mtemp, design)
head(fit$coefficients)
fit <- eBayes(fit)
topTable(fit,coef=2,genelist = rownames(mtemp))
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

# cat(red("Treat ...\n"))
# fit2 <- treat(fit,lfc=1.5)
# topTreat(fit2, coef=2, sort.by="logFC")

cat(red("Volcano ...\n"))
# volcanoplot(fit,coef=2,highlight=3)
FCthres = 0.75
pCut = 0.01
x <- topTable(fit, coef=2, adjust="fdr", sort.by="P",number=nrow(mtemp),
              genelist = rownames(mtemp))
v <- EnhancedVolcano(x,
                     lab = x$ID,
                     x = 'logFC',
                     y = 'P.Value',
                     selectLab = c('CASTOR2','PCDH17','COL1A2'),
                     FCcutoff = FCthres,
                     pCutoff = pCut,
                     legendPosition = 'right',
                     legendLabSize = 3,
                     legendIconSize = 3,
                     drawConnectors = TRUE,
                     widthConnectors = 0.75,
                     ylim = c(0, -log10(x[order(x$P.Value),][1,5])))
v
