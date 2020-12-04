library(edgeR)
library(oligo)

dat <- read.celfiles(list.celfiles())
eset <- rma(dat)
annotation(eset)
dim(eset)
# con <- db(pd.hta.2.0)
# # z <- dbGetQuery(con, "select distinct core_mps.transcript_cluster_id, type from featureSet inner join core_mps")
# z <- dbGetQuery(con, "select distinct core_mps.transcript_cluster_id, type from featureSet inner join core_mps using(fsetid);")
# table(z$type)
# dim(getMainProbes(eset))
# reps <- z[duplicated(z[,1]) & !is.na(z$type),1]
# table(z[z[,1] %in% reps,2], useNA = "ifany")


library(affycoretools)
# library(hta20transcriptcluster.db)

# eset.main <- getMainProbes(eset, pd.hta.2.0)
eset.main <- getMainProbes(eset)
# eset.main <- annotateEset(eset.main, hta20transcriptcluster.db)
eset.main <- annotateEset(eset.main, pd.hta.2.0)
# eset.main <- annotateEset(eset.main, hta20probeset.db)
# eset.main <- annotateEset(eset.main, pd.hta.2.0, level = "probeset")

# featureNames(eset.main)[1:5]
# snames <- sampleNames(eset.main)
M <- exprs(eset.main)

# nfile <- sampleNames(dat)
# snames <- sampleNames(eset.main)
# replica <- substr(snames,3,3)
# muestra <- paste0(substr(snames, 1, 2),substr(snames, 4, 6))
# muestras <- data.frame(nfile, muestra, replica)
# muestras <- muestras[order(muestras$muestra),]
# DESIGN MATRIX
# group <- factor(muestras$muestra)
#Con intercepto
# dmatrix <- model.matrix(~ group)
#Sin intercepto
# dmatrix <-  model.matrix(~ group + 0)
M <- M[,c(-2,-4,-6,-8,-10,-12)]
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))


## AJUSTE ##
fit <- lmFit(M, design)
head(fit$coefficients)
fit <- eBayes(fit)
topTable(fit,coef=2)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

fit2 <- treat(fit,lfc=1.5)
topTreat(fit2,coef=2)

# Volcano plot
volcanoplot(fit,coef=2,highlight=2)






