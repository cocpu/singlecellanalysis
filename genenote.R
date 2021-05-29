#获取表达矩阵中的ERCC信息 , grep1返回逻辑值
#spike-in是已知浓度的外源RNA，在单细胞裂解后加入再进行反转录，可以更好地估计和消除系统误差
is.spike <- grepl("^ERCC", rownames(sce))#检查基因名称中是否包含ERCC
#summary(is.spike)
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC","gene"))#将spike-in分离出来

#删除掉SIRV信息
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,]
summary(is.sirv)

#读取细胞注释信息metadata，并将其添加到colData中
metadata <- read.delim(file.path('C:/Users/26068/Desktop/github/singlecellanalysis',
                                 "E-MTAB-5522.sdrf.txt"),
                       check.names = FALSE, header = TRUE)
#match的作用是返回前一个参数在后一个参数的位置
m <- match(colnames(counts(sce)), metadata[["Source Name"]])
stopifnot(all(!is.na(m))) #检查是否完整包含所有细胞
metadata <- metadata[m,]  #选择和表达矩阵相关的细胞注释
head(colnames(metadata))

#添加细胞注释到对象
#1、细胞来源板
colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
#2、细胞表型
pheno <- metadata[["Factor Value[phenotype]"]]
colData(sce)$Oncogene <- pheno
levels(colData(sce)$Oncogene) <-c("induced", "control")

table(colData(sce)$Oncogene, colData(sce)$Plate)

#基因ID转换
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys = rownames(sce),
               keytype="ENSEMBL", column="SYMBOL")
head(symb)
identical(names(symb),rownames(sce))
rowData(sce)$SYMBOL <- symb
rowData(sce)$ENSEMBL <- rownames(sce)
head(rowData(sce))
#整合行名
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL,
                                      rowData(sce)$SYMBOL)
head(rownames(sce))
#添加基因在染色体位置信息
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
columns(TxDb.Mmusculus.UCSC.mm10.ensGene)
head(transcripts(TxDb.Mmusculus.UCSC.mm10.ensGene,
                 columns=c('CDSCHROM')))
location <-mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene,
                             keys=rowData(sce)$ENSEMBL, 
                             column ="CDSCHROM",
                             keytype = "GENEID")
rowData(sce)$CHR <-location
summary(location=="chrM")