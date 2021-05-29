plate1 <- read.delim(file.path('C:/Users/26068/Desktop/github/singlecellanalysis',
                               "counts_Calero_20160113.tsv"),
                     header = TRUE, row.names = 1, check.names = FALSE) #check-name:防止修改-符号
plate2 <- read.delim(file.path('C:/Users/26068/Desktop/github/singlecellanalysis', 
                               "counts_Calero_20160325.tsv"),
                     header = TRUE, row.names = 1, check.names = FALSE)
gene.lengths <- plate1$Length #第一列为基因长度
plate1 <- as.matrix(plate1[,-1]) #将第一列舍去
plate2 <- as.matrix(plate2[,-1])
stopifnot(identical(rownames(plate1), rownames(plate2))) #检查行名是否一致
all.counts <- cbind(plate1, plate2) #合并两个样本
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts = all.counts)) #生成一个S4对象，这个对象可理解为由多个dataframe组成
rowData(sce)$GeneLength <- gene.lengths
