#利用一个做好的训练集，训练数据集中已经计算好两两基因的差异，并将属于不同细胞周期
#且存在差异的基因对作为一个marker pair，然后就在已知表达矩阵中对每个细胞测试
#这些marker pairs与训练数据集的相似程度，最后得到相似程度分值，对细胞进行归类
library(scran)#下载scran包后，在电脑系统文件中会有exdata
mm.pairs <- readRDS(system.file("exdata","mouse_cycle_markers.rds", 
                                package = "scran"))
head(mm.pairs$G1)#反映G1期的基因对
assignments <- cyclone(sce, mm.pairs, gene.names = rowData(sce)$ENSEMBL)
plot(assignments$score$G1, assignments$score$G2M,
     xlab = "G1 score", ylab = "G2/M score", pch = 16)

#细胞周期划分规则
#如果一个细胞再G1中得分大于0.5，并且高于G2M得分，为G1
#如果G2M得分大于0.5且得分高于G1，为G2M
#如果G1和G2M都小于0.5，则为S
colData(sce)$phases <- assignments$phases

#在细胞周期判断前不能过滤低丰度转录本，否则无法利用基因对进行判断