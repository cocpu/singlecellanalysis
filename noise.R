#模拟基因表达中的技术噪音
var.out <-modelGeneVarWithSpikes(sce,"ERCC",block = sce$Plate)#这个函数与下一个函数不同在于根据spike拟合技术误差曲线
var.out.spike <- modelGeneVar(altExp(sce,"ERCC"),block=sce$Plate)
#内源基因的总方差随表达量的变化，总方差包含生物因素和技术因素
plot(var.out$mean, var.out$total, pch = 16, cex = 0.6,
     xlab = "mean log-expression", ylab="variance of log-expression")
var.fit.spike <- fitTrendVar(rowMeans(logcounts(altExp(sce,"ERCC"))), rowVars(logcounts(altExp(sce,"ERCC"))))
#蓝线：根据spike-in表达量拟合得到的曲线，可理解成技术因素的方差随表达量的变化
curve(var.fit.spike$trend(x), col="dodgerblue", lwd = 2, add = TRUE)
#红点：Spike-in表达量的均值和方差散点
points(var.out.spike$mean, var.out.spike$total,col="red",pch = 16)

#为了证明不是离群点，而是一种细胞造成的基因表达差异，绘图观察
chosen.gene <- order(var.out$bio, decreasing = TRUE)[1:10]
plotExpression(sce, features = rownames(var.out)[chosen.gene]) + fontsize

#去除批次效应
library(limma)
assay(sce,"corrected") <- removeBatchEffect(logcounts(sce),
                                            design = model.matrix(~sce$Oncogene), batch = sce$Plate)
#可以改善下游的降维、聚类
#当每个批次中的细胞组成类型已知或不同批次都是同一种类型的细胞只是处理方式不同，这种情况
#用removeBatchEffect效果不错
#但大部分的scRNA数据中，每个批次的细胞类型我们可能不确定，这种情况下使用mnnCorrect()