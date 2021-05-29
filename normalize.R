#标准化normalization：normalize within the cell for the difference in sequencing depth/mRNA thruput
#即对文库进行处理，为了消除一些技术差异,比如细胞捕获效率，测序深度
#scaling:为了消除极值带来的影响

#size factor:represent the extent to which counts should be scaled in each library
#常规转录组也有size factor，但是由于单细胞转录组有太多低丰度，单个细胞计算误差大
#因此：首先计算所有细胞中每个基因的几何平均值，每个细胞的size factor
#是基因表达与基因几何平均值比值的中位数

#1\去卷积方法：根据内源基因计算size factor
#假设大部分基因在细胞间并非差异表达
#先将细胞couns混合起来，计算混合的sf，在去卷积，还原得到单细胞SF
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))
plot(sce$sum/1e6, sizeFactors(sce), log="xy", xlab = "Library size(millions)",
     ylab = "size factor", col = c("red", "black"), pch = 16)
legend("bottomright", col = c("red","black"), pch =16, cex=1.2,
       legend = levels(sce$Oncogene))
#SF与文库大小相关性好，说明大部分文库大小差异是由技术因素造成

#对于高度异质的细胞，为了保证假设，需要quickCluster()做一个聚类，再将clusters参数传递给
#computerSumFactors计算

#2\根据spike-in转录本计算size factor
#单细胞测序，每个细胞的文库大小不一样，但是spike-in的量是一样的
#sF目的是消除技术误差，但是spike-in与内源基因不同，因此不能用上面的方法计算spike-in的SF
sce <- computeSpikeFactors(sce,"ERCC")

#这个函数可以综合以上的内容
sce <- logNormCounts(sce, use.altexps = "ERCC")
#将每个细胞中基因的表达量除以这个细胞计算的size FACTOR,然后log(x+1)处理，方便展示数据全貌