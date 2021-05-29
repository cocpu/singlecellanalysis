#检查高表达基因，可以检验是否在文库制备、比对环节出现问题
#调整字体大小
fontsize <- theme(axis.text = element_text(size=12), axis.title = element_text(size = 16))
plotHighestExprs(sce, n = 50)+fontsize

#过滤低丰度基因，即表达量几乎为0，在下游分析中起不到统计作用的基因
#在假设检验中，它们不能提供足够的证据去推翻原假设，还会增加多重检验的计算量
#并且离散的counts也可能干扰统计过程
#我们需要根据每一步的需要去过滤
ave.counts <- calculateAverage(sce)#计算单个基因的平均表达量
hist(log10(ave.counts),breaks = 100, col = "grey80",
     xlab = expression(Log[10]~"average count"))

#1、根据均值进行过滤
demo.keep <- ave.counts >= 1 #一般过滤小于1的基因（log10 <0)
summary(demo.keep)
sce <- sce[demo.keep,]

#2、根据每个基因在多少个细胞中有表达进行过滤
num.cell <- nexprs(sce, byrow = TRUE) #计算非0的行
smoothScatter(log10(ave.counts), num.cell, ylab = "Number of cells",
              xlab = expression(Log[10]~"average count"))
to.keep <-num.cell > 0
summary(to.keep)
sce <- sce[to.keep, ]
