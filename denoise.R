#PCA去除表达量中的技术噪音
#生物差异更能代表大部分的总体差异，这大部分的总体差异又是体现在前几个PCs
#技术差异也有影响，但是会体现在后面的几个主成分中
#因此可以利用denoisePCA去除,原理是去掉主成分中方差之和等于var.out中技术因素的和
sce <- denoisePCA(sce, technical = var.out, assay.type = "corrected")

#可视化PCA
plotReducedDim(sce, dimred="PCA",ncomponents = 3,colour_by = "Oncogene")
#看看去除批次后的效果，如果仍然分开，则plate差异仍然很大，会影响生物学差异
plotReducedDim(sce, dimred="PCA",ncomponents = 3,colour_by = "Plate")

#可视化t-SNE,基于在邻域图随机游走的概率分布，因为tSNE能在高维空间直接捕获细胞非线性关系
#需要设置随机种子，同时需要设置perplexity参数，会影响最后显示的细胞分类状况
#perplexity一般在5-50之间
#perplexity:the number of close neighbours each point has
set.seed(100)
sce <- runTSNE(sce,dimred="PCA",perplexity = 20)
plotTSNE(sce, colour_by="Oncogene")
