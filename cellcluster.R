#细胞聚类：利用去除噪音的PCA结果进行细胞聚类，推断细胞亚群
pcs<-reducedDim(sce,"PCA")[,1:2]
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree,distM=as.matrix(my.dist),
                                    minClusterSize = 10, verbose=0))
table(my.clusters, sce$Plate)#观察批次效应
table(my.clusters, sce$Oncogene)#观察生物差异
sce$cluster <- factor(my.clusters)
#观察结果
plotTSNE(sce, colour_by = "cluster") + fontsize
#检查分群结果：silhouette width轮廓图
#一个细胞的横坐标width越大,表示相比其他cluster细胞更接近这个cluster
library(cluster)
clust.col <- scater:::.get_palette("tableau10mediu")
sil <- silhouette(my.clusters, dist=my.dist)
sil.cols<-clust.col[ifelse(sil[,3]>0, sil[,1], sil[,2])]
sil.cols<-sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)),"clusters"),
     border = sil.cols, col = sil.cols, do.col.sort = FALSE)

#检测marker基因
#两两cluster间对每个基因的log表达量进行Weltch t检验
markers <- findMarkers(sce, my.clusters, block = sce$Plate)#包含排列好的基因
marker.set1 <- markers$"1"
top.markers <- rownames(marker.set1)[marker.set1$Top <= 10]
plotHeatmap(sce, features = top.markers, columns = order(sce$cluster),
            colour_columns_by = c("cluster","Plate","Oncogene"),
            cluster_cols = FALSE, center = TRUE, symmetric = TRUE,
            zlim = c(-5,5), show_colnames = F)

