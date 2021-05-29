##计算质控指标
mito <- which(rowData(sce)$CHR=="chrM")
sce<-addPerCellQC(sce,subsets = list(Mt = mito))
sce<-addPerFeatureQC(sce)
#原版：sce <- calculateQCMetrics(sce, feature_controls =list(Mt = mito))
head(colnames(colData(sce)),10)

#可视化
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
#multiplot(
#  plotColData(sce, y="sum", x="PlateOnco"),
# plotColData(sce, y="detected", x="PlateOnco"),
#  plotColData(sce, y="altexps_ERCC_sum", x="PlateOnco"),
#  plotColData(sce, y="subsets_Mt_sum", x="PlateOnco"),
#  cols = 2)

#文库大小，基因表达以及ERCC转录本或线粒体基因的比例
par(mfrow = c(1,3))
plot(sce$detected, sce$sum/1e6, 
     xlab="number of expressed gene",
     ylab="library size (millions")
plot(sce$detected, sce$altexps_ERCC_percent, 
     xlab="number of expressed gene",
     ylab="ERCC proportion(%)")
plot(sce$detected, sce$subsets_Mt_percent, 
     xlab="number of expressed gene",
     ylab="Mitochondrial proportion(%)")

##鉴定离群点
#过滤低文库大小、低表达基因、高spike-in
#（以低于中位数3倍MAD作为过滤条件）
libsize.drop <- isOutlier(sce$sum, nmads = 3, 
                          type = "lower",log = TRUE,
                          batch = sce$PlateOnco)
feature.drop <- isOutlier(sce$detected, nmads = 3, 
                          type = "lower",log = TRUE,
                          batch = sce$PlateOnco)
spike.drop <- isOutlier(sce$altexps_ERCC_percent, nmads = 3, 
                          type = "higher",log = TRUE,
                          batch = sce$PlateOnco)
keep <- !(libsize.drop|feature.drop|spike.drop)
data.frame(Bylibsize=sum(libsize.drop), Byfeature=sum(feature.drop),
           BySpike = sum(spike.drop), remaining=sum(keep))

sce$PassQC <- keep
#saveRDS(sce, file="4sce16B_preQC.rds")
sce <- sce[,keep]
dim(sce)


#PCA离群值检测
#sce_tmp <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE)
#table(sce_tmp$outlier)