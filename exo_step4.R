
##WGCNA分析不同riskscore的模块 进行分析  GO 三套
library(WGCNA)
library(NMF)
library(MOVICS)

#data$ID <- rownames(data)
##C1 C2
##使用前5000个gene
colnames(KIRC_tpm) <- str_sub(colnames(KIRC_tpm),1,15)
KIRC_tpm[1:4,1:4]
KIRC_tpm <- log2(KIRC_tpm+1)
mydata <- KIRC_tpm[mRNA_anno$gene_name,sub$Sample]
dim(mydata)
#此处用head(mydata)看出行名是基因名，第1列至第20列是基因在20个sample里的表达量
head(mydata)
mydata[1:4,1:4]
##转置
datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr1<-data.matrix(expro.upper)
##检查样本质量
gsg = goodSamplesGenes(datExpr1, verbose = 3)
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 35000”和“cutHeight = 35000”中的“500”换成你的cutoff
  abline(h = 90000, col = "red") 

clust = cutreeStatic(sampleTree, cutHeight = 90000, minSize = 10)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
#运行下面的代码，用全部20个样本进行后续的分析
#datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##选择阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
#构建网络，找出gene module
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)

mergedColors = labels2colors(net$colors)

pdf("2module.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}

##载入表型数据
colnames(sub)
samples=sub[,c("group","gpr_score")]

datTraits <- sub[,c("group","gpr_score")]
datTraits$group <- factor(datTraits$group,levels = c("Low","High"))

### 注意一点。。。
design=model.matrix(~0+ datTraits$group)
design <- as.data.frame(design)
colnames(design)=levels(datTraits$group)
design
design$risk_score <- as.numeric(datTraits$gpr_score)
rownames(design) <- rownames(sub)
samples <- design
dim(MEs)
samples <- samples[rownames(MEs),]

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)

modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()

###选择模块进行分析  GO KEGG 分析
#yellow模块进行分析
gene_yellow <- fread("~/TCGA/exosome_ccRCC/yellow.csv",header = T,data.table = F)
gene_yellow[1:4]

module = "yellow"
moduleGenes = moduleColors==module
datExpr <- as.data.frame(datExpr)
dd <- datExpr[,moduleGenes]
dd[1:4,1:4]

#相关R包载入（还有部分R包使用前再载入）：
library(dplyr)#数据清洗
library(stringr)
library(org.Hs.eg.db)#物种注释(Homo sapiens)
library(clusterProfiler)#富集分析
library(ggplot2)#个性化绘图
library(RColorBrewer)#配色调整
#根据研究目的筛选差异基因(仅上调、下调或者全部)：
yellow_gene <-colnames(dd)#所有差异基因
#ID转换：
#查看可转换的ID类型：
columns(org.Hs.eg.db)

##使用clusterProfiler包自带ID转换函数bitr(基于org.Hs.eg.db)：
#diff：
diff_entrez <- bitr(yellow_gene,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
head(diff_entrez)
#GO(Gene Ontology)富集分析：
##MF(我们以总差异基因的GO富集为例):
GO_MF_diff <- enrichGO(gene = diff_entrez$ENTREZID, #用来富集的差异基因
                       OrgDb = org.Hs.eg.db, #指定包含该物种注释信息的org包
                       ont = "MF", #可以三选一分别富集,或者"ALL"合并
                       pAdjustMethod = "BH", #多重假设检验矫正方法
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE) #是否将gene ID映射到gene name
#提取结果表格：
GO_MF_result <- GO_MF_diff@result
#View(GO_MF_result)
#CC:
GO_CC_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
#提取结果表格：
GO_CC_result <- GO_CC_diff@result

#BP:
GO_BP_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
#提取结果表格：
GO_BP_result <- GO_BP_diff@result

#MF、CC、BP三合一：
GO_all_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "ALL", #三合一选择“ALL”
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
#提取结果表格：
GO_all_result <- GO_all_diff@result

##保存GO富集结果：
save(GO_MF_diff,GO_CC_diff,GO_BP_diff,GO_all_diff,file = c("GO_diff_yellow.rds"))

#KEGG富集分析，我们以总差异基因(diff_entrez)为例：
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa",#物种，Homo sapiens (human)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

#从KEGG富集分析结果中提取结果表格：
KEGG_result <- KEGG_diff@result
#View(KEGG_result)
#使用ggplot2进行更加自由的可视化呈现：
#先提取富集结果表前Top20：
KEGG_top20 <- KEGG_result[1:20,]

#指定绘图顺序（转换为因子）：
KEGG_top20$pathway <- factor(KEGG_top20$Description,levels = rev(KEGG_top20$Description))

#Top20富集数目条形图：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

p <- ggplot(data = KEGG_top20,
            aes(x = Count,
                y = pathway,
                fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdPu",direction = 1) +
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "KEGG enrichment barplot") + mytheme
p


library(topGO)
library(enrichplot)
####使用ggplot2进行可视化:
#取前top20，并简化命名：
MF <- GO_MF_result[1:20,]
CC <- GO_CC_result[1:20,]
BP <- GO_BP_result[1:20,]

#自定义主题
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))
#在MF的Description中存在过长字符串，我们将长度超过50的部分用...代替：
MF2 <- MF
MF2$Description <- str_trunc(MF$Description,width = 50,side = "right")
MF2$Description
MF2$term <- factor(MF2$Description,levels = rev(MF2$Description))
CC$term <- factor(CC$Description,levels = rev(CC$Description))
BP$term <- factor(BP$Description,levels = rev(BP$Description))

#GO富集柱形图：
GO_bar <- function(x){
  y <- get(x)
  ggplot(data = y,
         aes(x = Count,
             y = term,
             fill = -log10(pvalue))) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 50) ) + #label换行，部分term描述太长
    geom_bar(stat = "identity",width = 0.8) +
    labs(x = "Gene Number",
         y = "Description",
         title = paste0(x," of GO enrichment barplot")) +
    theme_bw() +
    mytheme
}

#MF:
p1 <- GO_bar("MF2")+scale_fill_distiller(palette = "Blues",direction = 1)
#CC:
p2 <- GO_bar("CC")+scale_fill_distiller(palette = "Reds",direction = 1)
#BP:
p3 <- GO_bar("BP")+scale_fill_distiller(palette = "Oranges",direction = 1)

##免疫热图###


###fig255免疫热图
setwd("~/TCGA/exosome_ccRCC/immnue_heatmap/")
library(utils)
library(GSVA)
library(ComplexHeatmap) # 用到里面的pheatmap画图然后拼图，需安装2.8以上版本的ComplexHeatmap
source("/home/data/vip39/huitu/FigureYa255TIME/pheatmap_translate.R") # 如果不想放弃较老版本的R及其对应的老ComplexHeatmap，可以只加载新版ComplexHeatmap里的这个函数，该脚本出自2.9.4版ComplexHeatmap
library(circlize)
library(viridis)
library(gplots)
library(data.table)
library(estimate)
source("/home/data/vip39/huitu/FigureYa255TIME/annTrackScale.R") # 加载函数，对数值进行标准化并对极端值做截断
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 加载自定义函数
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


# 读入用MOVICS获取的muscle-invasive bladder cancer (MIBC)亚型

annCol.tcga <- data.frame(CMOIC=paste0("CS",str_sub(sub$Cluster,2,2)))

rownames(annCol.tcga) <- sub$Sample
annCol.tcga <- annCol.tcga[rownames(annCol.tcga)%in%colnames(meth),]
table(colnames(meth)%in%rownames(annCol.tcga))

# 读取450K甲基化数据
# 甲基化数据只要用到5个探针的marker就可以，所以这里我的输入数据是简化的，方便传输
load("/home/data/vip39/TCGA/Immune_MOVICS_sub/KIRC_meth_metil.rds")
meth <- meth.metil
meth[1:4,1:4]
# 匹配亚型
#colnames(meth) <- substr(colnames(meth), start = 1,stop = 16)
#meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]
meth <- meth[,rownames(annCol.tcga)] 

MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
meth.metil <- meth[MeTIL.marker,]
MeTIL <- meth[MeTIL.marker,]
MeTIL <- t(scale(t(MeTIL)))

# 计算MeTIL得分
pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
MeTIL <- pca.MeTIL$rotation[,1]
annCol.tcga$MeTIL <- MeTIL



# 加载表达数据
tpm <- KIRC_mRNA_fpkm[,rownames(annCol.tcga)]
immune.signature <- read.table("/home/data/vip39/huitu/FigureYa255TIME/Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)

# 构建计算GSVA的列表
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}

# 免疫检查点相关基因
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9") 

# 免疫细胞的排序
immune.sig.ccr.order <- c("T.cells.CD8", 
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# 计算immune/stromal得分
indata <- log2(tpm + 1)
# 保存到文件
write.table(indata,file = "TCGA_log2TPM_hugo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

filterCommonGenes(input.f = "TCGA_log2TPM_hugo.txt" , output.f="TCGA_log2TPM_hugo_ESTIMATE.txt", id="GeneSymbol")

estimateScore("TCGA_log2TPM_hugo_ESTIMATE.txt","TCGA_log2TPM_hugo_estimate_score.txt", platform="affymetrix")

est.tcga <- read.table(file = "TCGA_log2TPM_hugo_estimate_score.txt",header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(indata)

# 对数值进行标准化并对极端值做截断
est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga)) 
rownames(est.tcga) <- colnames(tpm)

tcga.immune.gsva <- gsva(as.matrix(log2(tpm + 1)),
                         immune.sig.ccr,
                         method = "gsva")

# 设置颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
blue <- "#5bc0eb"
gold <- "#ECE700"
cyan <- "#00B3D0"

annCol.tcga$ImmuneScore <- as.numeric(est.tcga[rownames(annCol.tcga),"ImmuneScore"])
annCol.tcga$StromalScore <- as.numeric(est.tcga[rownames(annCol.tcga),"StromalScore"])
annCol.tcga <- annCol.tcga[order(annCol.tcga$CMOIC),]
annColors.tcga <- list() # 构建热图的图例颜色
annColors.tcga[["CMOIC"]] <- c("CS1" = clust.col[1],
                               "CS2" = clust.col[2]
                               #"CS3" = clust.col[3],
                               #"CS4" = clust.col[4]
)
annColors.tcga[["ImmuneScore"]] <- inferno(64)
annColors.tcga[["StromalScore"]] <- viridis(64)

## 热图1：免疫检查点基因表达
indata <- log2(tpm[intersect(rownames(tpm),imm.targets),] + 1)
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.tcga)],halfwidth = 2), # 表达谱数据标准化
                border_color = NA, # 热图单元格无边框
                annotation_col = annCol.tcga[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.tcga[c("CMOIC","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T, # 显示行名
                show_colnames = F, # 不显示列名
                cellheight = 12, # 热图高度固定
                cellwidth = 0.6, # 热图宽度固定
                name = "ICI", # 图例名字
                cluster_rows = F, # 行不聚类
                cluster_cols = F) # 列不聚类

#pdf("CheckPoint_heatmap.pdf",width = 10,height = 10)
hm1

dev.off()

## 热图2：肿瘤免疫微环境富集得分
hm2 <- pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(annCol.tcga)],halfwidth = 1), # 富集得分标准化并排序
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22), # 根据不同类别的免疫细胞分割
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "TIME",
                cluster_rows = F,
                cluster_cols = F)

#pdf("TIMEheatmap.pdf",width = 10,height = 10)
hm2

#dev.off()

## 热图3：MeTIL得分
hm3 <- pheatmap(standarize.fun(t(annCol.tcga[,"MeTIL",drop = F]),halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(c(cyan,"black","#F12AFE"),64),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "MeTIL",
                cluster_rows = F,
                cluster_cols = F)

#pdf("MeTILheatmap.pdf",width = 10,height = 10)
hm3
#dev.off()

# 合并热图并输出
pdf("TIME.pdf",width = 10,height = 10)
draw(hm1 %v% hm2 %v% hm3, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
invisible(dev.off())





#####突变大图248
##准备 突变2分类矩阵  tmb  突变签名  拷贝数 gistic2.0  亚型数据

library(ComplexHeatmap)
library(RColorBrewer)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor



load("/home/data/vip39/TCGA/Immune_MOVICS_sub/mutu248fig/KIRC_mut_and_tmb.rds")
## 突变二值矩阵
mut <- mut[,sub$Sample]
mut <- mut[,allid_need_fig]
## 突变频数
rownames(tmb) <- tmb$samID
tmb <- tmb["log10TMB"]
tmb <- tmb[allid_need_fig,]
tmb <- data.frame(log10TMB=tmb)
rownames(tmb) <- allid_need_fig
## 突变签名
#mutsig <- read.csv("mutationsignature.cosmic2013.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 或者换成下面这行，读取FigureYa110mutationSignature的输出文件mutsig.weightMatrix.txt
mutsig <- read.table("/home/data/vip39/TCGA/Immune_MOVICS_sub/mutu248fig/mutsig.weightMatrix.txt", row.names = 1, header = T) 
mutsig <- mutsig[allid_need_fig,]

## 拷贝数GISTIC2.0结果（利用SNP segment文件在GenePattern上获得）
cna.region <- read.table("/home/data/vip39/TCGA/Immune_MOVICS_sub/mutu248fig/all_lesions.conf_90.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dim(cna.region)
cna.region[1:4,1:4]
cna.gene <- read.table("/home/data/vip39/TCGA/Immune_MOVICS_sub/mutu248fig/all_thresholded.by_genes.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dim(cna.gene)
cna.gene[1:4,1:4]

## 亚型数据
subt <- annCol.tcga["CMOIC"]
subt <- data.frame(CMOIC=subt[allid_need_fig,])
rownames(subt) <- allid_need_fig
#allid_need_fig <- intersect(rownames(mutsig),colnames(cna.modified))


# 设置亚型颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
blue   <- "#5bc0eb"
red    <- "#f25f5c"

# 处理突变签名数据
mutsig <- mutsig[,c("Signature.1","Signature.2","Signature.5","Signature.13")] # 文章中用到3种类型的signature， SBS1 (age-related), SBS2 and SBS13 (APOBEC activity-related) and SBS5 (ERCC2 mutation-related)
mutsig$APOBEC <- mutsig$Signature.2 + mutsig$Signature.13 # APOBEC相关的签名由签名2和3叠加
mutsig$CMOIC <- subt[rownames(mutsig),"CMOIC"] # 添加亚型结果
mutsig <- mutsig[order(mutsig$CMOIC,-mutsig$APOBEC,decreasing = F),] # 确定整个热图的排序，按照亚型升序以及APOBEC降序排列


# 挑选要展示的基因
mutgene <- c("VHL",
             "PBRM1",
             "TTN",
             "SETD2",
             "BAP1",
             "MTOR",
             "MUC16",
             "KDM5C",
             "HMCN1",
             "DNAH9",
             "LRP2",
             "ATM",
             "ARID1A",
             "CSMD3",
             "DST",
             "KMT2C",
             "ERBB4",
             "SMARCA4",
             "USH2A")

# 制作oncoprint的输入数据
onco.input <- mut[mutgene,rownames(mutsig),]
onco.input[onco.input == 1] <- "Mutated" # 二值矩阵中1记为突变
onco.input[onco.input != "Mutated"] <- "" # 非“突变”给予空值
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Mutated = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A60000", col = "#A60000")) 
  }
)
col = c("Mutated" ="#A60000") # 突变颜色，注意这里只给了主图像的图例

my_ann <- subt[rownames(mutsig),,drop = F]
my_annotation = HeatmapAnnotation(df = my_ann, 
                                  col = list(CMOIC = c("CS1" = clust.col[1],
                                                       "CS2" = clust.col[2]
                                                       #"CS3" = clust.col[3],
                                                       #"CS4" = clust.col[4]
                                  )))

# 突变主区域的上部注释（突变负荷柱状图）
top_anno <- anno_barplot(as.numeric(tmb[rownames(mutsig),"log10TMB"]),
                         border = FALSE,
                         gp = gpar(fill = "#3379B4",border =NA,lty="blank"), 
                         height = unit(2.5, "cm"))

# 突变主区域的上部注释（突变签名柱状图）
tmp <- mutsig[,c("Signature.2","Signature.13","Signature.1","Signature.5")] # 只取和APOBEC有关的签名
tmp$Others <- 1 - rowSums(tmp) # 计算其他签名的比例
top_anno2 <- anno_barplot(as.matrix(tmp),
                          border = FALSE,
                          gp = gpar(fill = c(brewer.pal(6,"Paired")[c(2,1,6,5)],"grey90"), 
                                    border = NA, # 无边框
                                    lty = "blank"),
                          height = unit(2, "cm")) # 高度

tmp <- as.data.frame(t(mut[mutgene,rownames(mutsig),]))
mut.order <- names(sort(colSums(tmp),decreasing = T)) # 根据突变频数高低排序展示突变的顺序
tmp$CMOIC <- subt[rownames(tmp),"CMOIC"]
pct <- NULL # 计算各个基因突变的百分比
for (i in mut.order) {
  tmp1 <- tmp[,c(i,"CMOIC")]
  tmp1 <- as.data.frame.array(table(tmp1[,1],tmp1$CMOIC))[2,]/sum(tmp1[,1])
  pct <- rbind.data.frame(pct,tmp1)
}
rownames(pct) <- mut.order

# 添加右侧百分比堆积柱状图
right_anno <- anno_barplot(as.matrix(pct),
                           which = "row",
                           border = FALSE,
                           gp = gpar(fill = clust.col,border=NA,lty="blank"), 
                           bar_width = 0.6,
                           width = unit(1.8, "cm"),
                           height = unit(1, "cm"))

op1 <- oncoPrint(onco.input[mut.order,rownames(my_ann)], # 排序的突变矩阵
                 alter_fun = alter_fun,  # 主区域的函数，包括各单元格大小、背景颜色等等
                 col = col, # 突变颜色
                 bottom_annotation = NULL, # 无底部注释
                 top_annotation = c(HeatmapAnnotation(TMB = top_anno), # 顶部第一个注释：TMB
                                    my_annotation, # 顶部第二个注释：亚型
                                    HeatmapAnnotation(MutSig = top_anno2)), # 顶部第三个注释：突变签名
                 column_order = rownames(my_ann), # 样本的排序，根据突变签名的顺序
                 right_annotation = rowAnnotation(PCT = right_anno), # 右侧堆叠柱状图注释
                 show_pct = T, # 展示左侧的百分比
                 column_title = "", # 不显示主题
                 show_heatmap_legend = T, # 展示图例
                 column_split = my_ann$CMOIC, # 根据亚型切分热图
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))


op1

# 选择要展示的拷贝数区域
lesion.sig <- c(#"3q26.33-Amp",
  "3q26.2-Amp",
  #"5q35.1-Amp",
  #"7p22.2-Amp",
  #"12p11.23-Amp",
  "17q25.3-Amp",
  "1p36.11-Del",
  "2q37.3-Del",
  "3p25.3-Del",
  #"3q13.32-Del",
  #"6q25.2-Del",
  "6q26-Del",
  "9p21.3-Del",
  #"10q23.2-Del",
  "19p13.3-Del")

cna <- cna.region[1:(nrow(cna.region)/2),c(1,8,9:(ncol(cna.region)-1))] # 选取有效列
cna <- cna[!duplicated(cna$Descriptor),]
rownames(cna) <- paste0(gsub(" ","",cna$Descriptor),"-", substr(rownames(cna),1,3)) # 重命名行以确定扩增和缺失的位点
cna.modified <- cna[1:nrow(cna),3:ncol(cna)]
colnames(cna.modified) <- str_sub(colnames(cna.modified),1,15)
colnames(cna.modified)
#onco.input2 <- cna.modified[lesion.sig,rownames(mutsig)] # 选取要展示的拷贝数变异
table(rownames(mutsig)%in%colnames(cna.modified))

onco.input2 <- cna.modified[lesion.sig,rownames(subt)]
tmp1 <- onco.input2[1:2,] # 前6个为扩增
tmp1[tmp1 == 1] <- "Gain" # 数值大于0即为Gain
tmp1[tmp1 == 2] <- "Gain"
tmp1[tmp1 == 0] <- ""

tmp2 <- onco.input2[3:8,] # 后9个为缺失
tmp2[tmp2 == 1] <- "Loss"
tmp2[tmp2 == 2] <- "Loss"
tmp2[tmp2 == 0] <- ""
onco.input2 <- rbind.data.frame(tmp1,tmp2)

alter_fun2 = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = red, col = red)) 
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = blue, col = blue)) 
  }
)
col2 = c("Gain" = red,
         "Loss" = blue)
# 确定展示的顺序（看自己喜好，我这里是按照臂的顺序来的）
# lesion.order <- c("3q26.33-Amp","3q26.2-Amp","5q35.1-Amp","7p22.2-Amp","12p11.23-Amp","17q25.3-Amp",
#                   "1p36.11-Del","2q37.3-Del","3p25.3-Del","3q13.32-Del","6q25.2-Del","6q26-Del","9p21.3-Del","10q23.2-Del","19p13.3-Del")

lesion.order <- c(#"3q26.33-Amp",
  "3q26.2-Amp",
  #"5q35.1-Amp",
  #"7p22.2-Amp",
  #"12p11.23-Amp",
  "17q25.3-Amp",
  "1p36.11-Del",
  "2q37.3-Del",
  "3p25.3-Del",
  #"3q13.32-Del",
  #"6q25.2-Del",
  "6q26-Del",
  "9p21.3-Del",
  #"10q23.2-Del",
  "19p13.3-Del")
table(lesion.order%in%rownames(cna.modified))
lesion.order[lesion.order%in%rownames(cna.modified)]

tmp <- as.data.frame(t(cna.modified[lesion.order,rownames(mutsig),]))
tmp[tmp > 0] <- 1 # 所有大于1的均改为1以便计算变异频数
tmp$CMOIC <- as.character(subt[rownames(tmp),"CMOIC"])
pct <- NULL
for (i in lesion.order) {
  tmp1 <- tmp[,c(i,"CMOIC")]
  tmp1 <- as.data.frame.array(table(tmp1[,1],tmp1$CMOIC))[2,]/sum(tmp1[,1])
  pct <- rbind.data.frame(pct,tmp1)
}
rownames(pct) <- lesion.order

# 右侧堆叠百分比柱状图
right_anno2 <- anno_barplot(as.matrix(pct),
                            which = "row",
                            border = FALSE,
                            gp = gpar(fill = clust.col,
                                      border = NA,
                                      lty = "blank"), 
                            bar_width = 0.6,
                            width = unit(1.8, "cm"),
                            height = unit(1, "cm"))

# 同样的方式绘制热图
op2 <- oncoPrint(onco.input2[lesion.order,rownames(my_ann)], 
                 alter_fun = alter_fun2, 
                 col = col2, 
                 bottom_annotation = NULL, 
                 top_annotation = NULL,
                 column_order = rownames(my_ann),
                 right_annotation = rowAnnotation(PCT = right_anno2),
                 row_order = lesion.order, 
                 show_pct = T,
                 column_title = "", 
                 show_heatmap_legend = T, 
                 column_split = my_ann$CMOIC,
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))


op2
colnames(cna.gene) <- str_sub(colnames(cna.gene),1,15)

cna <- cna.gene
cna <- cna[c("VHL","PBRM1","TTN","SETD2"),rownames(mutsig)] # 文章筛选了4个基因
table(c("VHL","PBRM1","TTN","SETD2")%in%rownames(cna.gene))
onco.input3 <- cna
# 由于上面的分析中缺失的部分不存在High balanced loss，所以直接让2也属于Gain而没有添加High balanced gain
# 这里小伙伴如果自己的数据满足要求，上面的拷贝数也可以分为4类
onco.input3[onco.input3 == 1] <- "Gain"
onco.input3[onco.input3 == 2] <- "High_balanced_gain"
onco.input3[onco.input3 == 0] <- ""
onco.input3[onco.input3 == -1] <- "Loss"
onco.input3[onco.input3 == -2] <- "High_balanced_loss"

alter_fun3 = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[5], col = brewer.pal(6,"Paired")[5])) 
  },
  High_balanced_gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[6], col = brewer.pal(6,"Paired")[6]))
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[1], col = brewer.pal(6,"Paired")[1])) 
  },
  High_balanced_loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = brewer.pal(6,"Paired")[2], col = brewer.pal(6,"Paired")[2]))
  }
)
col3 = c("Gain" = brewer.pal(6,"Paired")[5],
         "High_balanced_gain" =brewer.pal(6,"Paired")[6],
         "Loss" = brewer.pal(6,"Paired")[1],
         "High_balanced_loss" =brewer.pal(6,"Paired")[2])

op3 <- oncoPrint(onco.input3[,rownames(my_ann)], 
                 alter_fun = alter_fun3,  
                 col = col3, 
                 bottom_annotation = NULL, 
                 top_annotation = NULL,
                 column_order = rownames(my_ann), 
                 right_annotation = NULL,
                 show_pct = T, 
                 column_title = "", 
                 show_heatmap_legend=T, 
                 column_split = my_ann$CMOIC,
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))

op3

# 构建额外图例
lgd.mutsig = Legend(labels = c("SBS2","SBS13","SBS1","SBS5","Others"), 
                    title = "MutSig", 
                    legend_gp = gpar(fill = c(brewer.pal(6,"Paired")[c(2,1,6,5)],"grey90")))

lgd.cna.region = Legend(labels = c("Gain","Loss"), 
                        title = "CNA (arm-level)", 
                        legend_gp = gpar(fill = c(red,blue)))

lgd.cna.gene = Legend(labels = c("Gain","High_balanced_gain","Loss","High_balanced_loss"), 
                      title = "CNA (gene-level)", 
                      legend_gp = gpar(fill = brewer.pal(6,"Paired")[c(5,6,1,2)]))              

lgd_list <- list(lgd.mutsig, lgd.cna.region, lgd.cna.gene)

# 合并热图
pdf("mutational landscape in TCGA.pdf", width = 10,height = 10)
draw(op1 %v% op2 %v% op3, # 垂直叠加热图
     annotation_legend_list = lgd_list) # 添加自定义的图例
invisible(dev.off())




