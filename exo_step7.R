

###计算分数 risk score 238突变 228蝴蝶  227药敏
.libPaths(c("/home/data/refdir/Rlib","/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1"))

.libPaths(c("/home/data/refdir/Rlib",
            "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1",
            "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/"))

.libPaths(c("/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/",
            "/home/data/refdir/Rlib",
            "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1"))

library(data.table)
library(GSVA)
library(ggplot2)
library(ggcorrplot)
#library(ggcor)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 读取gmt文件为列表形式，满足pathifier的输入要求
gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

# 加载TCGA-BLCA表达谱FPKM
expr <- read.table("tcga_blca_fpkm.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tumsam <- colnames(expr)[substr(colnames(expr),11,13) == "01A"] # 取出肿瘤样本
# 取出SIGLEC15表达量
#SIGLEC15 表达量改为 risk score

siglec15 <- log2(expr["SIGLEC15",tumsam] + 1)
siglec15 <- sub$gpr_score

#write.csv(siglec15, "easy_input_gene.csv")



# 加载immunotherapy-predicted pathways
immPath <- read.table("~/huitu/FigureYa228linkCor/immunotherapy predicted pathways.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
immPath.list <- list()
for (i in rownames(immPath)) {
  tmp <- immPath[i,"Genes"]
  tmp <- toupper(unlist(strsplit(tmp,", ",fixed = T)))
  tmp <- gsub(" ","",tmp)
  immPath.list[[i]] <- tmp
}

# 加载膀胱癌相关签名
##这里改为step1-15的免疫步骤
blcaPath <- read.table("~/huitu/FigureYa228linkCor/bladder cancer signature.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
blcaPath.list <- list()
for (i in rownames(blcaPath)) {
  tmp <- blcaPath[i,"Genes"]
  tmp <- toupper(unlist(strsplit(tmp,", ",fixed = T)))
  tmp <- gsub(" ","",tmp)
  blcaPath.list[[i]] <- tmp
}




# 计算immunotherapy-predicted pathways的富集得分
KIRC_tpm
range(KIRC_tpm)
dim(KIRC_tpm)


immPath.score <- gsva(expr = as.matrix(KIRC_tpm[,rownames(sub)]),
                      immPath.list, 
                      method = "ssgsea")
immPath.score[1:4,1:4]
dim(immPath.score)
colnames(immPath.score) <- str_sub(colnames(immPath.score),1,15)
# 保存到文件，便于套用
write.table(immPath.score, "easy_input_immPath.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


# 计算膀胱癌相关签名的富集得分
blcaPath.score <- gsva(expr = as.matrix(log2(expr[,tumsam] + 1)),
                       blcaPath.list,
                       method = "ssgsea")

stepscore[1:4,1:4]
dim(sub)
sub$sample <- str_sub(rownames(sub),1,15)

blcaPath.score <- stepscore[,sub$sample]
blcaPath.score[1:4,1:4]
# 保存到文件，便于套用
write.table(blcaPath.score, "easy_input_blcaPath.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


# 从gmt文件提取感兴趣的通路
# 例如提取带有`METABOLISM`字样的通路
gset <- gmt2list("c5.all.v7.4.symbols.gmt") 
gset.list <- gset[names(gset) %like% "METABOLISM"]

# 计算通路富集得分
gset.score <- gsva(expr = as.matrix(log2(expr[,tumsam] + 1)),
                   gset.list,
                   method = "ssgsea")
# 保存到文件
write.table(gset.score, "easy_input_gset.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


##开始绘图
# 加载单个基因的表达量，也就是“1对多”的“1”
siglec15 <- read.csv("easy_input_gene.csv", row.names = 1, check.names = F)
colnames(siglec15)[1:4]


# 加载第一组富集得分
immPath.score <- read.table("easy_input_immPath.txt", check.names = F)
# 跟目标基因SIGLEC15的表达量合并
immPath.score <- rbind.data.frame(immPath.score,
                                  siglec15)

# 加载第二组富集得分
blcaPath.score <- read.table("easy_input_blcaPath.txt", check.names = F)
# 跟目标基因SIGLEC15的表达量合并
blcaPath.score <- rbind.data.frame(blcaPath.score,
                                   siglec15)
# 循环计算相关性并绘制左下角
immCorSiglec15 <- NULL
for (i in rownames(immPath.score)) {
  cr <- cor.test(as.numeric(immPath.score[i,]),
                 as.numeric(siglec15),
                 method = "pearson")
  immCorSiglec15 <- rbind.data.frame(immCorSiglec15,
                                     data.frame(gene = "Siglec15",
                                                path = i,
                                                r = cr$estimate,
                                                p = cr$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
immCorSiglec15$sign <- ifelse(immCorSiglec15$r > 0,"pos","neg")
immCorSiglec15$absR <- abs(immCorSiglec15$r)
immCorSiglec15$rSeg <- as.character(cut(immCorSiglec15$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
immCorSiglec15$pSeg <- as.character(cut(immCorSiglec15$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
immCorSiglec15[nrow(immCorSiglec15),"pSeg"] <- "Not Applicable"

immCorSiglec15$rSeg <- factor(immCorSiglec15$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCorSiglec15$pSeg <- factor(immCorSiglec15$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCorSiglec15$sign <- factor(immCorSiglec15$sign, levels = c("pos","neg"))

p1 <- quickcor(t(immPath.score), 
               type = "lower",
               show.diag = TRUE) + 
  geom_colour() +
  add_link(df = immCorSiglec15, 
           mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
           spec.key = "gene",
           env.key = "path",
           diag.label = FALSE) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#19A078","#DA6003","#7570B4","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "navy",mid = "white",high = "orange",midpoint=0) +
  remove_axis("x")
p1
##修改热图颜色

ggsave(filename = "ggcor plot in bottom left.pdf", width = 10,height = 8)

# 循环计算相关性并绘制右上角
blcaCorSiglec15 <- NULL
for (i in rownames(blcaPath.score)) {
  cr <- cor.test(as.numeric(blcaPath.score[i,]),
                 as.numeric(siglec15),
                 method = "pearson")
  blcaCorSiglec15 <- rbind.data.frame(blcaCorSiglec15,
                                      data.frame(gene = "Siglec15",
                                                 path = i,
                                                 r = cr$estimate,
                                                 p = cr$p.value,
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
}
blcaCorSiglec15$sign <- ifelse(blcaCorSiglec15$r > 0,"pos","neg")
blcaCorSiglec15$absR <- abs(blcaCorSiglec15$r)
blcaCorSiglec15$rSeg <- as.character(cut(blcaCorSiglec15$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
blcaCorSiglec15$pSeg <- as.character(cut(blcaCorSiglec15$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
blcaCorSiglec15[nrow(blcaCorSiglec15),"pSeg"] <- "Not Applicable"

blcaCorSiglec15$rSeg <- factor(blcaCorSiglec15$rSeg, levels = c("0.25","0.50","0.75","1.00"))
blcaCorSiglec15$pSeg <- factor(blcaCorSiglec15$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
blcaCorSiglec15$sign <- factor(blcaCorSiglec15$sign, levels = c("pos","neg"))

p2 <- quickcor(t(blcaPath.score), 
               type = "upper",
               show.diag = TRUE) + 
  geom_colour() +
  add_link(df = blcaCorSiglec15, 
           mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
           spec.key = "gene",
           env.key = "path",
           diag.label = FALSE) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#19A078","#DA6003","#7570B4","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "navy",mid = "white",high = "orange",midpoint=0) +
  remove_axis("x")
p2

ggsave(filename = "ggcor plot in top right.pdf", width = 10,height = 8)



### riskscore 趋势图  


##risk score nomogram  

#time roc 竞争曲线

##不同队列的 ROC 1 3 5年 优于其他队列的表现
##tcga 1,3,5 roc


##japan
load("/home/data/vip39/database/KIRC_immune_therapy/japan_cohort/japan_RCCall.rds")




##cancer cell 
load("/home/data/vip39/database/KIRC_immune_therapy/cancercellpaper/cancercellRCC.rds")

##gse29609
load("/home/data/vip39/database/KIRC_immune_therapy/GSE29609/GSE29609_RCCall.rds")




####免疫相关性分析 riskscore#####
#fig56
library(xCell)
library(circlize)


##TIME图#####
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
library(ComplexHeatmap)
library(grid)
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
table(sub$cluster)

annCol.tcga <- data.frame(CMOIC=paste0("CS",str_sub(sub$cluster,2,2)))

rownames(annCol.tcga) <- str_sub(rownames(sub),1,15)
inter_id <- intersect(rownames(annCol.tcga),colnames(meth))

# 读取450K甲基化数据
# 甲基化数据只要用到5个探针的marker就可以，所以这里我的输入数据是简化的，方便传输
load("/home/data/vip39/TCGA/Immune_MOVICS_sub/KIRC_meth_metil.rds")

meth <- meth.metil
meth <- meth[,inter_id]

annCol.tcga <- annCol.tcga[inter_id,]
annCol.tcga <- data.frame(CMOIC=annCol.tcga)
rownames(annCol.tcga) <- inter_id
# 匹配亚型
#colnames(meth) <- substr(colnames(meth), start = 1,stop = 16)
#meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]
#meth <- meth[,rownames(annCol.tcga)] 

MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
meth.metil <- meth[MeTIL.marker,]
MeTIL <- meth[MeTIL.marker,]
MeTIL <- t(scale(t(MeTIL)))

# 计算MeTIL得分
pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
MeTIL <- pca.MeTIL$rotation[,1]
annCol.tcga$MeTIL <- MeTIL



# 加载表达数据
KIRC_cancer_tpm <- KIRC_tpm[,str_sub(colnames(KIRC_tpm),14,15)=="01"]
colnames(KIRC_cancer_tpm) <- str_sub(colnames(KIRC_cancer_tpm),1,15)
tpm <- KIRC_cancer_tpm[,rownames(annCol.tcga)]
immune.signature <- read.table("/home/data/vip39/huitu/FigureYa255TIME/Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)

# 构建计算GSVA的列表
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}

# 免疫检查点相关基因
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4") 

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
range(tpm)
indata <- tpm
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
library(GSVA)
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
library(grid)
library(ComplexHeatmap)
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






#####热图补充 low risk and normal tissue##### 
#加上正常组织三组绘图
#dir.create("../complexheatmap")
#setwd("/home/data/vip39/TCGA/met_stage1/complexheatmap/")
library(ComplexHeatmap)
library(stringr)
library(pheatmap)
library(gplots)
library(grid)
library(circlize)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
###准备两个数据 mygene_data 和 Subtype 分组信息
range(KIRC_tpm)
KIRC_tpm[1:4,1:4]
#load("~/vip39/TCGA/exosome_ccRCC/exosome_sub.rds")
head(sub)
sub$Sample <- paste0(sub$Sample,"A")
normal_id <- colnames(KIRC_tpm)[str_sub(colnames(KIRC_tpm),14,15)=="11"]
head(sub)

c(normal_id,sub$Sample)

Subtype <- data.frame(Subtype=c(rep("N",length(normal_id)),as.character(sub$Cluster)),
                      id=c(normal_id,sub$Sample))
rownames(Subtype) <- Subtype$id
head(Subtype)
table(Subtype$Subtype)
class(Subtype)

#Subtype <- Subtype[,-2]
Subtype <- data.frame(Subtype=Subtype$Subtype,
                      row.names = rownames(Subtype))
RNA_gene <- Coxoutput$gene
dim(KIRC_tpm)
KIRC_tpm[1:4,1:4]

table(rownames(Subtype)%in%colnames(KIRC_tpm))
interid <- intersect(rownames(Subtype),colnames(KIRC_tpm))
Subtype <- Subtype[interid,]
Subtype <- data.frame(Subtype=Subtype,
                      row.names = interid)
rownames(Subtype) <- interid
table(rownames(Subtype)%in%colnames(KIRC_tpm))

mygene_data <- KIRC_tpm[RNA_gene,rownames(Subtype)]
#Subtype <- Subtype[com_sam,,drop = F]
head(Subtype)
table(Subtype$Subtype)
## 用前面的自定义函数计算组间统计差异
#setwd("./supplementaryresult/")

dim(mygene_data)
dim(Subtype)

comprTab <- cross_subtype_compr(expr = mygene_data, # 或log2(mygene_data + 1)，如果用参数检验，请注意对数转化；若非参均可
                                subt = Subtype,
                                #two_sam_compr_method = "wilcox", # 两组"t.test", "wilcox"
                                multi_sam_compr_method = "kruskal", # 多组"anova", "kruskal"
                                res.path = ".")


# 用全部基因来画
n.show_top_gene <- nrow(mygene_data)

# 按分组排序
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]

# 数据标准化和边界设置
indata <- log2(indata+1)
plotdata <- t(scale(t(indata)))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

# 调整行名
blank <- "    " # 行名和p值之间的间隔
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
p.label <- formatC(p.value, # 将p值变成保留两位小数的科学计数法
                   format = "e",
                   digits = 2)

add.label <- str_pad(paste0(rownames(plotdata),sig.label), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")

annCol <- subt.order # 获得排序后的亚型注释信息，这里只有一个变量需要注释
colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], # 注释列名补上"P-value"，宽度和刚才一致
                                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                                     side = "right"),
                             "P-value",
                             sep = blank)

annColors <- list(c( "C1"="#E31A1CFF", "C2"="#1F78B4FF","N"="#33A02CFF")) # 如果有多个变量要注释颜色请补充c()
names(annColors) <- colnames(annCol)[1] # 如果有多个变量要注释颜色请补充每张list的name

# 绘制热图
library(circlize)
library(stringr)
library(pheatmap)
library(gplots)
library(grid)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
table(sub$Cluster)
dim(plotdata)
pheatmap(cellheight = 9, cellwidth = 1,
         mat = plotdata, # 输入数据
         scale = "none", # 不标准化因为数据已经被标准化
         annotation_col = annCol, # 列注释信息
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_colnames = F, # 不显示列名
         show_rownames = T, # 显示基因名
         annotation_legend = F, # 不显示图例
         #gaps_col = c(323,599),
         color = colorRampPalette(c('#4D4398',"white",'#F18D00'))(100),
         labels_row = paste(add.label, p.label, sep=blank),
)

