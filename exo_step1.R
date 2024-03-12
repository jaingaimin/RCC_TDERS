
.libPaths()
.libPaths(c(
  #"~/vip39/R/x86_64-pc-linux-gnu-library/4.2/",
  "/home/data/refdir/Rlib/",
  "~/vip39/R/x86_64-pc-linux-gnu-library/4.0/",
  "~/vip39/R/x86_64-pc-linux-gnu-library/4.1/",
  "/home/data/t030432/R/x86_64-pc-linux-gnu-library/4.2"))


####外泌体分析 up signature
setwd("~/vip39/TCGA/exosome_ccRCC/")
library(data.table)
exo_deg <- fread("~/TCGA/exosome_ccRCC/emRNA_deg.csv",header = T,data.table = F)
colnames(exo_deg)
exo_deg_up <- exo_deg[exo_deg$log2FoldChange>0,]
exo_deg_up <- exo_deg_up$`gene name`


##开始构建模型 risk score
##这里可以模仿
#已经缩小了范围 外泌体来源的DEG
#后续就是构建模型  TCGA 和JAPAN 队列
library(data.table)
library(limma) # 芯片差异表达
library(impute,lib.loc = "/home/data/refdir/Rlib/") # 芯片缺失值多重填补
library(dplyr)
library(pheatmap)
library(gplots)
library(pROC)
library(ggplot2)
library(stringr)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 自定义函数显示进度
display.progress = function ( index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
} 

# 自定义函数标准化表达谱
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

##加载TCGA
load("/home/data/refdir/database/tcga_counts_fpkm_tpm/TCGA-KIRC_tpm_gene_symbol.Rdata")
KIRC_tpm <- as.data.frame(tpms)
KIRC_tpm[1:4,1:4]
KIRC_tpm <- KIRC_tpm[,str_sub(colnames(KIRC_tpm),14,15)!=11]
colnames(KIRC_tpm) <- str_sub(colnames(KIRC_tpm),1,12)
clin <- fread("/home/data/refdir/database/TCGA临床数据/肾透明细胞癌/肾透明细胞癌.txt",header = T,data.table = F)
rownames(clin) <- clin$id
clin[1:4,1:4]

comsam <- intersect(colnames(KIRC_tpm), rownames(clin))
kirc.expr <- log2(KIRC_tpm[,comsam] + 1)
kirc.sinfo <- clin[comsam,,drop = F]
kirc.sinfo$fustat <- ifelse(kirc.sinfo$status=="Alive",0,1)




##加载JAPAN
load("/home/data/vip39/database/KIRC_validation_data/japan_input.rds")
japan_clin <- fread("/home/data/vip39/database/KIRC_validation_data/japan_clin_clean.csv",header=T,data.table=F)
rownames(japan_clin) <- japan_clin$sample_ID
colnames(japan_clin)
sub5 <- japan_clin[,c("sample_ID","outcome","month")]
colnames(sub5) <- c("ID","OS","OS.time")
str(sub5)
sub5$OS <- ifelse(sub5$OS=="alive",0,1)
rownames(sub5) <- sub5$ID
sub5 <- sub5[colnames(japan_mRNA_fpkm),]
japan_mRNA_fpkm <- japan_mRNA_fpkm%>%
  column_to_rownames("gene_name")

japan_mRNA_fpkm[1:4,1:4]
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

japan_tpms<-as.data.frame(apply(japan_mRNA_fpkm,2,fpkmToTpm))%>%
  rownames_to_column(var='gene_id')%>%
  mutate(gene_id1=substr(.$gene_id,1,15))%>%
  dplyr::select(gene_id1,everything(),-gene_id)%>%
  column_to_rownames("gene_id1")
japan_tpms[1:4,1:4]

comsam <- intersect(colnames(japan_tpms), rownames(sub5))
japan.expr <- log2(japan_tpms[,comsam] + 1)
japan.sinfo <- sub5[comsam,,drop = F]
japan.sinfo$fustat <- japan.sinfo$OS



###在两个测试集里寻找预后判断力都好的基因集合
# 计算auc
comgene <- intersect(intersect(rownames(japan.expr),
                               rownames(kirc.expr)),
                     exo_deg_up)
auc1 <- auc2 <- c()
for (i in comgene) {
  # TCGA-KIRC
  tmp <- data.frame(expr = as.numeric(kirc.expr[i,]),
                    fustat = kirc.sinfo$fustat,
                    stringsAsFactors = F)
  auc1 <- c(auc1, roc(predictor = tmp$expr, response = tmp$fustat)$auc)
  
  # GSE10186
  tmp <- data.frame(expr = as.numeric(japan.expr[i,]),
                    fustat = japan.sinfo$fustat,
                    stringsAsFactors = F)
  auc2 <- c(auc2,roc(predictor = tmp$expr,response = tmp$fustat)$auc)
}
names(auc1) <- names(auc2) <- comgene

# 原文在Human Protemome Project上去找哪些定义的基因是有人类蛋白表达的，这里我就随便取了几个基因，只是为了绘图方便而已。
# 找到满足条件的基因
#需要确定蛋白
##可以选择分泌蛋白 或者 叶定伟数据中差异蛋白  

comgene_df <- data.frame(gene=comgene)
write.csv(comgene_df,file = "~/TCGA/exosome_ccRCC/comgene.csv",quote = F)
protein_evi <- fread('~/TCGA/exosome_ccRCC/protein_evidence.tsv',header = T,data.table = F)

FAKE.protein.positive <- intersect(protein_evi$Gene,comgene) #随机取了一些基因作为蛋白表达基因
FAKE.protein.negtive <- setdiff(comgene, FAKE.protein.positive)

auc1.cutoff <- 0.55
auc2.cutoff <- 0.55
finalgene <- intersect(comgene[auc1 > auc1.cutoff & auc2 > auc2.cutoff], FAKE.protein.positive) # 既满足阈值又是蛋白表达的基因
avgauc <- c(auc1[finalgene] + auc2[finalgene])/2

pickgene <- "VHL" # 因为和原文结果不同所以随便选的

##开始绘图
xrange <- pretty(range(auc1))
yrange <- pretty(range(auc2))

par(mfrow = c(1,2)) # 把画布分成左右两个
par(bty="l", mgp = c(2.4,.33,0), mar=c(4.1,4.6,2.1,2.1)+.1, las=1, tcl=-.25)

# 绘制“假的”有蛋白表达的散点
plot(x = auc1[FAKE.protein.positive],
     y = auc2[FAKE.protein.positive],
     xlab = "AUC for disease event\nin TCGA-KIRC",
     ylab = "AUC for disease event\nin JAPAN-KIRC",
     pch = 19,
     col = "#5D7ABE",
     cex = 1.2,
     xlim = c(xrange[1],xrange[length(xrange)]),
     ylim = c(yrange[1],yrange[length(yrange)]),
     xaxt = "n",
     yaxt = "n")
# 添加“假的”无蛋白表达的散点
points(x = auc1[FAKE.protein.negtive],
       y = auc2[FAKE.protein.negtive],
       pch = 19,
       cex = 1.2,
       col = ggplot2::alpha("grey70", 0.8))
# 添加“假的”有蛋白表达且AUC满足阈值且感兴趣的目标基因
points(x = auc1[pickgene],
       y = auc2[pickgene],
       pch = 19,
       cex = 1.2,
       col = ggplot2::alpha("red", 0.8))
axis(side = 1,at = xrange)
axis(side = 2,at = yrange)
# 添加阈值区域框
rect(xleft = auc1.cutoff,xright = max(auc1) + 0.01,
     ybottom = auc2.cutoff,ytop = max(auc2) + 0.01,
     lwd = 1.5,
     lty = 2)

# 右侧绘制AUC均值
par(bty="l", mgp = c(1.9,.33,0), mar=c(4.1,5.1,2.1,2.1)+.1, las=1, tcl=-.25)
a <- barplot(sort(avgauc,decreasing = F),
             horiz = T,
             col = ifelse(names(sort(avgauc,decreasing = F)) == pickgene,"red","#5D7ABE"),
             xlab = "Mean AUC")
axis(side = 2,at = a,labels = F)
dev.copy2pdf(file = "combined pairwise auc.pdf", width = 9, height = 5)


###
save(finalgene,file = '~/TCGA/exosome_ccRCC/finalgene50_input.rds')

##开始构建模型
#lasso联合bostrap






