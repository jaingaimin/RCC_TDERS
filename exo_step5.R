
##其他维度的分析 基于TCGA-KIRC 


##riskscore 相关性计算 figureya228

##



##药敏 使用四套方法 交叉结果验证 fig212  fig105 fig282  fig213

## 突变 增加figure238相关性 计算

#DNA甲基化数据

##模型多重验证  roc c-index 矫正曲线 DCA曲线

##fig161 分组和多个临床特征的相关性


##药敏增加部分
library(PharmacoGx)
library(parallel)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 以下生成signature基因的代码只适用于TCGA来源的数据
load("/home/data/refdir/database/tcga_counts_fpkm_tpm/TCGA-KIRC_tpm_gene_symbol.Rdata")
tpms <- as.data.frame(tpms)
colnames(tpms)
tpms <- tpms[,str_sub(colnames(tpms),14,15)=="01"]
colnames(tpms) <- str_sub(colnames(tpms),1,15) 
data <- tpms[,sub$Sample]

# 得到癌与癌旁比较的logfc
tumor <- data[, sub$Sample[sub$group=="High"]]
normal <- data[, sub$Sample[sub$group=="Low"]]

dis_sig <- data.frame(id=rownames(data),
                      fc=log2(rowMeans(tumor)/rowMeans(normal)))
dis_sig <- dis_sig[dis_sig$fc != "Inf" & dis_sig$fc != "-Inf",]
dis_sig <- na.omit(dis_sig)

# 提取变化倍数最大的基因
# 例文（Yang et al）提到疾病分子特征数量选择100时可以获得较好的预测性能。但这篇研究是基于LINCS数据，cmap数据的维度与LINCS差别明显。这里我们建议疾病分子特征数量可以稍微多一点，以下演示使用top300基因进行XSum分析
dis_sig <- rbind(top_n(dis_sig, 150, fc), top_n(dis_sig, -150, fc))

# 将logfc转成二分类变量（logfc>0的基因用1表示，logfc小于0的基因用-1表示）
# 使用XSum时不需要考虑差异基因的差异倍数，这步分析是为了让大家更好的理解，并不是并要的

dis_sig$fc[dis_sig$fc>0] <- 1; dis_sig$fc[dis_sig$fc<0] <- -1
rownames(dis_sig) <- NULL

# 保存结果
write.table(dis_sig, "dis_sig.csv", sep=",", quote=F, row.names=F, col.names=T)

# dis_sig只要整理成以下格式都可以在后面的分析中使用
head(dis_sig)




##GSEA美化

