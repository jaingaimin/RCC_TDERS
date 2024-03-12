
.libPaths(c(
  #"~/vip39/R/x86_64-pc-linux-gnu-library/4.2/",
  "~/vip39/R/x86_64-pc-linux-gnu-library/4.1/",
  "/home/data/refdir/Rlib/",
  "~/vip39/R/x86_64-pc-linux-gnu-library/4.0/",
  
  "/home/data/t030432/R/x86_64-pc-linux-gnu-library/4.2"))

#### ifelse(sub_exo$group=="High","C1","C2")

setwd("~/TCGA/exosome_ccRCC/subcompare")
load("~/vip39/TCGA/exosome_ccRCC/tcgasub.rds")
table(data$group)
# Low High 
# 263  263 
dim(sub)
#  [1] 511  22
#这一步最后只剩下511个样本


###直接load   进行常规分析
sub <- data
sub$cluster <- ifelse(sub$group=="High","C2","C1")
table(sub$cluster)

##进行movics分析
# 显示进程
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
} 

# 计算组间统计差异
cross_subtype_compr <- function(expr = NULL,
                                subt = NULL,
                                subt.label = "Subtype",
                                two_sam_compr_method = "wilcox",
                                multi_sam_compr_method = "kruskal",
                                res.path = NULL) {
  
  if (!is.element(two_sam_compr_method, c("t.test", "wilcox"))) {stop("Two samples comparison should be t.test or wilcox!\n") }
  if (!is.element(multi_sam_compr_method, c("anova", "kruskal"))) {stop("multiple samples comparison should be kruskal or anova!\n") }
  
  subt.name <- unique(subt[,subt.label])
  n.subt <- length(subt.name)
  if(n.subt < 2) {stop("The number of subtype should be greater than 2!\n")}
  
  comprTab <- NULL
  
  # 两个亚型且为非参数检验
  if(n.subt == 2 & two_sam_compr_method == "wilcox") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp1 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[1]),,drop = F])])
      tmp2 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[2]),,drop = F])])
      wt <- wilcox.test(tmp1,tmp2)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = wt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 两个亚型且为参数检验
  if(n.subt == 2 & two_sam_compr_method == "t.test") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp1 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[1]),,drop = F])])
      tmp2 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[2]),,drop = F])])
      tt <- t.test(tmp1,tmp2)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = tt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 多个亚型且为非参数检验
  if(n.subt > 2 & multi_sam_compr_method == "kruskal") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp.list <- list()
      for (n in 1:n.subt) {
        tmp.list[[n]] <- data.frame(value = as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[n]),,drop = F])]),
                                    subt = subt.name[n],
                                    stringsAsFactors = F)
      }
      tmp <- do.call(rbind,tmp.list)
      kt <- kruskal.test(value ~ subt,data = tmp)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = kt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 多个亚型且为参数检验
  if(n.subt > 2 & multi_sam_compr_method == "anova") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp.list <- list()
      for (n in 1:n.subt) {
        tmp.list[[n]] <- data.frame(value = as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[n]),,drop = F])]),
                                    subt = subt.name[n],
                                    stringsAsFactors = F)
      }
      tmp <- do.call(rbind,tmp.list)
      at <- summary(aov(value ~ subt,data = tmp))
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = at[[1]][1,5],
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 调整p值
  comprTab$adjusted.p.value = p.adjust(comprTab$nominal.p.value,method = "BH")
  # 按p值排序
  #comprTab <- comprTab[order(comprTab$adjusted.p.value, decreasing = F),] 
  
  write.table(comprTab,file.path(res.path,"comprTab.txt"),sep = "\t",row.names = F,quote = F)
  return(comprTab)
}

##和多个预后 免疫指标关联
#CD8, PD-1 and PD-L1
#相关性分析 绘图 高级版

##外扩到其他队列 免疫治疗队列尤其

##肾癌其他免疫治疗队列  N NR


setwd("~/vip39/TCGA/ex")
#####热图补充 low risk and normal tissue##### 
#加上正常组织三组绘图
#dir.create("../complexheatmap")
setwd("/home/data/vip39/TCGA/met_stage1/complexheatmap/")
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
load("~/vip39/TCGA/exosome_ccRCC/exosome_sub.rds")
head(sub)
sub$Sample <- paste0(sub$Sample,"A")
normal_id <- colnames(KIRC_tpm)[str_sub(colnames(KIRC_tpm),14,15)=="11"]
head(sub)

c(normal_id,sub$Sample)

Subtype <- data.frame(Subtype=c(rep("N",length(normal_id)),sub$Cluster),
                      id=c(normal_id,sub$Sample))
rownames(Subtype) <- Subtype$id
Subtype <- Subtype[-2]
RNA_gene <- gene
mygene_data <- KIRC_tpm[RNA_gene,rownames(Subtype)]
#Subtype <- Subtype[com_sam,,drop = F]
head(Subtype)
table(Subtype$Subtype)
## 用前面的自定义函数计算组间统计差异
setwd("./supplementaryresult/")

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
         gaps_col = c(258,511),
         color = colorRampPalette(c('#4D4398',"white",'#F18D00'))(100),
         labels_row = paste(add.label, p.label, sep=blank),
)


######cox分析#####
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)

dim(sub)
Coxoutput <- data.frame(OS=sub$fustat,
                        OS.time=sub$futime)
Coxoutput <- data.frame(OS=sub$PFI,
                        OS.time=sub$PFI.time)
rownames(Coxoutput) <- sub$Sample

Coxoutput2 <- as.data.frame(t(KIRC_tpm[unique(RNA_gene),sub$Sample]))
range(Coxoutput2)
#Coxoutput2 <- log2(Coxoutput2+1)
Coxoutput <- cbind(Coxoutput,Coxoutput2)

realdata <- Coxoutput
realdata[1:3,1:6]
#setwd("~/vip39/TCGA/KIRC_MLcelldeath/result/")
dir.create("COX")
setwd("./COX/")
str(realdata)

Coxoutput=data.frame()
for(i in colnames(realdata[,3:ncol(realdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ realdata[,i], data = realdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
# Coxoutput <- arrange(Coxoutput,pvalue)  %>% #按照p值排序
#   filter(pvalue < 0.01) 
Coxoutput <- arrange(Coxoutput,pvalue)
#Coxoutput <- Coxoutput[Coxoutput$pvalue<0.01,]

dim(Coxoutput)
#[1] 12  6
#保存到文件
write.csv(Coxoutput,'cox_output.csv', row.names = F)
Coxoutput <- read.csv("cox_output.csv")
head(Coxoutput)

#plotCoxoutput <- filter(Coxoutput,HR <=0.92 | HR>= 1.15)  #选择合适的top genes
plotCoxoutput <- Coxoutput
#采用那篇nature的配色
ggplot(data=plotCoxoutput,aes(x=HR,y=gene,color=pvalue))+
  geom_errorbarh(aes(xmax=upper,xmin=lower),color="black",height=0,size=1.2)+
  geom_point(aes(x=HR,y=gene),size=3.5,shape=18)+ #画成菱形
  geom_vline(xintercept = 1,linetype='dashed',size=1.2)+
  scale_x_continuous(breaks = c(0.75,1,1.30))+
  coord_trans(x='log2')+ 
  ylab("Gene")+  #标签
  xlab("Hazard ratios of TDERS in ccRCC of OS")+ 
  labs(color="P value",title ="Univariate Cox regression analysis" )+
  scale_color_viridis()+  #nature配色
  theme_ggstatsplot()+  #好看的主题，同原文一致
  theme(panel.grid =element_blank()) #去除网格线
ggsave('coxplot.pdf',width = 8,height = 8)
ggsave('coxplotPFI.pdf',width = 8,height = 8)


####Riskscore相关性分析蝴蝶图228####
setwd("~/vip39/TCGA/exosome_ccRCC/supplementaryresult/")
library(data.table)
library(GSVA)
library(ggplot2)
library(ggcor)
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


# 加载TCGA-KIRC表达谱FPKM
#expr <- read.table("tcga_blca_fpkm.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
load("/home/data/refdir/database/tcga_counts_fpkm_tpm/TCGA-KIRC_tpm_gene_symbol.Rdata")
expr <-as.data.frame(tpms)
range(expr)
expr[1:4,1:4]
#tumsam <- colnames(expr)[substr(colnames(expr),11,13) == "01A"] # 取出肿瘤样本
# 取出RISKSCORE 分数
#siglec15 <- log2(sub$riskscore + 1)
siglec15 <- scale(sub$gpr_score)

##plod2表达分数 进行分析
siglec15 <- log2( expr["PLOD2",rownames(data)]+ 1)

dim(siglec15)
siglec15 <-data.frame(SIGLEC15=as.numeric(siglec15),
                      row.names = sub$Sample)
siglec15 <- as.data.frame(t(siglec15))
head(siglec15)
setwd("~/vip39/TCGA/exosome_ccRCC/supplementaryresult/PLOD2/")
write.csv(siglec15, "easy_input_gene.csv")

dim(immPath.score) #[1]  18 535

# 加载immunotherapy-predicted pathways
immPath <- read.table("~/vip39/huitu/FigureYa228linkCor/immunotherapy predicted pathways.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
immPath.list <- list()
for (i in rownames(immPath)) {
  tmp <- immPath[i,"Genes"]
  tmp <- toupper(unlist(strsplit(tmp,", ",fixed = T)))
  tmp <- gsub(" ","",tmp)
  immPath.list[[i]] <- tmp
}

# 加载膀胱癌相关签名
blcaPath <- read.table("~/vip39/huitu/FigureYa228linkCor/bladder cancer signature.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
blcaPath.list <- list()
for (i in rownames(blcaPath)) {
  tmp <- blcaPath[i,"Genes"]
  tmp <- toupper(unlist(strsplit(tmp,", ",fixed = T)))
  tmp <- gsub(" ","",tmp)
  blcaPath.list[[i]] <- tmp
}

# 计算immunotherapy-predicted pathways的富集得分
expr[1:4,1:4]

# sub=data
# sub$Sample <- rownames(sub)
immPath.score <- gsva(expr = as.matrix(log2(expr[,sub$Sample] + 1)),
                      immPath.list, 
                      method = "ssgsea")
# 保存到文件，便于套用
write.table(immPath.score, "easy_input_immPath.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 计算膀胱癌相关签名的富集得分
blcaPath.score <- gsva(expr = as.matrix(log2(expr[,sub$Sample] + 1)),
                       blcaPath.list,
                       method = "ssgsea")
# 保存到文件，便于套用
write.table(blcaPath.score, "easy_input_blcaPath.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# 从gmt文件提取感兴趣的通路
# 例如提取带有`METABOLISM`字样的通路
gset <- gmt2list("~/vip39/huitu/FigureYa228linkCor/c5.all.v7.4.symbols.gmt") 
###读取所有通路的文件
gset <- gmt2list("/home/data/t030432/vip39/database/GSEA/human/msigdb.v2023.1.Hs.symbols.gmt") 
gset.list <- gset[names(gset) %like% "METABOLISM"]

gset.list <- gset[names(gset) %like% "METABOLISM"]
gset.list <- gset.list[names(gset.list) %like% "KEGG"]
# 计算通路富集得分
gset.score <- gsva(expr = as.matrix(log2(expr[,sub$Sample] + 1)),
                   IOBR::signature_tumor,
                   method = "ssgsea")
# 保存到文件
write.table(gset.score, "easy_input_gset.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

library(IOBR)
IOBR::signature_metabolism
length(IOBR::signature_metabolism) #114
length(IOBR::signature_tme)#[1] 119
length(IOBR::signature_tumor) #[1] 16
length(IOBR::signature_collection)  #264
# 加载单个基因的表达量，也就是“1对多”的“1”
siglec15 <- read.csv("easy_input_gene.csv", row.names = 1, check.names = F)

# 加载第一组富集得分
immPath.score <- read.table("easy_input_immPath.txt", check.names = F)
# 跟目标基因SIGLEC15的表达量合并
immPath.score <- rbind.data.frame(immPath.score,
                                  siglec15)
dim(immPath.score)


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
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#F18D00",midpoint=0) +
  remove_axis("x")
p1

c('#4D4398','#F18D00')
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
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#F18D00",midpoint=0) +
  remove_axis("x")
p2
ggsave(filename = "ggcor plot in top right.pdf", width = 10,height = 8)



# 加载第3组富集得分
blcaPath.score <- read.table("easy_input_gset.txt", check.names = F)
# 跟目标基因SIGLEC15的表达量合并
blcaPath.score <- rbind.data.frame(blcaPath.score,
                                   siglec15)
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

p3 <- quickcor(t(blcaPath.score), 
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
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#F18D00",midpoint=0) +
  remove_axis("x")
p3
ggsave(filename = "ggcor plot in top right2.pdf", width = 10,height = 8)

####TIP读入
stepscore <- fread("~/vip39/database/TIP/KIRC/ssGSEA.normalized.score.txt",header = T,data.table = F)

head(stepscore[,1:4])
stepscore <- stepscore%>%
  column_to_rownames("Steps")
colnames(stepscore) <- str_sub(colnames(stepscore),1,15)
rownames(stepscore)[1:3] <- c("Step1.Release of cancer cell antigens",
                              "Step2.Cancer antigen presentation",
                              "Step3.Priming and activation")
rownames(stepscore)[c(21,22,23)] <- c("Step5.Infiltration of cancer cells into tumors",
                                      "Step6.Recognitioin of cancer cells by T cells",
                                      "Step7.Killing of cancer cells")
stepscore[1:4,1:4]
colnames(stepscore) <- paste0(colnames(stepscore),"A")
blcaPath.score <- stepscore[,sub$Sample]
blcaPath.score <- rbind.data.frame(blcaPath.score,
                                   siglec15)
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

p3 <- quickcor(t(blcaPath.score), 
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
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#F18D00",midpoint=0) +
  remove_axis("x")
p3

ggsave(filename = "ggcor plot in top rightstepscores.pdf", width = 10,height = 8)



####经典富集分析 GO term####
###使用easyTCGA版本
library(msigdbr)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

# gene_entrezid <- merge(gene_entrezid,diff_limma,by.x = "SYMBOL", by.y = "genesymbol")
# genelist <- gene_entrezid$logFC
# names(genelist) <- gene_entrezid$ENTREZID
# genelist <- sort(genelist,decreasing = T)

head(genelist)
gsea_res <- GSEA(geneList, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 seed = 456
)
gsea_res_symbol <- setReadable(gsea_res,"org.Hs.eg.db",keyType = "ENTREZID")
library(GseaVis)

terms <- gsea_res@result$ID[1:4]

gseaplot_list <- lapply(terms, function(x){
  gseaNb(object = gsea_res,
         geneSetID = x,
         termWidth = 30,
         addPval = T,
         pvalX = 0.75,
         pvalY = 0.6
  )
})

# 可以直接拼
cowplot::plot_grid(plotlist=gseaplot_list, ncol = 2)


terms <- gsea_res@result$ID[57:60]

gseaplot_list <- lapply(terms, function(x){
  gseaNb(object = gsea_res,
         geneSetID = x,
         termWidth = 30,
         addPval = T,
         pvalX = 0.75,
         pvalY = 0.6
  )
})

# 可以直接拼
cowplot::plot_grid(plotlist=gseaplot_list, ncol = 2)

terms <- gsea_res@result$ID[1:3]

gseaNb(object = gsea_res,
       geneSetID = terms,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05,pvalY = 0.05)

gseaNb(object = gsea_res,
       geneSetID = "GOMF_OLFACTORY_RECEPTOR_ACTIVITY")

#####reactome 分析补充上去
#aaa<-aaa[,-1]
#aaa<-na.omit(aaa)
#aaa$logFC<-sort(aaa$logFC,decreasing = T)
#geneList = aaa[,2]
#names(geneList) = as.character(aaa[,1])
geneList
#GSEA分析——GO
Go_gseresult <- gseGO(geneList, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont="all",  minGSSize = 10, maxGSSize = 500, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.5)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#波浪图
ridgeplot(KEGG_gseresult,10) #输出前十个结果
ridgeplot(Go_gseresult,10) #输出前十个结果
ridgeplot(Go_Reactomeresult,10)


gseaplot2(Go_Reactomeresult, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

gseaplot2(Go_Reactomeresult, geneSetID = 1:5, pvalue_table = TRUE,
          color = paletteer_c("scico::berlin", n = 5), ES_geom = "dot")

# go <- enrichGO(gene = geneList, OrgDb = "org.Hs.eg.db", ont="all")
# barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")
gseaplot2(Go_Reactomeresult, geneSetID = c("R-HSA-381753","R-HSA-418555","R-HSA-1461957",
                                           "R-HSA-8851680","R-HSA-176974"), pvalue_table = TRUE,
          color = paletteer_c("scico::berlin", n = 5), ES_geom = "dot")

###GO 分析
library(dplyr)#数据清洗
library(stringr)
library(org.Hs.eg.db)#物种注释(Homo sapiens)
library(clusterProfiler)#富集分析
library(ggplot2)#个性化绘图
library(RColorBrewer)#配色调整
#根据研究目的筛选差异基因(仅上调、下调或者全部)：
#yellow_gene <-colnames(dd)#所有差异基因
#ID转换：
#查看可转换的ID类型：
columns(org.Hs.eg.db)

##使用clusterProfiler包自带ID转换函数bitr(基于org.Hs.eg.db)：
#diff：

diff_entrez <- bitr(geneID = gene,
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
#save(GO_MF_diff,GO_CC_diff,GO_BP_diff,GO_all_diff,file = "GO_diff_DEG.rds")

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
KEGG_top20 <- na.omit(KEGG_top20)
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
#MF2 <- MF
#MF2$Description <- str_trunc(MF$Description,width = 50,side = "right")
#MF2$Description

MF$term <- factor(MF$Description,levels = rev(MF$Description))
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
p1 <- GO_bar("MF")+scale_fill_distiller(palette = "Blues",direction = 1)
p1

#CC:
p2 <- GO_bar("CC")+scale_fill_distiller(palette = "Reds",direction = 1)
p2

#BP:
p3 <- GO_bar("BP")+scale_fill_distiller(palette = "Oranges",direction = 1)
p3


######单细胞补充分析####
library(Seurat)
library(monocle)
library(DDRTree)
library(CellChat)
library(SCENIC)

sce <- readRDS("~/vip39/baiduyunpan/BNTL9/all_cluster.rds")
#sce <- BNTL9
sce <- SeuratObject::UpdateSeuratObject(sce)
###只能addmoudlescore
####使用其中一组的分数 TCGA-KIRC
DefaultAssay(sce) <- "RNA"
DefaultAssay(sce) <- "integrated"
sce@assays$integrated[1:4,1:4]
dim(sce@assays$integrated)
dim(sce@assays$RNA)
colnames(meta)
###细胞分群
Idents(sce) <- "celltype"
DimPlot(sce,reduction = "tsne",group.by = "celltype",label = T)

DimPlot(sce,reduction = "umap",group.by = "celltype",label = T)
DimPlot(sce,reduction = "pca",group.by = "clusterAnn",label = T)

###marker细胞的鉴定 top3  绘图


###美化
meta <- sce@meta.data
meta$celltypemino <- str_split_fixed(meta$clusterAnn,":",2)[,2]

col_df = data.frame(name=unique(meta$celltype)) %>% 
  arrange(name)
col_df
mycol = c(
  "#409079","#52a5c1","#c65341","#425785","#d6873b","#92b8da","#9f2b39",
  "#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84","#817f7e"
)

mycol =setNames(mycol,col_df$name)
UMAP <- sce@reductions[["umap"]]@cell.embeddings%>%as.data.frame()
meta <- cbind(meta,UMAP)
mid_coord_type = meta %>% 
  dplyr::select(c('celltype','UMAP_1','UMAP_2')) %>% 
  group_by(celltype) %>% 
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
pumap1<-ggplot(meta,aes( UMAP_1 ,UMAP_2))+
  geom_point(size=0.01,aes(color=celltype))+
  geom_text(data = mid_coord_type,size=5,
            aes(UMAP_1,UMAP_2,label=celltype))+
  theme_classic()+
  theme(legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values = mycol)

png(paste0(fig_path,'umap_maintype2.png'),width = 8,height = 6,
    units = 'in',res = 300)
print(pumap1)
dev.off()

###tsne绘图
TSNE <- sce@reductions[["tsne"]]@cell.embeddings%>%as.data.frame()
head(TSNE)

meta <- cbind(meta,TSNE)
mid_coord_type = meta %>% 
  dplyr::select(c('celltype','tSNE_1','tSNE_2')) %>% 
  group_by(celltype) %>% 
  summarise(
    TSNE_1 = median(tSNE_1),
    TSNE_2 = median(tSNE_2)
  )
head(meta)

pTSNE1<-ggplot(meta,aes( tSNE_1 ,tSNE_2))+
  geom_point(size=0.01,aes(color=celltype))+
  geom_text(data = mid_coord_type,size=5,
            aes(TSNE_1,TSNE_2,label=celltype))+
  theme_classic()+
  theme(legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values = mycol)

png(paste0(fig_path,'TSNE_maintype2.png'),width = 8,height = 6,
    units = 'in',res = 300)
print(pTSNE1)
dev.off()



table(sce$orig.ident)
DimPlot(sce, 
        reduction = "umap",
        group.by = "orig.ident")
sce$Diagnosis <- ifelse(str_detect(sce$orig.ident,"P"),"Primary","Metastasis")

DimPlot(sce, 
        reduction = "umap", 
        group.by = "Diagnosis", # metadata的Diagnosis列
        cols = c("royalblue3", "yellow2"),
        pt.size = 0.1)
table(sce2$clusterAnn)
length(unique(sce2$clusterAnn))
# 定义足够多的颜色，可以自己改
cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
DimPlot(sce, 
        reduction = "umap", 
        group.by = "clusterAnn", # metadata的Diagnosis列
        cols = cb_palette,label = T,label.size = 2,
        pt.size = 0.1)

df <- data.frame(x=1:4,y=2:5,z=rep(1:2,2))

# 先看看list和as.list函数的结果是什么样的
as.list(df) # 每一列对应list的一个元素
list(df) # 一整个数据框成为list的一个元素
split(df, 1:4) # 每一行作为list的一个元素
split(df, df$z) # 按照z列进行分组


###制备gmt
data1 <- data.frame(exo=gene)
head(data1)
exosomelist <- as.list(data1)
names(exosomelist) <- "exoccRCC"
tme <- exosomelist
sink("./exo.gmt")
yourlist <- tme
for (i in 1:length(yourlist)){
  cat(names(yourlist)[i])
  cat('\tNA\t')
  cat(paste(yourlist[[i]], collapse = '\t'))
  cat('\n')
}
sink()

GSET.FILEexom <- "/home/data/t030432/vip39/TCGA/exosome_ccRCC/supplementaryresult/COX/exo.gmt"

DimPlot(sce,raster = T)
colnames(sce@meta.data)
FeaturePlot(sce,features = "risk_score")
VlnPlot(sce,features = "risk_score")
VlnPlot(sce,features = "PLOD2")


####featureplot#####
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
gene<-unique(gene)
i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = sce, 
                                  colors_use = viridis_magma_dark_high, 
                                  features = gene[i])+NoLegend()+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
library(patchwork)
fig_path="~/vip39/TCGA/exosome_ccRCC/supplementaryresult/singlecell/"
p<-wrap_plots(plots, ncol = 4);p
png(paste0(fig_path,'umap_12pic.png'),width = 12,height = 9,
    units = 'in',res = 300)
print(p)
dev.off()

ggsave(p,file="~/vip39/TCGA/exosome_ccRCC/supplementaryresult/singlecell/featureplot2.pdf",width = 12,height = 9)

###再来个dotplot
library(stringr)  
genes_to_check=str_to_upper(unique(gene))
genes_to_check

th=theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(sce, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype')  + coord_flip()+th
p_all_markers

ggsave(p_all_markers, filename="check_all_marker_by_seurat_celltype.pdf",height = 7,width = 9)

####提取数据####
data<-p_all_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)

####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features");p


p1 = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black"),
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features")

p1
# 保存
library(patchwork)
p+p1
ggsave(file="dotplot_celltype.pdf",height= 10,width = 14)
######自己设计genelist进行打分#####
##制备gmt文件 list 直接addmoudule打分

##fig231
AUC
###
sce2 <- AddModuleScore(sce,features = exosomelist,name = "Riskscore")
colnames(sce2@meta.data)
VlnPlot(sce2,features = "Riskscore1")
FeaturePlot(sce2,features = "Riskscore1")



####换一个打分方式
library(VISION)
library(Seurat)
library(BiocParallel)
library(AUCell)
library(GSEABase)
library(GSVA)
library(GSVA)
library(UCell)
library(singscore)
library(magrittr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(scales)
library(methods)
library(purrr)
library(reshape2)
library(stringr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggridges)
library(viridis)
      
####gene表达dotplot绘制
marker_list=list()
marker_list[['B-cell']]= c('CD79A', 'MS4A1')
marker_list[['Dendritic']]=c('FCER1A', 'CST3')
marker_list[['Monocytes']]=c('FCGR3A')
marker_list[['NK']]=c('GNLY', 'NKG7')
marker_list[['Other']]= c('IGLL1')
marker_list[['Plasma']]= c('IGJ')
marker_list[['T-cell']]= c('CD3D')

sce_pbmc <- SetIdent(sce_pbmc, value = sce_pbmc@meta.data$bulk_labels)
p <- DotPlot(object = sce, features=gene)
p



DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(object = sce, normalization.method = "LogNormalize", scale.factor = 10000) 
sum(sce[["RNA"]]@data[,1])

#TPM <- 
# 这样就能得到单细胞中的TPM矩阵（每一列的和为10,000）
head(colSums(TPM))
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
# 默认情况会选择前2000个基因，nfeatures可以设置选择的基因数
top10 <- head(VariableFeatures(sce), 10)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)

coef=read.csv('~/vip39/TCGA/exosome_ccRCC/GPR_coef.csv',header = T,check.names = F)
coef=coef$Coef..boot_SD
###data就是单细胞表达矩阵 scale后的
sce[["RNA"]]@scale.data[1:4,1:4]

data <- sce[["RNA"]]@data
data <- as.data.frame(data)
range(data)
data[1:4,1:4]
data <- expm1(sce[["RNA"]]@data)
data <- log2(data+1)
data[1:4,1:4]
datasc <- data
data <- as.data.frame(t(data[gene,]))
table(gene%in%rownames(data))
data <- as.data.frame(datasc[gene,])

#coef=GPR_japan
#coef=coef$Coef..boot_SD

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}

data[1:4,1:4]
data=as.data.frame(t(data))

sce$risk_score=gpr_score
#save(sce,file = "~/vip39/TCGA/exosome_ccRCC/supplementaryresult/exoseurat.rds")
load("~/vip39/TCGA/exosome_ccRCC/supplementaryresult/exoseurat.rds")
gene <- c("COBLL1","PLOD2","MEGF9","EMCN","PLXNA2","ASRGL1","MICU1","TTC33","WDR11","UPB1","SLC5A12",
          "LSM14A")


signatures <- "/home/data/t030432/vip39/TCGA/exosome_ccRCC/supplementaryresult/COX/exo.gmt"
geneSets <- getGmt(signatures)

matrix <- Seurat::GetAssayData(object=sce, slot="counts", assay="RNA")
set.seed(1)
cells_rankings <- AUCell_buildRankings(matrix, nCores=10, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
signature_exp <- data.frame(t(getAUC(cells_AUC)))
sce2 <- Seurat::AddMetaData(sce, as.data.frame(signature_exp))
colnames(sce2@meta.data)
FeaturePlot(sce2,features="exoccRCC")
VlnPlot(sce2,features="exoccRCC")
####美化分数
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidydr)
library(ggsci)
library(paletteer)
#最常规的小提琴图
VlnPlot(sce2, features = "exoccRCC", group.by="celltype")+NoLegend()
VlnPlot(sce2, features = "exoccRCC", group.by="celltype",cols = color,flip = T)+NoLegend()
genes_to_check = c(gene,"exoccRCC")
#颜色
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
p1 <- VlnPlot(sce2, features = genes_to_check,
              group.by  = "celltype",
              flip = T,stack = T,cols = color
)
p2<-p1 + NoLegend()

ggsave(p2,file="vlnplot1.pdf",height = 8,width=6)
#分组
table(sce2$Diagnosis)
sce2$Diagnosis <- factor(sce2$Diagnosis,levels = c("Primary","Metastasis"))
#分组分半小提琴图
p<-VlnPlot(sce2, features = genes_to_check,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'Diagnosis',
           group.by = "celltype",
           cols = c("#78C2C4","#C73E3A"),
           split.plot = T)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')
p
ggsave(p,file="vlnplot2.pdf",height = 8,width=6)

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
####分数分布  正式分析 开始了
##fig290 bargraph#####
# 提取细胞分类信息
# Sample：由GSE151914_metadata.txt.gz提供的样本信息
# 规律为：样本组织来源_样本测序时使用通道_样本分组，如Tum_963_WT
# GSE151914 Series Matrix File(s)显示963对应R1，650对应R2
# 使用正则表达式(gsub)进行将Sample列进行分割
# 使用按位置替换(substr)或字符切割(strsplit)也可起到类似结果
colnames(sce2@meta.data)

cellinfo <- FetchData(sce2, vars = c("orig.ident", "seurat_clusters"))

cellinfo$Run <- gsub("(\\w+)(_)(\\d+)(_)(\\w+)", "\\3", cellinfo$Sample)
cellinfo$Run <- ifelse(test = cellinfo$Run %in% c("963", "650"),
                       yes = "R1", no = "R2")
cellinfo <- sce2@meta.data
cellinfo$Group <- ifelse(sce2$exoccRCC>median(sce2$exoccRCC),"High","Low")
cellinfo$Group <- sce
str(cellinfo)
str(cellinfo)
# 设置分组颜色
cell.col <- setNames(object = c("#1A63A8", "#FC6910"),
                     nm = c("WT", "KO"))

# 制作列联表
plot.data <- as.data.frame(table(cellinfo$seurat_clusters, cellinfo$Run, cellinfo$Group))

p1 <- ggplot(plot.data, aes(x = Var2, y = Freq, fill = Var3)) +
  geom_bar(stat = "identity", position = position_fill(reverse = T)) +                    # 绘制堆叠比例柱状图，将WT和KO的顺序倒过来
  scale_fill_manual(values = cell.col) +                                                  # 设置不同组对应的颜色
  facet_wrap(~Var1, nrow = 1) +                                                           # 设置柱状图按seurat_cluster分别显示
  theme_classic() +                                                                       # 去除不必要的元素
  xlab("Replicate") +                                                                     # 修改x和y轴列名
  ylab("Fraction of Cells") +                           
  theme(strip.background = element_rect(fill="grey", color = "white", size = 1),          # 设置顶部seurat_cluster显示为白框灰底
        strip.text = element_text(size = 12, colour="black"),                             # 设置顶部seurat_cluster字号和文字颜色
        legend.title = element_blank(),                                                   # 去除图例标题
        axis.title = element_text(size = 15))                                             # 调整坐标轴标题字号
p1


library(ggplot2)
library(rstatix)
# 构建seurat_clusters和Group的列联表
tbl <- table(cellinfo$seurat_clusters, cellinfo$Group)
tbl
tbl <- table(cellinfo$celltype, cellinfo$Group)
tbl

tbl <- table(cellinfo$clusterAnn, cellinfo$Group)
tbl
# 进行fisher精确检验和post hoc检验
# post hoc检验：对各组样本进行 one vs other的fisher检验，进行多重性校正，得到各组的p-adj
fisher.test(tbl, simulate.p.value = T) # fisher精确检验

post.hoc <- row_wise_fisher_test(tbl) # post hoc检验
post.hoc$sig.label <- ifelse(test = post.hoc$p.adj.signif == "ns", # 调整显著性显示标签
                             yes = "", no = post.hoc$p.adj.signif) # 不显示NS (No Significant)
str(post.hoc)

# 绘制柱状图
plot.data <- as.data.frame(tbl)
head(plot.data)

# 设置分组颜色
cell.col <- setNames(object = c("#1A63A8", "#FC6910"),
                     nm = c("Low", "High"))

p2 <- ggplot(plot.data, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity", position = position_fill(reverse = T)) +  # 绘制堆叠比例柱状图，将WT和KO的顺序倒过来
  scale_fill_manual(values = cell.col) +                                # 设置不同组对应的颜色
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +                   # 设置y坐标轴的刻度
  geom_text(data = post.hoc,                                            # 设置顶部的显著性标签
            aes(x = group, y = 1.1, label = sig.label),                 # 标签的位置和具体内容
            inherit.aes = F) +                                          
  xlab("Cluster") +                                                     # 修改x和y轴列名
  ylab("Fraction of Cells") +
  theme_classic() +                                                     # 去除不必要的元素
  theme(legend.title = element_blank(),                                 # 去除图例标题
        axis.text.x = element_text(colour = 'black',size = 10,angle = 45),
        axis.text = element_text(size = 10),                            # 调整坐标轴刻度字号大小
        axis.title = element_text(size = 15))                           # 调整坐标轴标题字号大小
p2

ggsave(filename = "p2.pdf", width = 8, height = 4)

###细胞gene分布

###cellchat分析#####
##fig267
library(CellChat)
save(sce2,file = "~/vip39/TCGA/exosome_ccRCC/supplementaryresult/singlecell/sceforcellchat.rds")

#pbmc3k里的seurat_annotations有一些NA注释，过滤掉
data.input = sce2@assays$RNA@data
meta.data =  sce2@meta.data
meta.data = meta.data[!is.na(meta.data$seurat_annotations),]
meta.data$Riskgroup=ifelse(meta.data$exoccRCC>median(meta.data$exoccRCC),"High","Low")
data.input = data.input[,row.names(meta.data)]

#设置因子水平
table(meta.data$celltype)
meta.data$seurat_annotations = factor(meta.data$celltype,
                                      levels = c("Mast cells", "Endothelial cells", "B cells", "T cells", "Neutrophils", 
                                                 "Epithelial cells (cancer cells)", "NK cells", "Macrophages", "Dendritic cells",
                                                 "Fibroblasts","Erythroblasts"))
meta.data$Riskgroup = factor(meta.data$Riskgroup,levels = c("Low","High"))
### 1.2 Create a CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "seurat_annotations")
### 1.3 可在cellchat对象的meta插槽中添加表型信息
# 添加meta.data信息
#cellchat <- addMeta(cellchat, meta = meta.data)

# 设置默认的labels
levels(cellchat@idents) # show factor levels of the cell labels
table(cellchat@meta$exoccRCC)

cellchat <- setIdent(cellchat, ident.use = "Riskgroup") 
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
### 1.4 加载CellChat受配体数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat.low@DB <- CellChatDB.use
### 1.5 对表达数据进行预处理，用于细胞间通讯分析
# subset the expression data of signaling genes for saving computation cost
cellchat.low <- subsetData(cellchat.low) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 10) # do parallel

cellchat.low <- identifyOverExpressedGenes(cellchat.low)
cellchat.low <- identifyOverExpressedInteractions(cellchat.low)
cellchat.low <- computeCommunProb(cellchat.low,population.size = F)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.low <- filterCommunication(cellchat.low, min.cells = 10)




cellchat.low = subsetCellChat(cellchat ,idents.use = c("Low"))
cellchat.high= subsetCellChat(cellchat ,idents.use = c("High"))

meta.data1=meta.data[meta.data$Riskgroup=="Low",]
data.input1 = data.input[,row.names(meta.data1)]
### 1.2 Create a CellChat object
cellchat.low <- createCellChat(object = data.input1, 
                           meta = meta.data1, 
                           group.by = "seurat_annotations")

meta.data2=meta.data[meta.data$Riskgroup=="High",]
data.input2 = data.input[,row.names(meta.data2)]
### 1.2 Create a CellChat object
cellchat.high <- createCellChat(object = data.input2, 
                               meta = meta.data2, 
                               group.by = "seurat_annotations")


object.list <- list(Low = cellchat.low, High = cellchat.high)
cellchat2 <- mergeCellChat(object.list, add.names = names(object.list))
cellchat2

gg1 <- compareInteractions(cellchat2, show.legend = F, group = c(1,2),color.use = c("#1A63A8", "#FC6910"))
gg2 <- compareInteractions(cellchat2, show.legend = F, group = c(1,2), color.use = c("#1A63A8", "#FC6910"),measure = "weight")
gg1 + gg2



#差异性circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat2, weight.scale = T,comparison = c(1,2))
netVisual_diffInteraction(cellchat2, weight.scale = T, measure = "weight",comparison = c(1,2))
#差异性heatmap
gg1 <- netVisual_heatmap(cellchat2,comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat2, measure = "weight",comparison = c(1,2))
gg1 + gg2
# 总览性circle plot
object.list <- list(group1 = cellchat.low, group2 = cellchat.high)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))


rankSimilarity(cellchat2, type = "functional",comparison1 =c(1,2),comparison2 =c(1,2))

gg1 <- rankNet(cellchat2, mode = "comparison", stacked = T, color.use = c("#1A63A8", "#FC6910"),do.stat = TRUE,comparison = c(1,2))
gg2 <- rankNet(cellchat2, mode = "comparison", stacked = F, color.use = c("#1A63A8", "#FC6910"), do.stat = TRUE,comparison = c(1,2))
gg1 + gg2


#run netAnalysis_computeCentrality
object.list = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})

object.list2 = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})

# merge data
cellchat3 <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)



library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

##outgoing
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height =13,font.size = 5);ht1
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height =13,font.size = 5)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

##incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(object.list)[i],
                                        width = 5, height =13,font.size = 5, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = pathway.union, title = names(object.list)[i+1],
                                        width = 5, height =13,font.size = 5, color.heatmap = "GnBu")
#ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", 
#                                        signaling = pathway.union, title = names(object.list)[i+2],
#                                        width = 5, height = 6, color.heatmap = "GnBu")
#draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

##all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", 
                                        signaling = pathway.union, title = names(object.list)[i], 
                                        width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all",
                                        signaling = pathway.union, title = names(object.list)[i+1], 
                                        width = 5, height = 6, color.heatmap = "OrRd")
#ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", 
#                                        signaling = pathway.union, title = names(object.list)[i+2], 
#                                        width = 5, height = 6, color.heatmap = "OrRd")
#draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat2, sources.use = 4, targets.use = c(1:11),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat2, sources.use = 4, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Low risk", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat2, sources.use = 4, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Low risk", angle.x = 45, remove.isolate = T)
gg1 + gg2



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Low"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat2 <- identifyOverExpressedGenes(cellchat2, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1,
                                       thresh.fc = 0.1,
                                       thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat2, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat2, net = net,
                              datasets = "Low",ligand.logFC = 0.2,
                              receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat2, net = net,
                                datasets = "High",
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat2)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat2)

pairLR.use.up = net.up[,"interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat2,
                        pairLR.use = pairLR.use.up,
                        sources.use = 4, targets.use = c(5:11),
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat2,
                        pairLR.use = pairLR.use.down,
                        sources.use = 4, targets.use = c(5:11),
                        comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

#####不同分组代谢改变 特定细胞类型#######





