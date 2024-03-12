
#save.image(file = "~/TCGA/exosome_ccRCC/exoalldata.rds")
load("~/vip39/TCGA/exosome_ccRCC/exoalldata.rds")

setwd("~/TCGA/exosome_ccRCC/validation/")
load("~/TCGA/exosome_ccRCC/validation/exsomevalidation.rds")
load("~/TCGA/exosome_ccRCC/exoalldata.rds")
##开始构建模型
#lasso联合bostrap
#参考 生信作曲家课程
#install.packages('boot')
library(boot)
###########寻找预后相关的G蛋白偶联受体##############

library(survival)      #引用包
pFilter= 0.05         #显著性过滤条件

library(limma)               #引用包

# 数据处理
## 载入TCGA胃癌表达谱
# load('STAD_tpm.Rdata')
# gpr=read.table('G_protein.txt')
# gpr=as.data.frame(t(gpr))
# gpr=gpr$V1
# 获取GPR基因
gpr=finalgene

# GPR表达谱的差异分析
KIRC_tpm <- log(KIRC_tpm+1)
kirc.expr[1:4,1:4]
data=KIRC_tpm[rownames(KIRC_tpm) %in% gpr,]
## 肿瘤和正常样品

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))
library(limma)
design=model.matrix(~ group_list)

fit=lmFit(data,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=0.05) 
write.csv(allDiff,file ='allDiff.csv',quote = F)


## 更新data
data=data[rownames(allDiff),]
#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

## 去除在大多数样本中都表达为0的基因
keep <- rowSums(data>0) >= floor(0.75*ncol(data))
table(keep)
## 30
data<- data[keep,]

data=as.data.frame(t(data))


# 单因素寻找预后影响的GPR
## 读取生存数据
#suv=read.table('TCGA-STAD.survival.tsv',row.names = 1,header = T,check.names = F)

cli=dplyr::select(kirc.sinfo,'survival_time','fustat')
colnames(cli)=c("futime", "fustat")

## 数据合并并输出结果
data1 <- data%>%
  rownames_to_column("ID")%>%
  select(ID,everything())
data1$ID <-   str_sub(data1$ID,1,12)
data1 <- data1[!duplicated(data1$ID),]
rownames(data1) <- data1$ID
data1 <- data1[,-1]
data <- data1
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
#你需要牢记这样的结构，用于cox
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
## 写出GPR的基础矩阵
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)

rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件


##单因素cox分析
outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  set.seed(123456)
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP))
    print(coxP)
  }
}


#输出单因素结果
write.table(outTab,file="uniCox_gpr.txt",sep="\t",row.names=F,quote=F)
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp_gpr.txt",sep="\t",row.names=F,quote=F)


#### 热图
##########热图
#正常和肿瘤数目、
#load('STAD_tpm.Rdata')
gene=read.table('uniCox_gpr.txt',header = T)
gene=gene$gene

data=KIRC_tpm[gene,]
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)

conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目

sampleType=ifelse(group=='1',1,2)

identical(colnames(data),colnames(KIRC_tpm))

#基因差异分析
sigVec=c()
allDiff=read.csv('allDiff.csv',header = T,row.names = 1)
alldiff_cox=allDiff[gene,]
pvalue=alldiff_cox$adj.P.Val
Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
sigVec=paste0(gene, Sig)


## 另起一个compare矩阵，避免破坏原来的
compare=data
# 修饰一下行名
row.names(compare)=sigVec

#调整顺序，保证出图美观
normal=compare[,sampleType==1]
tumor=compare[,sampleType==2]
compare=cbind(normal,tumor)

#热图可视化
Type=c(rep("Normal",conNum), rep("Tumor",treatNum))
names(Type)=colnames(compare)
Type=as.data.frame(Type)
library(pheatmap)

pheatmap::pheatmap(compare,
                   #color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
                   color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   annotation=Type,
                   breaks = c(seq(-3,3,length.out = 100)),
                   cluster_cols =F,
                   cluster_rows =T,
                   scale="row",
                   show_colnames=F,
                   show_rownames=T,
                   fontsize=6,
                   fontsize_row=7,
                   fontsize_col=6)



################先lasso,进行bootstrap_multicox回归#####################
gene=read.table('uniCox_gpr.txt',header = T)
gene=gene$gene
# 没有包先安装
#install.packages('survival')
library(survival)
library(survminer)
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt=rt[,c('futime','fustat',gene)]




## 先lasso筛基因!!!!!!!
set.seed(202209)   #设定随机种子 
rt <- rt[rt$futime>30,]
x=as.matrix(rt[,c(3:ncol(rt))]) 
y=data.matrix(Surv(rt$futime,rt$fustat)) 

# 没包先安装
# install.packages('glmnet')
library(glmnet)
fit=glmnet(x, y, family = "cox", alpha = 1) 
plot(fit, xvar = "lambda", label = F)

cvfit = cv.glmnet(x, y, family="cox",nfolds = 10,alpha=1) 
plot(cvfit) 
#其中两条虚线分别指示了两个特殊的λ值 
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

coef =coef(fit,s = cvfit$lambda.min)
index = which(coef !=0)
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef   #查看模型的相关系数

## 剩下基因
gene=read.table('uniCox_gpr.txt',header = T)
gene=gene[gene$gene %in% geneCoef[,1],]

write.table(gene,file = 'uniCox_lasso_gpr.txt',quote = F,sep = '\t',row.names = F)



########进行bootstrap_multicox回归####
#install.packages('boot')
library(boot)
gene=read.table('uniCox_lasso_gpr.txt',header = T)
# boot_coef=coef/Boot_sd
gene=gene$gene
library(survival)
library(survminer)
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt=rt[,c('futime','fustat',gene)]
rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt)
ggforest(cox)

#install.packages('boot')
library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
head(GPR)
write.csv(GPR,file= 'GPR_coef.csv',quote = F)

# 构建GPRscore
gene=read.table('uniCox_lasso_gpr.txt',header = T)
gene=gene$gene
data=KIRC_tpm[gene,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

# 读取系数
coef=read.csv('GPR_coef.csv',header = T,check.names = F)
coef=coef$Coef..boot_SD
gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据
#suv=read.table('TCGA-STAD.survival.tsv',row.names = 1,header = T,check.names = F)
colnames(kirc.sinfo)
cli=dplyr::select(kirc.sinfo,'survival_time','fustat')
colnames(cli)=c("futime", "fustat")
head(cli)
rownames(cli) <- paste0(rownames(cli),"-01A")
##数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]


## K-M生存分析
rt=cbind(cli,data)
rt$futime=rt$futime/30
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)

### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt$gpr_score), "Low", "High")
data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
save(data,file = "~/TCGA/exosome_ccRCC/tcgasub.rds")

diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    

save(data,file ='GPRscore_and_group.Rdata')


View(sub)
##ROC 和time ROC
#c-index
library(pROC)
library(timeROC)
riskRoc <- timeROC(T = sub$futime/12,delta = sub$fustat,
                   marker = sub$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




####在多个数据集中进行验证####
####japan####
load("/home/data/vip39/database/KIRC_immune_therapy/japan_cohort/japan_RCCall.rds")
library(boot)
gene=read.table('uniCox_lasso_gpr.txt',header = T)
# boot_coef=coef/Boot_sd
gene=gene$gene
library(survival)
library(survminer)
class(japan_mRNA_fpkm)

japan.expr2 <- as.data.frame(t(japan_mRNA_fpkm[gene,]))
colnames(japan_clin)

japan.cli <- japan_clin[,c("month","outcome")]
colnames(japan.cli) <- c("futime", "fustat")

comsam <- intersect(rownames(japan.expr2),rownames(japan.cli))
japan.expr2 <- japan.expr2[comsam,]
japan.cli <- japan.cli[comsam,]
rt2 <- cbind(japan.cli,japan.expr2)
rt2$fustat <- ifelse(rt2$fustat=="alive",0,1)
#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt2=rt2[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
head(rt2)

cox=coxph(Surv(futime, fustat) ~.,data = rt2)
ggforest(cox)

#install.packages('boot')
library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt2, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_japan=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_japan,file= 'GPR_coef_japan.csv',quote = F)

# 构建GPRscore


# 读取系数

coef=read.csv('GPR_coef_japan.csv',header = T,check.names = F)
coef=GPR_japan
coef=coef$Coef..boot_SD
class(rt2)
colnames(rt2)
data <- as.data.frame(t(rt2[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt2$gpr_score <- data$gpr_score
rt2$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)

### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt2$gpr_score), "Low", "High")
data=rt2
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    



riskRoc <- timeROC(T = rt2$futime/12,delta = rt2$fustat,
                   marker = rt2$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果

ROC <- timeROC(T = rt2$futime/12,   
               delta = rt2$fustat,   
               marker = rt2$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

##DCA
as.numeric(str_sub(japan_clin[rownames(rt2),"Stage"],3,3))
rt2$Stage <- as.numeric(factor(str_split_fixed(KIRC_cli[rownames(rt2),"pathologic_stage"]," ",2)[,2]))
rt2$Stage <-as.numeric(str_sub(japan_clin[rownames(rt2),"Stage"],3,3))
rt2$Fuhrman <- japan_clin[rownames(rt2),"Fuhrman"]
rt2$Age <- japan_clin[rownames(rt2),"Age"]
rt2$cancer <- rt2$fustat==1
rt2$ttcancer <- rt2$futime/12 ##年份

##取出不要的数据
rt2 <- rt2[rt2$Fuhrman!="undetermined",]
rt2$Fuhrman <- as.numeric(rt2$Fuhrman)

# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                    data = rt2)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt2)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt2)
Fuhrman <- coxph(Surv(ttcancer, cancer) ~ Fuhrman, data = rt2)
df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)

##准备候选基因和临床信息 标准输入矩阵即可
###cancer cell ####
load("/home/data/vip39/database/KIRC_immune_therapy/cancercellpaper/cancercellRCC.rds")
table(Cancercell_cli$ARM)
# atezo_bev sunitinib 
# 407       416 
table(Cancercell_cli$OBJECTIVE_RESPONSE)
# CR  NE  PD  PR  SD 
# 32  65 151 257 318
#完全缓解（CR）
#部分缓解（PR）
#疾病进展（PD）
#病稳定（SD）
#无法评估 （NE）
#控制率CR+PR+SD
32+318+257

Sunitinib_arm <- Cancercell_cli[Cancercell_cli$ARM=="sunitinib",]
Sunitinib_arm$risk_score <- rt3[rownames(Sunitinib_arm),"gpr_score"]
Sunitinib_arm <- Sunitinib_arm[Sunitinib_arm$OBJECTIVE_RESPONSE!="NE",]

Sunitinib_arm$outcome <- ifelse(Sunitinib_arm$OBJECTIVE_RESPONSE=="PD", "Poor","Good")

##ROC
# 先来之前的
library(pROC)
roc1 <- roc(Sunitinib_arm$outcome, Sunitinib_arm$risk_score,data=Sunitinib_arm,aur=TRUE,
            levels=c("Good", "Poor"),smooth=T,ci=T,
            boot.n=100) 

plot(roc1,print.auc=T,ci=T)
# 计算sensitivity(se)的CI,根据100次bootstrap，每0.1取个点
se.obj <- ci(roc1, of="se", boot.n=100,specificities=seq(0, 1, 0.1)) 

plot(se.obj,type='bars',col = 'black',print.auc=T) 
## 获取矩阵
roc1_df=data.frame(TPR=roc1$sensitivities,FPR=1-roc1$specificities)
roc1_df=roc1_df[order(roc1_df$TPR,decreasing = T),]
se.obj=as.data.frame(se.obj)

se.obj$FPR=1-as.numeric((rownames(se.obj)))
se.obj$TPR=se.obj$`50%`

g <- ggplot() + 
  geom_line(data=roc1_df,aes(x =FPR, y = TPR),size=0.6,color='#0072b5') + 
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size = 1,linetype=6 )+ 
  labs(x = "False positive rate", y = "Ture positive rate", title ="ROC curve")+
  geom_point(aes(x =FPR, y = TPR), size=3,data = 
               se.obj,shape=21,color='white',stroke=2,fill='#0072b5') + 
  geom_errorbar(aes(ymin=`2.5%`,ymax=`97.5%`,x=FPR,y=TPR),data = se.obj, 
                width=.03,size=1,color='#0072b5')+
  annotate("text",x = .75, y = .25,label = paste("Riskscore 
AUC=",round(roc1$auc,2)),color = "#0073c2",size=4)+
  theme_bw()

g




atezo_bev_arm <- Cancercell_cli[Cancercell_cli$ARM=="atezo_bev",]
atezo_bev_arm$risk_score <- rt3[rownames(atezo_bev_arm),"gpr_score"]
##ROC


atezo_bev_arm <- atezo_bev_arm[atezo_bev_arm$OBJECTIVE_RESPONSE!="NE",]

atezo_bev_arm$outcome <- ifelse(atezo_bev_arm$OBJECTIVE_RESPONSE=="PD", "Poor","Good")
atezo_bev_arm$risk_group <- ifelse(atezo_bev_arm$risk_score>median(atezo_bev_arm$risk_score),"high","low")
##ROC
# 先来之前的
library(pROC)
roc1 <- roc(atezo_bev_arm$outcome, atezo_bev_arm$risk_score,data=atezo_bev_arm,aur=TRUE,
            levels=c("Good", "Poor"),smooth=T,ci=T,
            boot.n=100) 

plot(roc1,print.auc=T,ci=T)
# 计算sensitivity(se)的CI,根据100次bootstrap，每0.1取个点
se.obj <- ci(roc1, of="se", boot.n=100,specificities=seq(0, 1, 0.1)) 

plot(se.obj,type='bars',col = 'black',print.auc=T) 
## 获取矩阵
roc1_df=data.frame(TPR=roc1$sensitivities,FPR=1-roc1$specificities)
roc1_df=roc1_df[order(roc1_df$TPR,decreasing = T),]
se.obj=as.data.frame(se.obj)

se.obj$FPR=1-as.numeric((rownames(se.obj)))
se.obj$TPR=se.obj$`50%`

g <- ggplot() + 
  geom_line(data=roc1_df,aes(x =FPR, y = TPR),size=0.6,color='#0072b5') + 
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size = 1,linetype=6 )+ 
  labs(x = "False positive rate", y = "Ture positive rate", title ="ROC curve")+
  geom_point(aes(x =FPR, y = TPR), size=3,data = 
               se.obj,shape=21,color='white',stroke=2,fill='#0072b5') + 
  geom_errorbar(aes(ymin=`2.5%`,ymax=`97.5%`,x=FPR,y=TPR),data = se.obj, 
                width=.03,size=1,color='#0072b5')+
  annotate("text",x = .75, y = .25,label = paste("Riskscore 
AUC=",round(roc1$auc,2)),color = "#0073c2",size=4)+
  theme_bw()

g




##sunitinib 组进行预测



##atezo_bev组进行预测


rt_cancercell <- as.data.frame(t(Cancercell_tpm[gene,]))
colnames(Cancercell_cli)
Cancercell_cli$futime <- as.numeric(Cancercell_cli$PFS_MONTHS)*30
Cancercell_cli$fustat <- ifelse(Cancercell_cli$PFS_CENSOR=="TRUE",1,0)
cancercell_cli <- Cancercell_cli[,c("futime","fustat")]


comsam <- intersect(rownames(rt_cancercell),rownames(cancercell_cli))
rt_cancercell <- rt_cancercell[comsam,]
cancercell_cli <- cancercell_cli[comsam,]
rt3 <- cbind(cancercell_cli,rt_cancercell)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt3=rt3[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt3)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt3, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_cancercell=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_cancercell,file= 'GPR_coef_cancercell.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_cancercell.csv',header = T,check.names = F)
coef <- GPR_cancercell
coef=coef$Coef..boot_SD
class(rt3)
colnames(rt3)
data <- as.data.frame(t(rt3[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt3$gpr_score <- data$gpr_score
rt3$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt3$gpr_score), "Low", "High")
data=rt3
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    
library(timeROC)
riskRoc <- timeROC(T = rt3$futime/12,delta = rt3$fustat,
                   marker = rt3$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt3$futime/12,   
               delta = rt3$fustat,   
               marker = rt3$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

##校准

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt2)
rt2$futime
range(rt2$futime)
class(rt2$time)
rt2$time <- rt2$futime*30

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt2,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=30,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt2,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f5, cmethod="KM", method="boot",u=730,m=30,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt2,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=30,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt2,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=30,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt2,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()


###美化ROC

##DCA


##cancerll DCA
colnames(Cancercell_cli)

rt3$IMDC <- Cancercell_cli[rownames(rt3),"IMDC_RISK_SCORE"]
rt3$MSKCC <- Cancercell_cli[rownames(rt3),"MSKCC_RISK_SCORE"]
rt3$Age <- as.numeric(Cancercell_cli[rownames(rt3),"AGE"])
rt3$TMB <- as.numeric(Cancercell_cli[rownames(rt3),"TMB"])

rt3$cancer <- rt3$fustat==1
rt3$ttcancer <- rt3$futime/30 ##月份

#rt3 <- rt3[rt3$TCGA_cluster!="Other",]
rt3 <- na.omit(rt3)
table(rt3$IMDC)
table(rt3$MSKCC)

rt3$IMDC <- as.numeric(factor(rt3$IMDC,levels = c("POOR","Intermediate","Favorable")))
rt3$MSKCC <- as.numeric(factor(rt3$MSKCC,levels = c("Low","Intermediate","High")))

# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                    data = rt3)
IMDC <- coxph(Surv(ttcancer, cancer) ~ IMDC, data = rt3)
MSKCC <- coxph(Surv(ttcancer, cancer) ~ MSKCC, data = rt3)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt3)
TMB <- coxph(Surv(ttcancer, cancer) ~ TMB, data = rt3)
range(rt3$ttcancer)

df3 <- ggDCA::dca(Risk_score,IMDC,MSKCC,Age,TMB,
                  times = c(1,2,3,6,12,18)
)
library(ggsci)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)



####GSE29609 不合适 缺少gene   ####
load("/home/data/vip39/database/KIRC_immune_therapy/GSE29609/GSE29609_RCCall.rds")
library(boot)
gene=read.table('~/TCGA/exosome_ccRCC/uniCox_lasso_gpr.txt',header = T)
# boot_coef=coef/Boot_sd
gene=gene$gene
library(survival)
library(survminer)

GSE29609.expr2 <- as.data.frame(t(GSE29609_exp[gene,]))
colnames(GSE29609_cli)

GSE29609_cli <- GSE29609_cli[,c("futime", "fustat")]
colnames(GSE29609_cli) <- c("futime", "fustat")

comsam <- intersect(rownames(GSE29609.expr2),rownames(GSE29609_cli))

GSE29609.expr2 <- GSE29609.expr2[comsam,]
GSE29609_cli <- GSE29609_cli[comsam,]
rt2 <- cbind(GSE29609_cli,GSE29609.expr2)
#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt2=rt2[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt2)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt2, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_japan=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_japan,file= 'GPR_coef_japan.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_japan.csv',header = T,check.names = F)
coef=coef$Coef..boot_SD
class(rt2)
colnames(rt2)
data <- as.data.frame(t(rt2[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt2$gpr_score <- data$gpr_score
rt2$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt2$gpr_score), "Low", "High")
data=rt2
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    

riskRoc <- timeROC(T = rt2$futime/12,delta = rt2$fustat,
                   marker = rt2$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果







###免疫治疗队列####

#####checkmate####
load("/home/data/vip39/database/KIRC_immune_therapy/PD1/CheckmateRCC.rds")


#os
Checkmate_clin$futime <- as.numeric(Checkmate_clin$OS)*30
Checkmate_clin$fustat <- ifelse(Checkmate_clin$OS_CNSR=="0",0,1)

Checkmate_clin <- Checkmate_clin[,c("futime","fustat")]

Checkmate.expr2 <- as.data.frame(t(Checkmate_exp[gene,]))
colnames(Checkmate_clin)

Checkmate_cli <- Checkmate_clin[,c("futime", "fustat")]
colnames(Checkmate_cli) <- c("futime", "fustat")

comsam <- intersect(rownames(Checkmate.expr2),rownames(Checkmate_cli))

Checkmate.expr2 <- Checkmate.expr2[comsam,]
Checkmate_cli <- Checkmate_cli[comsam,]
rt4 <- cbind(Checkmate_cli,Checkmate.expr2)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt4=rt4[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt4)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt4, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_Checkmate=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_Checkmate,file= 'GPR_coef_Checkmate.csv',quote = F)

# 读取系数
coef=read.csv('GPR_coef_Checkmate.csv',header = T,check.names = F)
coef <- GPR_Checkmate
coef=coef$Coef..boot_SD
class(rt4)
colnames(rt4)
data <- as.data.frame(t(rt4[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt4$gpr_score <- data$gpr_score
rt4$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt4$gpr_score), "Low", "High")
data=rt4
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    
library(timeROC)
riskRoc <- timeROC(T = rt4$futime/12,delta = rt4$fustat,
                   marker = rt4$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt4$futime/12,   
               delta = rt4$fustat,   
               marker = rt4$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

##DCA
colnames(Checkmate_clin)
rt4$IMDC <- Checkmate_clin[rownames(rt4),"IMDC"]
rt4$MSKCC <- Checkmate_clin[rownames(rt4),"MSKCC"]
rt4$Age <- as.numeric(Checkmate_clin[rownames(rt4),"Age"])
rt4$cancer <- rt4$fustat==1
rt4$ttcancer <- rt4$futime/30 ##月份

#rt4 <- rt4[rt4$TCGA_cluster!="Other",]
rt4 <- na.omit(rt4)
table(rt4$IMDC)
table(rt4$MSKCC)
rt4 <- rt4[rt4$IMDC!="NA",]
rt4 <- rt4[rt4$IMDC!="NOT REPORTED",]
rt4 <- rt4[rt4$MSKCC!="NOT REPORTED",]

rt4$IMDC <- as.numeric(factor(rt4$IMDC,levels = c("POOR","INTERMEDIATE","FAVORABLE")))
rt4$MSKCC <- as.numeric(factor(rt4$MSKCC,levels = c("POOR","INTERMEDIATE","FAVORABLE")))

# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                    data = rt4)
IMDC <- coxph(Surv(ttcancer, cancer) ~ IMDC, data = rt4)
MSKCC <- coxph(Surv(ttcancer, cancer) ~ MSKCC, data = rt4)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt4)

df3 <- ggDCA::dca(Risk_score,IMDC,MSKCC,Age,
                  times = c(1,2,3,6,12,18)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)
##校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt4)
rt4$futime
range(rt4$futime)
class(rt4$time)
rt4$time <- rt4$futime

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt4,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=50,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt4,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f5, cmethod="KM", method="boot",u=730,m=50,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt4,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=50,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt4,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=50,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt4,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()


####JAVELIN验证####
load("/home/data/vip39/database/KIRC_immune_therapy/JAVELIN_Renal_101_trial/JAVELIN_RCCall.rds")
JAVELIN_cli$futime <- as.numeric(JAVELIN_cli$PFS_P)*30
JAVELIN_cli$fustat <- JAVELIN_cli$PFS_P_CNSR

JAVELIN_cli <- JAVELIN_cli[,c("futime","fustat")]

JAVELIN.expr2 <- as.data.frame(t(JAVELIN_exp[gene,]))
colnames(JAVELIN_cli)
colnames(JAVELIN_cli) <- c("futime", "fustat")

comsam <- intersect(rownames(JAVELIN.expr2),rownames(JAVELIN_cli))

JAVELIN.expr2 <- JAVELIN.expr2[comsam,]
JAVELIN_cli <- JAVELIN_cli[comsam,]
rt5 <- cbind(JAVELIN_cli,JAVELIN.expr2)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt5=rt5[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt5)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt5, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_JAVELIN=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_JAVELIN,file= 'GPR_coef_JAVELIN.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_JAVELIN.csv',header = T,check.names = F)
coef <- GPR_JAVELIN
coef=coef$Coef..boot_SD
class(rt5)
colnames(rt5)
data <- as.data.frame(t(rt5[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt5$gpr_score <- data$gpr_score
rt5$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt5$gpr_score), "Low", "High")
data=rt5
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    
library(timeROC)
riskRoc <- timeROC(T = rt5$futime/12,delta = rt5$fustat,
                   marker = rt5$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt5$futime/12,   
               delta = rt5$fustat,   
               marker = rt5$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

#校准


library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt5)
rt5$futime
range(rt5$futime)
class(rt5$time)
rt5$time <- rt5$futime

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt5,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=50,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt5,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f5, cmethod="KM", method="boot",u=730,m=50,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt5,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=50,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt5,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=50,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt5,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()


plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year"), #图例文字
       col =c("#2166AC","#B2182B"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框

###DCA
colnames(JAVELIN_cli)
rt5$Stage <- as.numeric(factor(str_split_fixed(KICH_cli[rownames(rt5),"pathologic_stage"]," ",2)[,2]))
rt5$PDL1 <- JAVELIN_cli[rownames(rt5),"PDL1FL"]

rt5$TCGA_cluster <- JAVELIN_cli[rownames(rt5),"TCGA_cluster"]
rt5$Age <- JAVELIN_cli[rownames(rt5),"AGE"]
rt5$cancer <- rt5$fustat==1
rt5$ttcancer <- rt5$futime/30 ##月份
rt5 <- rt5[rt5$TCGA_cluster!="Other",]
rt5 <- na.omit(rt5)
rt5$TCGA_cluster <- as.numeric(factor(rt5$TCGA_cluster))
# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                    data = rt5)
#Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt5)
TCGA_cluster <- coxph(Surv(ttcancer,cancer) ~ TCGA_cluster, data = rt5)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt5)

df3 <- ggDCA::dca(Risk_score,TCGA_cluster,Age,
                  times = c(1,2,3,6,12)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)
###GSE22541####
#没有其他信息 无法做DCA
load("/home/data/vip39/database/KIRC_immune_therapy/GSE22541/GSE22541_RCCall.rds")

GSE22541_cli$futime <- as.numeric(GSE22541_cli$futime)
GSE22541_cli$fustat <- GSE22541_cli$fustat
GSE22541_exp <- GSE22541_exp[,rownames(GSE22541_cli)]

GSE22541_cli <- GSE22541_cli[,c("futime","fustat")]

GSE22541.expr2 <- as.data.frame(t(GSE22541_exp[gene,]))
colnames(GSE22541_cli)
colnames(GSE22541_cli) <- c("futime", "fustat")

comsam <- intersect(rownames(GSE22541.expr2),rownames(GSE22541_cli))

GSE22541.expr2 <- GSE22541.expr2[comsam,]
GSE22541cli <- GSE22541_cli[comsam,]
rt6 <- cbind(GSE22541_cli,GSE22541.expr2)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt6=rt6[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt6)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt6, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_GSE22541=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_GSE22541,file= 'GPR_coef_GSE22541.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_GSE22541.csv',header = T,check.names = F)
coef <- GPR_GSE22541
coef=coef$Coef..boot_SD
class(rt6)
colnames(rt6)
data <- as.data.frame(t(rt6[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt6$gpr_score <- data$gpr_score
rt6$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt6$gpr_score), "Low", "High")
data=rt6
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    

library(timeROC)
riskRoc <- timeROC(T = rt6$futime/12,delta = rt6$fustat,
                   marker = rt6$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt6$futime/365,   
               delta = rt6$fustat,   
               marker = rt6$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

###校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt6)
rt6$futime
range(rt6$futime)
class(rt6$time)
rt6$time <- rt6$futime

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt6,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=10,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt6,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f5, cmethod="KM", method="boot",u=730,m=10,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt6,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=10,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt6,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=10,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt6,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()



###MTAB3267####

load("/home/data/vip39/database/KIRC_immune_therapy/E-MTAB-3267/MTAB3267_RCCall.rds")

MTAB3267_cli$futime <- as.numeric(MTAB3267_cli$`Characteristics[progression free survival]`)*30
MTAB3267_cli$fustat <- MTAB3267_cli$`Characteristics[progression]`
MTAB3267_exp <- MTAB3267_exp[,rownames(MTAB3267_cli)]
table(gene%in%rownames(MTAB3267_exp))
# FALSE  TRUE 
# 2    10 

###Fudancohort####
load("/home/data/vip39/database/KIRC_immune_therapy/Fudan_RCC/FudanRCCall.rds")
group <- str_detect(colnames(FDRCC_exp2),"TA")

FDRCC_exptumor <- FDRCC_exp2[,!str_detect(colnames(FDRCC_exp2),"TA")]
rownames(FDRCC_clin) <- paste0(rownames(FDRCC_clin),"_T")
table(paste0(rownames(FDRCC_clin),"_T")%in%colnames(FDRCC_exptumor))
interid_FD <- intersect(rownames(FDRCC_clin),colnames(FDRCC_exptumor))
FDRCC_exp <- FDRCC_exptumor[,interid_FD]
FDRCC_clin <- FDRCC_clin[interid_FD,]

FDRCC_clin$futime <- as.numeric(FDRCC_clin$`OS(month)`)*30
FDRCC_clin$fustat <- ifelse(FDRCC_clin$`Live Status`=="living",0,1)
table(gene%in%rownames(FDRCC_exp))
# FALSE  TRUE 
# 5     7
####TCGA_KIRC####
library(data.table)
load("/home/data/refdir/database/tcga_counts_fpkm_tpm/TCGA-KIRC_tpm_gene_symbol.Rdata")
tpms <- as.data.frame(tpms)
tpms[1:4,1:4]
KIRC_exp <- tpms[,str_sub(colnames(tpms),14,15)=="01"]
colnames(KIRC_exp)=str_sub(colnames(KIRC_exp),1,12)
KIRC_cli <- fread("/home/data/refdir/database/TCGA临床数据/肾透明细胞癌/肾透明细胞癌.txt",header = T,data.table = F)
rownames(KIRC_cli) <- KIRC_cli$id
KIRC_cli <- KIRC_cli[colnames(KIRC_exp),]

KIRC_cli$futime <- as.numeric(KIRC_cli$survival_time)
KIRC_cli$fustat <- ifelse(KIRC_cli$status=="Alive",0,1)

KIRC_cli <- KIRC_cli[,c("futime","fustat")]

KIRC.expr2 <- as.data.frame(t(KIRC_exp[gene,]))
rownames(KIRC.expr2) <- str_replace_all(rownames(KIRC.expr2),"\\.","\\-")

comsam <- intersect(rownames(KIRC.expr2),rownames(KIRC_cli))

KIRC.expr2 <- KIRC.expr2[comsam,]
KIRC_cli <- KIRC_cli[comsam,]
rt9 <- cbind(KIRC_cli,KIRC.expr2)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt9=rt9[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt9)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt9, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_KIRC=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_KIRC,file= 'GPR_coef_KIRC.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_KIRC.csv',header = T,check.names = F)
coef <- GPR_KIRC
coef=coef$Coef..boot_SD
class(rt9)
colnames(rt9)
data <- as.data.frame(t(rt9[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt9$gpr_score <- data$gpr_score
rt9$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt9$gpr_score), "Low", "High")
data=rt9
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


rt9 <- rt9[rt9$futime!=0,]
class(rt9)
class(rt9$futime)
range(rt9$futime)

library(timeROC)
riskRoc <- timeROC(T = rt9$futime/365,delta = rt9$fustat,
                   marker = rt9$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt9$futime/365,   
               delta = rt9$fustat,   
               marker = rt9$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

##DCA
rt9$Stage <- as.numeric(factor(str_split_fixed(KIRC_cli[rownames(rt9),"pathologic_stage"]," ",2)[,2]))
rt9$Age <- KIRC_cli[rownames(rt9),"age"]
rt9$cancer <- rt9$fustat==1
rt9$ttcancer <- rt9$futime/365 ##年份

# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                    data = rt9)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt9)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt9)

df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)

##校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt9)
rt9$futime
rt9$time <- rt9$futime

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

###这个地方差点出错了 是f2!!!!!
#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt9,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.2,1),ylim= c(0.2,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()



###TCGA_KIRP####
library(data.table)
load("/home/data/refdir/database/tcga_counts_fpkm_tpm/TCGA-KIRP_tpm_gene_symbol.Rdata")
tpms <- as.data.frame(tpms)
tpms[1:4,1:4]
KIRP_exp <- tpms[,str_sub(colnames(tpms),14,15)=="01"]
colnames(KIRP_exp)=str_sub(colnames(KIRP_exp),1,12)
KIRP_cli <- fread("/home/data/refdir/database/TCGA临床数据/乳头状肾细胞癌/乳头状肾细胞癌.txt",header = T,data.table = F)
rownames(KIRP_cli) <- KIRP_cli$id
KIRP_cli <- KIRP_cli[colnames(KIRP_exp),]

KIRP_cli$futime <- as.numeric(KIRP_cli$survival_time)
KIRP_cli$fustat <- ifelse(KIRP_cli$status=="Alive",0,1)

KIRP_cli <- KIRP_cli[,c("futime","fustat")]

KIRP.expr2 <- as.data.frame(t(KIRP_exp[gene,]))


comsam <- intersect(rownames(KIRP.expr2),rownames(KIRP_cli))

KIRP.expr2 <- KIRP.expr2[comsam,]
KIRP_cli <- KIRP_cli[comsam,]
rt7 <- cbind(KIRP_cli,KIRP.expr2)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt7=rt7[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt7)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt7, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_KIRP=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_KIRP,file= 'GPR_coef_KIRP.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_KIRP.csv',header = T,check.names = F)
coef <- GPR_KIRP
coef=coef$Coef..boot_SD
class(rt7)
colnames(rt7)
data <- as.data.frame(t(rt7[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt7$gpr_score <- data$gpr_score
rt7$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt7$gpr_score), "Low", "High")
data=rt7
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


rt7 <- rt7[rt7$futime!=0,]
class(rt7)
class(rt7$futime)
range(rt7$futime)

library(timeROC)
riskRoc <- timeROC(T = rt7$futime/10950,delta = rt7$fustat,
                   marker = rt7$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt7$futime/10950,   
               delta = rt7$fustat,   
               marker = rt7$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

##DCA
rt7$Stage <- as.numeric(factor(str_split_fixed(KIRP_cli[rownames(rt7),"pathologic_stage"]," ",2)[,2]))
rt7$Age <- KIRP_cli[rownames(rt7),"age"]
rt7$cancer <- rt7$fustat==1
rt7$ttcancer <- rt7$futime/10950 ##年份

# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                    data = rt7)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt7)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt7)

df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)
###校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt7)
rt7$futime
rt7$time <- rt7$futime/30

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt7,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt7,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f5, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt7,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt7,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt7,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.2,1),ylim= c(0.2,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()


###KICH####

library(data.table)
load("/home/data/refdir/database/tcga_counts_fpkm_tpm/TCGA-KICH_tpm_gene_symbol.Rdata")
tpms <- as.data.frame(tpms)
tpms[1:4,1:4]
KICH_exp <- tpms[,str_sub(colnames(tpms),14,15)=="01"]
colnames(KICH_exp)=str_sub(colnames(KICH_exp),1,12)
KICH_cli <- fread("/home/data/refdir/database/TCGA临床数据/肾嫌色细胞癌/肾嫌色细胞癌.txt",header = T,data.table = F)
rownames(KICH_cli) <- KICH_cli$id
KICH_cli <- KICH_cli[colnames(KICH_exp),]

KICH_cli$futime <- as.numeric(KICH_cli$survival_time)
KICH_cli$fustat <- ifelse(KICH_cli$status=="Alive",0,1)

KICH_cli <- KICH_cli[,c("futime","fustat")]

KICH.expr2 <- as.data.frame(t(KICH_exp[gene,]))


comsam <- intersect(rownames(KICH.expr2),rownames(KICH_cli))

KICH.expr2 <- KICH.expr2[comsam,]
KICH_cli <- KICH_cli[comsam,]
rt8 <- cbind(KICH_cli,KICH.expr2)

#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt8=rt8[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
cox=coxph(Surv(futime, fustat) ~.,data = rt8)
##KICH有两列的表达值为0 需要去除  "UPB1"    "SLC5A12"
rt8 <- rt8[,-c(12,13)]
cox=coxph(Surv(futime, fustat) ~.,data = rt8)
ggforest(cox)

library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt8, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_KICH=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_KICH,file= 'GPR_coef_KICH.csv',quote = F)

# 读取系数

coef=read.csv('GPR_coef_KICH.csv',header = T,check.names = F)
coef <- GPR_KICH
coef=coef$Coef..boot_SD
class(rt8)
colnames(rt8)
data <- as.data.frame(t(rt8[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt8$gpr_score <- data$gpr_score
rt8$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)
### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt8$gpr_score), "Low", "High")
data=rt8
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


rt8 <- rt8[rt8$futime!=0,]
class(rt8)
class(rt8$futime)
range(rt8$futime)

library(timeROC)
riskRoc <- timeROC(T = rt8$futime/365,delta = rt8$fustat,
                   marker = rt8$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt8$futime/365,   
               delta = rt8$fustat,   
               marker = rt8$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

###多个时间点 多个模型比较  Grade Stage Age riskscore
#library(rlang,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1")
library(ggDCA)
library(dcurves)
library(survival)
library(ggsci)
colnames(rt8)
colnames(KICH_cli)
#rt8$Stage <- KICH_cli[rownames(rt8),"pathologic_stage"]
str_split_fixed(KICH_cli[rownames(rt8),"pathologic_stage"]," ",2)[,2]
as.numeric(factor(str_split_fixed(KICH_cli[rownames(rt8),"pathologic_stage"]," ",2)[,2]))
rt8$Stage <- as.numeric(factor(str_split_fixed(KICH_cli[rownames(rt8),"pathologic_stage"]," ",2)[,2]))
rt8$Age <- KICH_cli[rownames(rt8),"age"]
rt8$cancer <- rt8$fustat==1
rt8$ttcancer <- rt8$futime/365 ##年份
# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score, 
                  data = rt8)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt8)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt8)

df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)
###矫正曲线绘制


#单因素多因素 tabel绘制
#install.packages("gtsummary")
.libPaths(c("/home/data/refdir/Rlib/","/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1"))
#install.packages("autoReg",dependencies = T)
library(rlang,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.1")
library(gtsummary)
library(autoReg)
###校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt8)
rt8$futime
rt8$time <- rt8$futime

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt8,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt8,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f5, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt8,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt8,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt8,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.4,1),ylim= c(0.4,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.4,1),ylim= c(0.4,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.4,1),ylim= c(0.4,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.4,1),ylim= c(0.4,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()



#cox和logistics两种 单因素多因素

###IGGC####
load("/home/data/vip39/database/KIRC_validation_data/ICGC/ICGC_KIRC_input.rds")
library(boot)
gene=read.table('uniCox_lasso_gpr.txt',header = T)
# boot_coef=coef/Boot_sd
gene=gene$gene
library(survival)
library(survminer)
class(ICGC_fpkm)

ICGC_fpkm <- log2(ICGC_fpkm+1)
ICGC.expr2 <- as.data.frame(t(ICGC_fpkm[gene,]))
colnames(ICGC_cli)

ICGC.cli <- ICGC_cli[,c("OS","status")]
colnames(ICGC.cli) <- c("futime", "fustat")
head(ICGC.cli)

comsam <- intersect(rownames(ICGC.expr2),rownames(ICGC.cli))
ICGC.expr2 <- ICGC.expr2[comsam,]
ICGC.cli <- ICGC.cli[comsam,]
rt1 <- cbind(ICGC.cli,ICGC.expr2)
#rt1$fustat <- ifelse(rt1$fustat=="alive",0,1)
#rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt1=rt1[,c('futime','fustat',gene)]
#rt <- rt[rt$futime>30,]
# 初始HR
head(rt1)

cox=coxph(Surv(futime, fustat) ~.,data = rt1)
ggforest(cox)

#install.packages('boot')
library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，稍等待
set.seed(123456)
boot_results <- boot(data=rt1, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)

## 单因素不良预后的这边可能变好，但是因为是多因素，不管
print(boot_results)

## 获取参数
coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)

## 定义coef/sd为新的参数
ratio=coef/sd

GPR_ICGC=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(GPR_ICGC,file= 'GPR_coef_ICGC.csv',quote = F)

# 构建GPRscore
# 读取系数

coef=read.csv('GPR_coef_ICGC.csv',header = T,check.names = F)
coef=GPR_ICGC
coef=coef$Coef..boot_SD
class(rt1)
colnames(rt1)
data <- as.data.frame(t(rt1[,-c(1,2)]))

gpr_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  gpr_score=c(gpr_score,score)
}
data[1:4,1:4]
data=as.data.frame(t(data))
data$gpr_score=gpr_score

#读取生存数据

## K-M生存分析
rt1$gpr_score <- data$gpr_score
rt1$futime
#rt <- rt[rt$futime>0,]
library(survival)
library(survminer)

### 中位值划分
Type=ifelse(data[,'gpr_score']<= median(rt1$gpr_score), "Low", "High")
data=rt1
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='GPR_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Days)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    



riskRoc <- timeROC(T = rt1$futime/365,delta = rt1$fustat,
                   marker = rt1$gpr_score,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果

ROC <- timeROC(T = rt1$futime/365,   
               delta = rt1$fustat,   
               marker = rt1$gpr_score,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p

##DCA

#as.numeric(str_sub(ICGC_clin[rownames(rt1),"Stage"],3,3))
#rt1$Stage <- as.numeric(factor(str_split_fixed(KIRC_cli[rownames(rt1),"pathologic_stage"]," ",2)[,2]))
colnames(ICGC_cli)
rt1$Stage <-as.numeric(str_sub(ICGC_cli[rownames(rt1),"T_stage"],2,2))
#rt1$Fuhrman <- ICGC_clin[rownames(rt1),"Fuhrman"]
rt1$Age <- ICGC_cli[rownames(rt1),"age"]
rt1$cancer <- rt1$fustat==1
rt1$ttcancer <- rt1$futime/365 ##年份

##取出不要的数据
#rt1 <- rt1[rt1$Fuhrman!="undetermined",]
#rt1$Fuhrman <- as.numeric(rt1$Fuhrman)

# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ gpr_score + Stage, 
                    data = rt1)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt1)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt1)
#Fuhrman <- coxph(Surv(ttcancer, cancer) ~ Fuhrman, data = rt1)
df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)

###校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt1)
rt1$futime
rt1$time <- rt1$futime

##1年
f1 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt1,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt1,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt1,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  gpr_score,data=rt1,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt1,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.2,1),ylim= c(0.2,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
###和其他预后模型进行比较####


