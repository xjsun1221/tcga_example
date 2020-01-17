rm(list=ls())
options(stringsAsFactors = F)

Rdata_dir='./Rdata/'
Figure_dir='./figures/'
# 加载上一步从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。
load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
table(group_list)


# 这里做生存分析，已经不需要正常样本的表达矩阵了，所以需要过滤。
# 而且临床信息，有需要进行整理。
### survival analysis only for patients with tumor.
if(F){
  exprSet=na.omit(expr)
  exprSet=exprSet[,group_list=='tumor']
  dim(exprSet)
  str(meta)
  colnames(meta)=c('ID','event','death','last_followup','race','age','gender','stage')
  #调整ID和表达矩阵内容和顺序都一样
  head(meta$ID)
  meta$ID=str_to_upper(meta$ID) 
  meta=meta[match(substr(colnames(exprSet),1,12),meta$ID),]
  head(meta$ID)
  head(colnames(exprSet))
  #1.整理临床信息
  #1.1由随访时间和死亡时间计算生存时间
  meta[,3][is.na(meta[,3])]=0
  meta[,4][is.na(meta[,4])]=0
  meta$days=as.numeric(meta[,3])+as.numeric(meta[,4])
  #时间以月份记
  meta$time=meta$days/30 
  #1.2 根据生死定义event，活着是0，死的是1
  meta$event=ifelse(meta$event=='alive',0,1)
  table(meta$event)
  #1.3 年龄和年龄分组
  meta$age=as.numeric(meta$age)
  meta$age_group=ifelse(meta$age>median(meta$age),'older','younger')
  table(meta$age_group)
  #1.4 stage
  library(stringr) 
  meta$stage=str_split(meta$stage,' ',simplify = T)[,2]
  table(meta$stage)
  #1.5 race
  table(meta$race)

  
  save(exprSet,meta,
       file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
      )
}
# 上面被关闭的代码，就是在整理临床信息和生存分析的表达矩阵。
# 整理好的数据，直接加载即可
load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)

#生存分析----
library(survival)
library(survminer)
#1.性别
sfit <- survfit(Surv(time, event)~gender, data=meta)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
#2.性别和年龄
## 多个 ggsurvplots作图生存曲线代码合并 
sfit1=survfit(Surv(time, event)~gender, data=meta)
sfit2=survfit(Surv(time, event)~age_group, data=meta)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2,pval =TRUE, data = meta, risk.table = TRUE)
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
dev.off()
# 可以很明显看到，肿瘤病人的生存受着诊断癌症的年龄的影响，却与性别无关。

#3.挑选感兴趣的（多个）基因做生存分析----
# 来自于文章：2015-TCGA-ccRCC-5-miRNAs-signatures
# Integrated genomic analysis identifies subclasses and prognosis signatures of kidney cancer
# miR-21,miR-143,miR-10b,miR-192,miR-183
gs=c('hsa-mir-21','hsa-mir-143','hsa-mir-192',
     'hsa-mir-183') 
splots <- lapply(gs, function(g){
  meta$gene=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(time, event)~gene, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 2, risk.table.height = 0.4)
dev.off()


#4.批量生存分析 使用  logrank test 方法----
#（就是批量计算了p值）
mySurv=with(meta,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  # gene=exprSet[1,]
  meta$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=meta)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p=sort(log_rank_p)
table(log_rank_p<0.01) #88个
# 可以看到，文章里面挑选出来的生存分析相关的miRNA基因，在我们的分析里面都是显著的。

c('hsa-mir-21','hsa-mir-143','hsa-mir-192',
  'hsa-mir-183','hsa-mir-10b')  %in% names(log_rank_p[log_rank_p<0.01])

library(pheatmap)
choose_gene=names(log_rank_p[log_rank_p<0.01])
choose_matrix=expr[choose_gene,]
source("3-plotfunction.R")

draw_heatmap(choose_matrix,group_list)
## 主成分分析
draw_pca(log(choose_matrix+1),group_list)
