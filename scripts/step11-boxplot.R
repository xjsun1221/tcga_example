rm(list=ls())
#1.任意miRNA在tumor和normal样本中的表达量----
#输入数据只要表达矩阵即可
load(file = './Rdata/TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
table(group_list)

library(ggstatsplot)
dat = data.frame(gene = expr["hsa-mir-143",],
                 group = group_list)
ggbetweenstats(data = dat, x = group,  y = gene,title = "hsa-mir-143")

#2.任意miRNA在任意两个分组中的表达量对比----
#只要肿瘤样本,522个，只要是可以根据临床信息查到的分组，例如生死、人种、阶段，都可以拿来做分组。
#需要注意调整样本顺序，一一对应。
exprSet = expr[,as.numeric(substr(colnames(expr),14,15)) < 10]
library(stringr)
x1=str_sub(colnames(exprSet),1,12)
x2=str_to_upper(meta$patient.bcr_patient_barcode)

table(x2 %in% x1)
table(x1 %in% x2)
length(unique(x1))
#发现一个问题，样本的前12位代表病人的编号，列名是有重复的，为了一对一关系，去重一下
exprSet = exprSet[,!duplicated(str_sub(colnames(exprSet),1,12))]
x1=str_sub(colnames(exprSet),1,12)
x2=str_to_upper(meta$patient.bcr_patient_barcode)
table(x2 %in% x1)
meta = meta[x2 %in% x1 ,]
#按照生死、人种、阶段分组看看
table(meta$patient.vital_status)
table(meta$patient.stage_event.pathologic_stage)
table(meta$patient.race)

dat = data.frame(gene = exprSet["hsa-mir-143",],
                 vital_status = meta$patient.vital_status,
                 stage = meta$patient.stage_event.pathologic_stage,
                 race = meta$patient.race)
ggbetweenstats(data = dat, x = vital_status,  y = gene,title = "hsa-mir-143")
ggbetweenstats(data = dat, x = stage,  y = gene,title = "hsa-mir-143")
ggbetweenstats(data = dat, x = race,  y = gene,title = "hsa-mir-143")

#3.根据某个基因是否突变分组比较某miRNA的表达量----
load(file = './Rdata/TCGA_KIRC_mut.Rdata')
dim(exprSet)
head(mut)
length(unique(str_sub(mut$Tumor_Sample_Barcode,1,12)))
table(x1 %in% unique(str_sub(mut$Tumor_Sample_Barcode,1,12)))
#522个样本中有331个有突变信息记录,将这些样本对应的表达矩阵取出来。
expm = exprSet[,x1 %in% unique(str_sub(mut$Tumor_Sample_Barcode,1,12))]

VHL_mut=substr(as.character(
  as.data.frame( mut[mut$Hugo_Symbol=='VHL','Tumor_Sample_Barcode'])[,1] ),
  1,12)

library(dplyr)
mut  %>% 
  filter(Hugo_Symbol=='VHL')  %>%   
  as.data.frame()  %>% 
  pull(Tumor_Sample_Barcode)   %>%  
  as.character()   %>%   
  substr(1,12)


#false 是未突变样本，true是突变样本

tail(rownames(expm))
dat=data.frame(gene=log2(expm['hsa-mir-98',]),
               mut= substr(colnames(expm),1,12) %in% VHL_mut)

ggbetweenstats(data = dat, x = mut,  y = gene)

#可以换个方法，直接计算根据某基因突变与否划分的两组之间表达量是否有差异。
res.aov <- aov(gene ~ as.factor(mut), data = dat)
summary(res.aov)
TukeyHSD(res.aov)
summary(res.aov)[[1]]$`Pr(>F)`[1]
#可以批量计算，哗啦哗啦的。
