rm(list=ls())
require(maftools) 
options(stringsAsFactors = F) 
# 理论上可以针对任何癌症的突变数据结果进行可视化和适当的统计。
laml = read.maf(maf = './GDC/TCGA.KIRC.mutect.somatic.maf.gz')
laml 
project='TCGA_KIRC'
tmp = laml@data
laml@data=laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),]
#laml@data=laml@data[!grepl('^MUC',laml@data$Hugo_Symbol),]
length(unique(laml@data$Tumor_Sample_Barcode))
#有336个病人

#了解基因、样本的突变情况
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)  
#可视化
plotmafSummary(maf = laml, rmOutlier = TRUE,showBarcodes = FALSE,
               addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#突变频谱图
oncoplot(maf = laml, top = 30, fontSize = 1)


#条形图是每个样本突变个数的可视化
va = laml@variants.per.sample
#热图每一行是一个基因在所有样本里的突变情况
dplyr::count(tmp,Hugo_Symbol,sort = T)
#结果显示VHL在169样本中突变，所以是49%
#条形图是突变种类的可视化:
table(tmp[tmp$Hugo_Symbol=="VHL",]$Variant_Classification)

tmp %>% filter(Hugo_Symbol=="VHL") %>%
  count(Variant_Classification,sort = T)


#算了个啥统计量
getFields(laml)  
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
mut=laml@data[laml@data$Variant_Type == "SNP",c("Hugo_Symbol","Chromosome","Start_Position","Tumor_Sample_Barcode","t_vaf")]
mut$pos=paste(mut$Chromosome,mut$Start_Position,sep=':')

save(mut,file = '../Rdata/TCGA_KIRC_mut.Rdata')


