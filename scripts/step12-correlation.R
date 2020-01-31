
rm(list=ls())
load(file = './Rdata/TCGA_KIRC_mut.Rdata')
load(file = './Rdata/TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')

table(group_list)
load(file='./Rdata/survival_input.Rdata')

dat=data.frame(gene1=log2(exprSet['hsa-mir-10b',]+1),
               gene2=log2(exprSet['hsa-mir-143',]+1),
               stage=phe$stage)
save(dat,file = 'for_scatter.Rdata')
library(ggpubr)
# google search : ggpubr boxplot add p-value
# http://www.sthda.com/english/rpkgs/ggpubr/reference/stat_cor.html 
dat$stage=as.factor(dat$stage)
sp <- ggscatter(dat, x = "gene1", y = "gene2",
                add = "reg.line",  # Add regressin line 
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) 
# Add correlation coefficient
sp
sp + stat_cor(method = "pearson", label.x = 15, label.y = 20)

# Color by groups and facet
#::::::::::::::::::::::::::::::::::::::::::::::::::::
sp <- ggscatter( dat, x = "gene1", y = "gene2",
                color = "stage", palette = "jco",
                add = "reg.line", conf.int = TRUE)
sp + stat_cor(aes(color = stage),label.x = 15 )


