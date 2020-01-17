rm(list=ls())
require(maftools) 
options(stringsAsFactors = F) 
load("forest.Rdata")
laml = read.maf(maf = './GDC/TCGA.KIRC.mutect.somatic.maf.gz')
laml 
project='TCGA_KIRC'

suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
options(stringsAsFactors = F)
mut=laml@data
head(mut)
mut=mut[mut$Variant_Type=='SNP',]
a=mut[,c(16,5,6,12,13)]
colnames(a)=c( "Sample","chr", "pos","ref",  "alt")

a$Sample=as.character(a$Sample)

plot(table(a$Sample),las=2)
#制作signatures的输入数据（96频谱）
sigs.input <- mut.to.sigs.input(mut.ref = a, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
class(sigs.input)
#第一个样本的突变绘图
barplot(as.numeric(sigs.input[1,])~colnames(sigs.input))
sigs.input[1:4,1:4]
if(F){
  w=lapply(unique(a$Sample), function(i){
    ## signatures.cosmic signatures.nature2013
    sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                               signatures.ref = signatures.cosmic, 
                               sample.id =  i, 
                               contexts.needed = TRUE,
                               tri.counts.method = 'default')
    print(c(i,which(unique(a$Sample)==i)))
    return(sample_1$weights)
  })
  w=do.call(rbind,w)
  save(w,file = paste0(project,"w.Rdata"))
}
load(paste0(project,"w.Rdata"))
library(pheatmap)
pheatmap(t(w),cluster_rows = F,cluster_cols = T,show_colnames = F)
#好看一点的热图,顺序还没调整
x = colnames(mat)[match(laml@variants.per.sample$Tumor_Sample_Barcode,colnames(mat))]
w = w[x,]
identical(x,as.character(laml@variants.per.sample$Tumor_Sample_Barcode))
library(ComplexHeatmap)

mat = t(w)
column_ha = HeatmapAnnotation(bar1 = anno_barplot(laml@variants.per.sample$Variants))
Heatmap(mat, name = "mat", 
        top_annotation = column_ha, 
        #right_annotation = row_ha,
        cluster_rows = F,
        cluster_columns = F,
        #show_row_names = F,
        show_column_names = F)

# Determine the signatures contributing to the two example samples
if(F){
  lapply(unique(a$Sample), function(i){
    i = "TCGA-B8-4622-01A-02D-1553-08"
    ## signatures.cosmic signatures.nature2013
    sample_1 = whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures.cosmic, 
                               sample.id =  i, 
                               contexts.needed = TRUE,
                               tri.counts.method = 'default')
    #pdf(paste0(i,'.sig.pdf'))
    plotSignatures(sample_1, sub = i)
    dev.off()
  })
}


