rm(list=ls())
options(stringsAsFactors = F)
getwd()
Rdata_dir='./Rdata/'
Figure_dir='./figures/'
if(!require(RTCGA))BiocManager::install("RTCGA")
if(!require(RTCGA.miRNASeq))BiocManager::install("RTCGA.miRNASeq")
if(!require(RTCGA.clinical))BiocManager::install("RTCGA.clinical")
if(F){
  library(RTCGA.miRNASeq) 
  ??RTCGA.miRNASeq
  ?KIRC.miRNASeq
  s=rownames(KIRC.miRNASeq)[seq(1,nrow(KIRC.miRNASeq),by=3)]
  expr <- expressionsTCGA(KIRC.miRNASeq)
  dim(expr)
  expr[1:40,1:4]
  expr=as.data.frame(expr[seq(1,nrow(expr),by=3),3:ncol(expr)])
  dim(expr)
  mi=colnames(expr)
  expr=apply(expr,1,as.numeric) 
  colnames(expr)=s
  rownames(expr)=mi
  expr[1:4,1:4]
  expr=na.omit(expr)
  #10个样本>1才留下
  dim(expr)
  expr=expr[apply(expr, 1,function(x){sum(x>1)>10}),]
  dim(expr)
  #临床信息
  library(RTCGA.clinical) 
  meta <- KIRC.clinical
  dim(meta)
  tmp=as.data.frame(colnames(meta))
  meta[(grepl('patient.bcr_patient_barcode',colnames(meta)))]
  meta[(grepl('patient.days_to_last_followup',colnames(meta)))]
  meta[(grepl('patient.days_to_death',colnames(meta)))]
  meta[(grepl('patient.vital_status',colnames(meta)))]

  meta=as.data.frame(meta[c('patient.bcr_patient_barcode','patient.vital_status',
                            'patient.days_to_death','patient.days_to_last_followup',
                            'patient.race',
                            'patient.age_at_initial_pathologic_diagnosis',
                            'patient.gender' ,
                           'patient.stage_event.pathologic_stage')])
  #meta[(grepl('patient.stage_event.pathologic_stage',colnames(meta)))]
  ## 每次运行代码，就会重新生成文件。
  save(expr,meta,
       file = file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
         )
}

load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

