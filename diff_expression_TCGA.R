library(pacman)
p_load(TCGAbiolinks)
p_load(stringi)
rm(list = ls())
options(stringsAsFactors = F)

#--------------------------数据准备和预处理-------------------------------------
#可以通过TCGAbiolinks:::getGDCprojects()$project_id得到每个项目的ID
TCGAbiolinks:::getGDCprojects()$project_id
project_id <- "TCGA-BRCA"


#可以通过TCGAbiolinks:::getProjectSummary(project_id)得到可以下载的数据类型
TCGAbiolinks:::getProjectSummary(project_id)
data_category <- "Transcriptome Profiling"

#data.type
data_type <- "Isoform Expression Quantification"

#查询数据
query <- GDCquery(project = project_id,
                  data.category = data_category,
                  data.type = data_type)
#下载
GDCdownload(query)

#整合，这里得到的dataAssy是一个数据框，并且每一个样本的每一个miRNA一行，后面需要处理一下
dataAssy <- GDCprepare(query = query,
                       summarizedExperiment=F)


sample_id <- as.data.frame(table(dataAssy$barcode))

miRNA_id <- as.data.frame(table(dataAssy$miRNA_region))

miRNA_matrue_RPM <- matrix(NA,ncol = nrow(sample_id),nrow = nrow(miRNA_id))
colnames(miRNA_matrue_RPM) <- sample_id$Var1
rownames(miRNA_matrue_RPM) <- as.character(miRNA_id$Var1)

for(i in 1:nrow(sample_id)){
  temp1 <- dataAssy[which(dataAssy$barcode==as.character(sample_id[i,1])),]
  for(j in 1:nrow(miRNA_id)){
    loc <- which(temp1$miRNA_region==as.character(miRNA_id[j,1]))
    if(length(loc)>0){
      miRNA_matrue_RPM[j,i] <- sum(temp1[loc,4])
    }else{
      miRNA_matrue_RPM[j,i] <- 0
    }
  }
  print(i)
}
dim(miRNA_matrue_RPM)
miRNA_matrue_RPM1 <- miRNA_matrue_RPM[1:2236,]


name <- substr(rownames(miRNA_matrue_RPM1),8,nchar(rownames(miRNA_matrue_RPM1)))


p_load(miRBaseVersions.db)
items <- select(miRBaseVersions.db,
                keys = name,
                keytype = "MIMAT",
                columns = c("ACCESSION","NAME","VERSION"))

id_name <- items[items$VERSION == 21.0, c("ACCESSION","NAME")]
miRNA_matrue_RPM2 <- cbind(id_name,miRNA_matrue_RPM1)

save(miRNA_matrue_RPM2,file = paste(project_id,"_miRNA_matrue_RPM.RData",sep=""))

write.table(miRNA_matrue_RPM2,file = paste(project_id,"_miRNA_matrue_RPM.xls",sep=""),sep="\t",row.names = F,quote=F)

# 导入数据
rt <- read.table("TCGA-BRCA_miRNA_matrue_RPM.xls",sep="\t",header=T,check.names=F)

# 数据处理
rownames(rt) <- rt[,2]
rt <- rt[,-1]
rt <- rt[,-1]
dim(rt)
write.table(rt, file="rt.xls", sep="\t", quote=F, col.names=T)

# 用substr函数在TCGA数据中提取样本信息
group <- factor(ifelse(as.integer(substr(colnames(rt),14,15))<10,'tumor','normal'),
                levels = c('normal','tumor'))
table(group)


#-----------------------------差异表达分析--------------------------------------------

# 差异性分析edgeR
library(edgeR)

logFC_cutoff=1
padj=0.05

dge <- DGEList(counts=rt,group=group)   # 将数据打包为edgeR包识别的数据
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 
design <- model.matrix(~0+factor(group))
colnames(design) <- c('normal','tumor')
rownames(design)<-colnames(dge)
colnames(design)<-levels(group)



dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit1 <- glmLRT(fit, contrast=c(-1,1)) 
DEG=topTags(fit1, n=nrow(rt))
DEG=as.data.frame(DEG)

write.csv(design, file="design.xls")


# 输出有显著差异的结果
significant = DEG[(
  (DEG$FDR < padj) & (abs(DEG$logFC) > logFC_cutoff)
),]  
write.table(significant, file="significant.xls",sep="\t",quote=F)

# 输出表达上调的结果
upregulation = DEG[(DEG$FDR < padj & (DEG$logFC>logFC_cutoff)),]
write.table(upregulation, file="upregulation.xls",sep="\t",quote=F)

# 输出表达下调的结果
downregulation = DEG[(DEG$FDR < padj & (DEG$logFC<(-logFC_cutoff))),]
write.table(downregulation, file="downregulation.xls",sep="\t",quote=F)


# 添加显著性标签
DEG$change <- factor(ifelse(DEG$PValue < padj & abs(DEG$logFC) > logFC_cutoff,
                            ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),
                            'NOT'))
up = DEG$change[DEG$change == 'UP']
down = DEG$change[DEG$change == 'DOWN']

# 输出结果
write.table(DEG,file="miRNA-edgeR.csv",sep="\t")


#-----------------------------------火山图绘制---------------------------------- 
# ggplot2绘制火山图 volcano
library(ggplot2)
g <- ggplot(DEG,aes(logFC,-log10(as.numeric(PValue)),color=change))

g + geom_point()

g + geom_point(size = 0.5) + 
  labs(x = 'log2 Fold Change', y = '-log10 p-value', title = 'volcano_plot') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20,color = 'firebrick'), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2, color = 'black') +  #添加阈值线
  geom_hline(yintercept = -log10(padj), lty = 2, color = 'black') 

# 添加点标签
library(ggrepel)

# 将感兴趣的miRNA在图中标记出来
significant = subset(DEG[c('hsa-miR-21-5p','hsa-miR-96-5p','hsa-miR-183-5p','hsa-miR-592','hsa-miR-3923','hsa-miR-885-5p','hsa-miR-142-5p','hsa-miR-339-3p','hsa-miR-184'),])

g + geom_point(size = 1) + 
  labs(x = 'log2 Fold Change', y = '-log10 p-value') + 
  theme(plot.title = element_text(hjust = 0.5, size = 15,color = 'firebrick'), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) + 
  scale_color_manual(values =c('UP' = 'red', 'DOWN' = 'blue', 'NOT' = 'gray'))+
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2, color = 'black') +  #添加阈值线
  geom_hline(yintercept = -log10(padj), lty = 2, color = 'black') +
  geom_text_repel(data = significant,aes(label = rownames(significant)),
                  max.overlaps = 10000, # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖
                  size=3,col = 'black',box.padding = unit(1.5, 'lines'),segment.color = 'black', show.legend = FALSE)

ggsave('volcano_plot.eps',height = 20,width = 30,units = 'cm')