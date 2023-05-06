library(tibble)

# 定义两个简单的函数
up=function(df){
  rownames(df)[df$change=="UP"]
}
down=function(df){
  rownames(df)[df$change=="DOWN"]
}

GSE = read.table('GSE59247.top.table-1.tsv',sep="\t",header=T,check.names=F)
rownames(GSE) <- GSE[,1]
GSE <- GSE[,-1]



logFC_cutoff=1
padj=0.05

# 添加显著性标签
GSE$change <- factor(ifelse(GSE$PValue < padj & abs(GSE$logFC) > logFC_cutoff,
                            ifelse(GSE$logFC > logFC_cutoff ,'UP','DOWN'),
                            'NOT'))

TCGA = read.table('miRNA-TCGA.csv',sep="\t",header=T,check.names=F)


# 查看上调及下调基因数目
daframe <- data.frame(gse=as.integer(table(GSE$change)),
                      tcg=as.integer(table(TCGA$change)),
                      row.names = c('DOWN','NOT','UP'))




# 韦恩图
# 上调三维维恩图
#install.packages('ggVennDiagram')
library(ggVennDiagram)
library(ggplot2)
mydata<-list(gse=up(GSE),tcga=up(TCGA))
ggVennDiagram(mydata)

# 换个配色
#install.packages('paletteer')
library(paletteer)
View(palettes_c_names)
ggVennDiagram(mydata) + 
  scale_fill_paletteer_c('ggthemes::Orange-Blue Light Diverging', direction = 1)

ggVennDiagram(mydata,
              category.names = c("GSE59247","TCGA-BRCA"),
              set_color="black",
              set_size=5,
              label = "both",
              label_percent_digit=0,
              label_color="black",
              label_size=5,
              label_alpha = 0)+
  scale_fill_paletteer_c('ggthemes::Orange-Blue Light Diverging', direction = 1)

up(GSE)[which(up(GSE) %in% up(TCGA))]
