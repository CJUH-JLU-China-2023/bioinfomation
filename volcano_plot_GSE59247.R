library(ggplot2)
rm(list = ls())
# 添加点标签
library(ggrepel)
logFC_cutoff=1
padj=0.05

DEG = read.table('GSE59247.top.table-1.tsv',sep="\t",header=T,check.names=F)
rownames(DEG) <- DEG[,1]
DEG <- DEG[,-1]

# 添加显著性标签
DEG$change <- factor(ifelse(DEG$PValue < padj & abs(DEG$logFC) > logFC_cutoff,
                            ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),
                            'NOT'))

significant = subset(DEG[c('hsa-miR-21-5p','hsa-miR-96-5p','hsa-miR-183-5p','hsa-miR-592','hsa-miR-3923','hsa-miR-885-5p','hsa-miR-142-5p','hsa-miR-339-3p','hsa-miR-184'),])



g <- ggplot(DEG,aes(logFC,-log10(as.numeric(PValue)),color=change))

g + geom_point()

g + geom_point(size = 0.5) + 
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'volcano_plot') + 
  theme(plot.title = element_text(hjust = 0.5, size = 20,color = 'firebrick'), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2, color = 'black') +  #添加阈值线
  geom_hline(yintercept = -log10(padj), lty = 2, color = 'black') 


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
                  size=3,col = 'black',box.padding = unit(2, 'lines'),segment.color = 'black', show.legend = FALSE)


ggsave('volcano_plot_GSE59247.eps',height = 20,width = 30,units = 'cm')
?ggsave
