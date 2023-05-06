#清空环境变量
rm(list = ls())

#加载包
library(pROC)
getwd()
#读取miRNA表达数据
miRNA <- read.csv("mirna_data.csv",sep="\t", header = TRUE)

mirna_data <- t(miRNA)

# 获取行名
row_names <- row.names(mirna_data)

# 根据行名的第14-15位是否小于10添加标签
labels <- ifelse(as.numeric(substr(row_names, 14, 15)) < 10, "Positive", "Negative")



# 选择需要绘制ROC曲线的miRNA
mirna_name_21 <- "hsa-miR-21-5p"
mirna_name_183 <- "hsa-miR-183-5p"
mirna_name_96 <- "hsa-miR-96-5p"
mirna_name_592 <- "hsa-miR-592"
mirna_name_3923 <- "hsa-miR-3923"
mirna_name_885 <- "hsa-miR-885-5p"
mirna_name_184 <- "hsa-miR-184"
mirna_name_339 <- "hsa-miR-339-3p"
mirna_name_142 <- "hsa-miR-142-5p"
# 获取表达数据和阳性/阴性标签数据
mirna_expression_21 <- mirna_data[, mirna_name_21]
mirna_expression_183 <- mirna_data[, mirna_name_183]
mirna_expression_96 <- mirna_data[, mirna_name_96]
mirna_expression_592 <- mirna_data[, mirna_name_592]
mirna_expression_3923 <- mirna_data[, mirna_name_3923]
mirna_expression_885 <- mirna_data[, mirna_name_885]
mirna_expression_184 <- mirna_data[, mirna_name_184]
mirna_expression_339 <- mirna_data[, mirna_name_339]
mirna_expression_142 <- mirna_data[, mirna_name_142]
#创建ROC对象
roc_data_21 <- roc(labels, mirna_expression_21)
roc_data_183 <- roc(labels, mirna_expression_183)
roc_data_96 <- roc(labels, mirna_expression_96)
roc_data_592 <- roc(labels, mirna_expression_592)
roc_data_3923 <- roc(labels, mirna_expression_3923)
roc_data_885 <- roc(labels, mirna_expression_885)
roc_data_184 <- roc(labels, mirna_expression_184)
roc_data_339 <- roc(labels, mirna_expression_339)
roc_data_142 <- roc(labels, mirna_expression_142)
pdf('roc_miR.pdf')
plot(roc_data_21,legacy.axes = TRUE,
     col="red",identity.col="skyblue2",main = paste0("ROC curve for TCGA-BRCA"))
plot(roc_data_183,add=T,legacy.axes = TRUE,
     col="green",identity.col="skyblue2")
plot(roc_data_96,add=T,legacy.axes = TRUE,
     col="blue",identity.col="skyblue2")
plot(roc_data_592,add=T,legacy.axes = TRUE,
     col="orange",identity.col="skyblue2")
plot(roc_data_3923,add=T,legacy.axes = TRUE,
     col="gray",identity.col="skyblue2")
plot(roc_data_885,add=T,legacy.axes = TRUE,
     col="cyan",identity.col="skyblue2")
plot(roc_data_184,add=T,legacy.axes = TRUE,
     col="yellow",identity.col="skyblue2")
plot(roc_data_339,add=T,legacy.axes = TRUE,
     col="brown",identity.col="skyblue2")
plot(roc_data_142,add=T,legacy.axes = TRUE,
     col="black",identity.col="skyblue2")
legend("bottomright",
       legend = c("miR-21-5p","miR-183-5p","miR-96-5p","miR-592","miR-3923","miR-885-5p","miR-184","miR-339-3p","miR-142-5p"),
       col = c("red","green","blue","orange","gray","cyan","yellow","brown","black"),lty = 1)
dev.off()







