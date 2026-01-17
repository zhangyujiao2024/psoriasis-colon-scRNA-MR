###################################01.数据前处理和分析###################################

# 导入必要的库
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(ggplot2)      #
library(RColorBrewer) # 用于生成调色板

#========================== 参数设置区域 ==========================
# 工作目录设置
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\04.单细胞质量控制、聚类和PCA选择"  # 设置工作目录

# 质量控制阈值设置
min_cells_per_feature <- 3              # 每个基因至少在多少个细胞中表达
min_features_per_cell <- 200            # 每个细胞至少表达多少个基因
max_mt_percent <- 15                    # 线粒体基因百分比上限
max_features_per_cell <- 5000           # 每个细胞最大基因数（可选，用于去除doublets）

# 归一化参数
normalization_method <- "LogNormalize"   # 归一化方法
scale_factor <- 10000                   # 归一化因子
n_variable_features <- 2000              # 选择的高变基因数量
variable_selection_method <- "vst"       # 高变基因选择方法

# PCA分析参数
n_pcs <- 20                             # PCA主成分数量
n_pcs_for_jackstraw <- 20               # JackStraw分析的主成分数量
jackstraw_replicates <- 100             # JackStraw重采样次数
cumulative_variance_threshold <- 0.90   # 累积方差贡献率阈值

# 可视化参数
top_genes_to_show <- 10                 # 显示的top基因数量
heatmap_cells <- 500                    # 热图显示的细胞数量
heatmap_features <- 30                  # 热图显示的基因数量
plot_width <- 10                        # 图片宽度
plot_height <- 8                        # 图片高度

# 图形标签和样式参数
axis_text_size <- 12                    # 坐标轴文字大小
axis_title_size <- 14                   # 坐标轴标题大小
plot_title_size <- 16                   # 图形标题大小
legend_text_size <- 10                  # 图例文字大小
legend_title_size <- 12                 # 图例标题大小
x_axis_angle <- 45                      # X轴标签倾斜角度

# PCA图形标签设置
pca_loading_features_to_show <- 20       # PCA负载图显示的基因数量

#================================================================

setwd(workDir)

# 读取单细胞数据文件 (single_cell_expression_matrix.rds)
single_cell_data <- readRDS("single_cell_expression_matrix.rds")

# 将数据转换为 Seurat 对象
Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
  counts = single_cell_data,
  min.cells = min_cells_per_feature,
  min.features = min_features_per_cell
)

# 使用 PercentageFeatureSet 函数计算线粒体基因的百分比
Peripheral_Blood_Mononuclear_Cells[["percent.mt"]] <- PercentageFeatureSet(
  object = Peripheral_Blood_Mononuclear_Cells,
  pattern = "^MT-"
)

# 保存过滤前的对象用于比较
before_filtering <- Peripheral_Blood_Mononuclear_Cells

# 导入patchwork包用于图表合并
library(patchwork)

# 绘制基因特征的小提琴图（过滤前）
pdf(file = "featureViolin.pdf", width = 15, height = 6.5)

# 创建过滤前的小提琴图
violin_before <- VlnPlot(object = before_filtering,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         pt.size = 0) +  # 去掉散点
  ggtitle("Before Filtering") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(violin_before)
dev.off()

# 对数据进行过滤
Peripheral_Blood_Mononuclear_Cells <- subset(
  x = Peripheral_Blood_Mononuclear_Cells,
  subset = nFeature_RNA > min_features_per_cell &
    nFeature_RNA < max_features_per_cell &
    percent.mt < max_mt_percent
)

# ================ 查看过滤后统计信息 ================
# 查看过滤后的细胞数量
filtered_cells <- ncol(Peripheral_Blood_Mononuclear_Cells)
message(paste("过滤后细胞数量:", filtered_cells))

# 查看过滤后的基因数量
filtered_genes <- nrow(Peripheral_Blood_Mononuclear_Cells)
message(paste("过滤后基因数量:", filtered_genes))

# 查看基因表达概况
genes_expressed <- rowSums(GetAssayData(Peripheral_Blood_Mononuclear_Cells, slot = "counts") > 0)
summary_stats <- summary(genes_expressed)
message("基因表达概况:")
print(summary_stats)

# 保存统计信息到文件
filtered_stats <- data.frame(
  指标 = c("细胞数量", "基因数量", 
         "平均每个基因表达的细胞数", 
         "中位数每个基因表达的细胞数"),
  数值 = c(filtered_cells, 
         filtered_genes,
         mean(genes_expressed),
         median(genes_expressed))
)
write.csv(filtered_stats, "filtered_cell_gene_stats.csv", row.names = FALSE)

# ====================================================

# 绘制基因特征的小提琴图（过滤后）
pdf(file = "featureViolin_after.pdf", width = 15, height = 6.5)

# 创建过滤后的小提琴图
violin_after <- VlnPlot(object = Peripheral_Blood_Mononuclear_Cells,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        pt.size = 0) +  # 去掉散点
  ggtitle("After Filtering") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(violin_after)
dev.off()

# 创建对比图（过滤前后）
pdf(file = "featureViolin_comparison.pdf", width = 20, height = 12)

# 重新创建小提琴图用于对比
violin_before_comp <- VlnPlot(object = before_filtering,
                              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                              pt.size = 0) +
  ggtitle("Before Filtering") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

violin_after_comp <- VlnPlot(object = Peripheral_Blood_Mononuclear_Cells,
                             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                             pt.size = 0) +
  ggtitle("After Filtering") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# 合并两个图（上下排列）
comparison_plot <- violin_before_comp / violin_after_comp

print(comparison_plot)
dev.off()

# 提取过滤后的原始表达矩阵
filtered_expression_matrix <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, assay = "RNA", slot = "counts")

# 保存原始表达矩阵为RDS文件 (用于后续分析脚本)
saveRDS(filtered_expression_matrix, file = "filtered_single_cell_data.rds")
message("过滤后的表达矩阵已保存为 'filtered_single_cell_data.rds'")

# 可选：同时保存完整的Seurat对象 (用于其他分析)
saveRDS(Peripheral_Blood_Mononuclear_Cells, file = "filtered_seurat_object.rds")
message("过滤后的完整Seurat对象已保存为 'filtered_seurat_object.rds'")

# 进行数据归一化和特征选择
Peripheral_Blood_Mononuclear_Cells <- NormalizeData(
  object = Peripheral_Blood_Mononuclear_Cells,
  normalization.method = normalization_method,
  scale.factor = scale_factor
)

Peripheral_Blood_Mononuclear_Cells <- FindVariableFeatures(
  object = Peripheral_Blood_Mononuclear_Cells,
  selection.method = variable_selection_method,
  nfeatures = n_variable_features
)

# 查看前N个变量基因
top_genes <- head(x = VariableFeatures(object = Peripheral_Blood_Mononuclear_Cells), top_genes_to_show)

# 绘制特征选择图
pdf(file = "featureVar.pdf", width = 8, height = 12)
plot1 <- VariableFeaturePlot(object = Peripheral_Blood_Mononuclear_Cells)
plot2 <- LabelPoints(plot = plot1, points = top_genes, repel = TRUE)
CombinePlots(plots = list(plot2))
dev.off()

summary(Peripheral_Blood_Mononuclear_Cells)

nFeature_RNA_summary <- summary(Peripheral_Blood_Mononuclear_Cells[["nFeature_RNA"]])
write.csv(nFeature_RNA_summary, "nFeature_RNA_summary.csv")

nCount_RNA_summary <- summary(Peripheral_Blood_Mononuclear_Cells[["nCount_RNA"]])
write.csv(nCount_RNA_summary, "nCount_RNA_summary.csv")

percent_mt_summary <- summary(Peripheral_Blood_Mononuclear_Cells[["percent.mt"]])
write.csv(percent_mt_summary, "percent_mt_summary.csv")
# 获取前N个变量基因
top_genes_list <- head(VariableFeatures(Peripheral_Blood_Mononuclear_Cells), top_genes_to_show)

# 查看前N个变量基因的内容
print(top_genes_list)

# 如果top_genes_list是基因名称而非索引，可以将其转换为数据框并保存
top_genes_df <- data.frame(Genes = top_genes_list)

# 将前N个变量基因保存为CSV文件
write.csv(top_genes_df, paste0("top", top_genes_to_show, "_variable_genes.csv"), row.names = FALSE)





# PCA 分析：标准化数据
Peripheral_Blood_Mononuclear_Cells <- ScaleData(Peripheral_Blood_Mononuclear_Cells)  # 对数据进行标准化，使每个基因的均值为 0，标准差为 1

# 运行PCA分析：检查是否有变异基因
if (length(VariableFeatures(Peripheral_Blood_Mononuclear_Cells)) == 0) {
  stop("没有选择到任何高变基因，请检查数据预处理！")
}

Peripheral_Blood_Mononuclear_Cells <- RunPCA(
  object = Peripheral_Blood_Mononuclear_Cells,
  npcs = n_pcs,
  features = VariableFeatures(object = Peripheral_Blood_Mononuclear_Cells)
)

# 可视化 PCA 的负载图
for (i in 1:4) {
  # 打印当前处理的主成分
  message(paste("正在生成PCA负载图，主成分：", i))

  # 生成图形对象
  pca_plot <- VizDimLoadings(
    object = Peripheral_Blood_Mononuclear_Cells,
    dims = i,
    reduction = "pca",
    nfeatures = pca_loading_features_to_show
  )

  # 判断图形是否生成成功
  if (is.null(pca_plot)) {
    warning(paste("PCA负载图生成失败，主成分：", i))
    next  # 跳过当前循环
  }

  # 修改样式，使用 ggplot2 自定义美化图形并添加坐标轴标签
  pca_plot <- pca_plot +
    theme_minimal() +  # 使用简洁的主题
    theme(
      plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5),  # 设置标题的样式
      axis.text = element_text(size = axis_text_size, color = "black"),  # 设置坐标轴标签的字体大小和颜色
      axis.title = element_text(size = axis_title_size, face = "bold"),  # 设置坐标轴标题的字体大小和粗体
      panel.grid.major = element_line(color = "grey90", size = 0.5),  # 设置网格线样式
      panel.grid.minor = element_blank(),  # 隐藏次网格线
      axis.text.x = element_text(angle = x_axis_angle, hjust = 1, size = axis_text_size),  # X轴标签倾斜显示
      axis.text.y = element_text(size = axis_text_size)  # Y轴标签
    ) +
    ggtitle(paste("PCA Loadings Plot - Principal Component", i)) +  # 标题
    xlab(paste("Genes (Top", pca_loading_features_to_show, "features for PC", i, ")")) +  # X轴标签
    ylab(paste("Loading Scores for PC", i)) +  # Y轴标签
    scale_color_brewer(palette = "Set2")  # 使用漂亮的颜色调色板

  # 保存PDF文件
  pdf(file = paste0("pcaGene_", i, ".pdf"), width = plot_width, height = plot_height)

  # 重新绘制图形到PDF文件
  print(pca_plot)  # 使用 print() 来确保图形被正确输出到PDF

  # 关闭PDF输出
  dev.off()

  # 打印完成信息
  message(paste("主成分", i, "的PCA负载图已保存"))
}

# 可视化 PCA 样本分布图
pdf(file = "PCA.pdf", width = 6.5, height = 6)
# 将PCA样本分布图保存为PDF文件

# 创建PCA样本分布图并添加坐标轴标签
pca_sample_plot <- DimPlot(object = Peripheral_Blood_Mononuclear_Cells, reduction = "pca") +
  NoGrid() +  # 移除网格线但保留坐标轴线
  theme(
    plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5),
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.title = element_text(size = axis_title_size, face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),  # 显示坐标轴线
    legend.text = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_title_size, face = "bold")
  ) +
  ggtitle("PCA Sample Distribution") +
  xlab("PC1 (Principal Component 1)") +  # X轴标签
  ylab("PC2 (Principal Component 2)")    # Y轴标签

print(pca_sample_plot)  # 绘制PCA样本分布图
dev.off()  # 关闭文件输出

# 可视化 PCA 热图
for (i in 1:4) {
  # 打印当前处理的主成分
  message(paste("正在生成PCA热图，主成分：", i))

  # 保存PCA热图为单独的PDF文件
  pdf(file = paste0("pcaHeatmap_", i, ".pdf"), width = plot_width, height = plot_height)

  # 绘制PCA热图，展示第 i 个主成分的细胞表达模式
  DimHeatmap(
    object = Peripheral_Blood_Mononuclear_Cells,
    dims = i,
    cells = heatmap_cells,
    balanced = TRUE,
    nfeatures = heatmap_features,
    ncol = 1
  )

  # 关闭PDF输出
  dev.off()

  # 打印完成信息
  message(paste("主成分", i, "的PCA热图已保存"))
}

# JackStraw 分析：评估每个主成分的重要性
Peripheral_Blood_Mononuclear_Cells <- JackStraw(
  object = Peripheral_Blood_Mononuclear_Cells,
  num.replicate = jackstraw_replicates
)

Peripheral_Blood_Mononuclear_Cells <- ScoreJackStraw(
  object = Peripheral_Blood_Mononuclear_Cells,
  dims = 1:n_pcs_for_jackstraw
)

# 可视化 JackStraw 分析结果
pdf(file = "pcaJackStraw.pdf", width = 9, height = 6)
JackStrawPlot(
  object = Peripheral_Blood_Mononuclear_Cells,
  dims = 1:15
)
dev.off()

# 获取PCA结果
pca_results <- Peripheral_Blood_Mononuclear_Cells[["pca"]]@stdev
pca_variance <- pca_results^2 / sum(pca_results^2)  # 方差贡献率
cumulative_variance <- cumsum(pca_variance)  # 累积方差贡献率

# 创建PCA结果表格
pca_table <- data.frame(PC = 1:length(pca_variance), 
                        StdDev = pca_results, 
                        Variance = pca_variance, 
                        CumulativeVariance = cumulative_variance)

# 显示PCA结果表
DT::datatable(pca_table)

# 将PCA结果表导出为CSV文件
write.csv(pca_table, file = "PCA_Results.csv", row.names = FALSE)

# 打印导出成功的信息
message("PCA results have been exported as 'PCA_Results.csv'")


# 访问 JackStraw 结果并提取 p 值
jackstraw_results <- Peripheral_Blood_Mononuclear_Cells[["pca"]]@jackstraw

# 提取每个主成分的整体 p 值
overall_p_values <- jackstraw_results$overall.p.values

# 打印每个主成分的 p 值
print(overall_p_values)

# 创建一个数据框显示每个主成分的p值和得分
overall_p_values_table <- data.frame(
  PC = 1:nrow(overall_p_values),
  PValue = overall_p_values[, 1],  # p值列
  Score = overall_p_values[, 2]    # 得分列
)

# 打印 p 值表格
print(overall_p_values_table)

# 导出 p 值表格为 CSV 文件
write.csv(overall_p_values_table, file = "JackStraw_Overall_PValues.csv", row.names = FALSE)

# 假设你已经计算了 pca_variance 和 cumulative_variance

# 打印各主成分的方差贡献率
print(pca_variance)

# 选择累积方差贡献率达到指定阈值的主成分数
selected_pcs <- which(cumulative_variance >= cumulative_variance_threshold)[1]
print(paste("选择的主成分数量：", selected_pcs, "（累积方差贡献率阈值：", cumulative_variance_threshold * 100, "%）"))

# 方法1：使用ElbowPlot（Seurat自带）
pdf(file = "pca_elbow_plot.pdf", width = 8, height = 6)
elbow_plot <- ElbowPlot(
  object = Peripheral_Blood_Mononuclear_Cells,
  ndims = n_pcs  # 显示所有主成分
) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5),
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.title = element_text(size = axis_title_size, face = "bold"),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("PCA Elbow Plot") +
  xlab("Principal Components") +
  ylab("Standard Deviation")

# 在图中添加累积方差贡献率阈值线（可选）
if (exists("cumulative_variance")) {
  # 找到达到阈值的主成分
  threshold_pc <- which(cumulative_variance >= cumulative_variance_threshold)[1]
  if (!is.na(threshold_pc)) {
    elbow_plot <- elbow_plot + 
      geom_vline(xintercept = threshold_pc, 
                 linetype = "dashed", 
                 color = "red", 
                 size = 1) +
      annotate("text", 
               x = threshold_pc + 1, 
               y = max(pca_results) * 0.9,
               label = paste("PC", threshold_pc, 
                             "\n", 
                             round(cumulative_variance[threshold_pc] * 100, 1), 
                             "% variance"),
               color = "red",
               size = 4)
  }
}

print(elbow_plot)
dev.off()






# 创建表格
pca_table <- data.frame(
  Principal_Component = 1:length(pca_variance),
  Variance = pca_variance,
  Cumulative_Variance = cumulative_variance
)

# 将结果导出为CSV文件
write.csv(pca_table, file = "PCA_Variance_Results.csv", row.names = FALSE)

# 还可以将选定的主成分数量也导出
selected_pcs_table <- data.frame(Selected_Principal_Components = selected_pcs)
write.csv(selected_pcs_table, file = "Selected_PCs.csv", row.names = FALSE)

# 输出文件路径
print("CSV文件已保存为 'PCA_Variance_Results.csv' 和 'Selected_PCs.csv'")

#

############################# 结束 ##################################
