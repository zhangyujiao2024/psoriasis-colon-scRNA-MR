#仅运行聚类步骤
# 核心单细胞分析包
library(Seurat)        # 单细胞数据分析核心包
library(dplyr)         # 数据处理和管道操作
library(magrittr)      # 管道操作符
library(tidyr)         # 数据整理
library(stringr)       # 字符串处理
library(ggplot2)       # 基础绘图
library(ggpubr)        # 统计绘图和发表级图形
library(RColorBrewer)  # 调色板
library(viridis)       # 现代配色方案
library(scales)        # 图形标度
library(patchwork)     # 图形组合

message("所有必需的R包已成功加载完成！")

# ==============================================================================
# 02. 参数设置和工作目录
# ==============================================================================

# 设置工作目录
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\32.聚类"

# 创建结果输出目录
output_dirs <- list(
  basic_analysis = "01.Basic_Analysis_Results",
  umap_plots = "02.UMAP_Clustering_Plots"
)

# 创建所有输出目录
for (dir_name in output_dirs) {
  full_path <- file.path(workDir, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    message(sprintf("创建目录: %s", dir_name))
  }
}

# 检查并设置工作目录
if (!file.exists(workDir)) {
  stop("工作目录不存在，请检查路径：", workDir)
} else {
  setwd(workDir)
  message("工作目录设置为：", getwd())
}

# 分析参数设置
analysis_params <- list(
  pcSelect = 12,                      # PCA主成分数选择
  cluster_resolution = 0.3,           # 聚类分辨率
  n_variable_features = 2000          # 高变基因数量
)

# ==============================================================================
# 03. 数据读取和Seurat对象创建
# ==============================================================================

message("正在读取单细胞数据...")

# 读取数据文件
single_cell_data <- readRDS("Extracted_Cells.rds")

# 数据质量检查
if (ncol(single_cell_data) == 0) {
  stop("数据文件为空，请检查文件内容！")
}

message(sprintf("数据包含 %d 个基因和 %d 个细胞",
                nrow(single_cell_data), ncol(single_cell_data)))

# 创建Seurat对象
Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
  counts = single_cell_data,
  min.cells = 3,
  min.features = 200
)

message("Seurat对象创建完成")

# ==============================================================================
# 04. 数据预处理和聚类
# ==============================================================================

message("开始数据预处理和聚类...")

# 数据标准化
Peripheral_Blood_Mononuclear_Cells <- NormalizeData(
  object = Peripheral_Blood_Mononuclear_Cells,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
message("1/7 数据标准化完成")

# 高变基因识别
Peripheral_Blood_Mononuclear_Cells <- FindVariableFeatures(
  object = Peripheral_Blood_Mononuclear_Cells,
  selection.method = "vst",
  nfeatures = analysis_params$n_variable_features
)
message("2/7 高变基因识别完成")

# 数据缩放
Peripheral_Blood_Mononuclear_Cells <- ScaleData(Peripheral_Blood_Mononuclear_Cells)
message("3/7 数据标准化完成")

# PCA分析
Peripheral_Blood_Mononuclear_Cells <- RunPCA(
  object = Peripheral_Blood_Mononuclear_Cells,
  npcs = 20,
  features = VariableFeatures(object = Peripheral_Blood_Mononuclear_Cells)
)
message("4/7 PCA分析完成")

# 邻居查找
Peripheral_Blood_Mononuclear_Cells <- FindNeighbors(
  object = Peripheral_Blood_Mononuclear_Cells,
  dims = 1:analysis_params$pcSelect
)
message("5/7 邻居查找完成")

# 聚类分析
Peripheral_Blood_Mononuclear_Cells <- FindClusters(
  object = Peripheral_Blood_Mononuclear_Cells,
  resolution = analysis_params$cluster_resolution
)
message("6/7 聚类分析完成")

# UMAP降维
Peripheral_Blood_Mononuclear_Cells <- RunUMAP(
  object = Peripheral_Blood_Mononuclear_Cells,
  dims = 1:analysis_params$pcSelect
)
message("7/7 UMAP降维完成")

message("数据预处理和聚类完成！")

# ==============================================================================
# 05. 基础可视化
# ==============================================================================

message("生成基础可视化图形...")

# UMAP聚类图
pdf(file = file.path(output_dirs$umap_plots, "UMAP_Clusters.pdf"), width = 10, height = 8)
print(
  DimPlot(
    object = Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    label = TRUE,
    pt.size = 0.5
  ) +
    ggtitle("UMAP Clustering Results") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
)
dev.off()
message("已保存聚类UMAP图: UMAP_Clusters.pdf")

# 分组UMAP图 - 提取对照组和疾病组信息
# 从细胞名称中提取分组信息 (格式: Control.xxx 或 Treat.xxx)
Type <- gsub("(.*?)\\..*", "\\1", colnames(Peripheral_Blood_Mononuclear_Cells))
names(Type) <- colnames(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
  object = Peripheral_Blood_Mononuclear_Cells,
  metadata = Type,
  col.name = "orig.ident"
)

message(sprintf("检测到分组: %s", paste(unique(Type), collapse = ", ")))

# 生成分组UMAP图 - 对照组vs疾病组
pdf(file = file.path(output_dirs$umap_plots, "UMAP_Clusters_Split_by_Group.pdf"), width = 14, height = 8)
print(
  DimPlot(
    object = Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    split.by = "orig.ident",
    label = TRUE,
    ncol = 2,
    pt.size = 0.5
  ) +
    ggtitle("UMAP Clusters: Control vs Treatment") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
)
dev.off()
message("已保存分组UMAP图: UMAP_Clusters_Split_by_Group.pdf")

# 保存聚类信息
write.table(
  Peripheral_Blood_Mononuclear_Cells$seurat_clusters,
  file = file.path(output_dirs$basic_analysis, "umapCluster.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE
)

message("聚类信息已保存")

# 保存Seurat对象
saveRDS(
  Peripheral_Blood_Mononuclear_Cells,
  file = file.path(output_dirs$basic_analysis, "seurat_clustering_result.rds")
)

message("Seurat对象已保存")

# ==============================================================================
# 完成
# ==============================================================================

message("\n========================================")
message("聚类分析完成！")
message("========================================")
message(sprintf("检测到 %d 个聚类", length(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))))
message(sprintf("细胞总数: %d", ncol(Peripheral_Blood_Mononuclear_Cells)))
message(sprintf("基因总数: %d", nrow(Peripheral_Blood_Mononuclear_Cells)))
message("\n生成的文件:")
message(sprintf("  - UMAP图: %s", file.path(output_dirs$umap_plots, "UMAP_Clusters.pdf")))
message(sprintf("  - 分组UMAP图: %s", file.path(output_dirs$umap_plots, "UMAP_Clusters_Split_by_Group.pdf")))
message(sprintf("  - 聚类信息: %s", file.path(output_dirs$basic_analysis, "umapCluster.txt")))
message(sprintf("  - Seurat对象: %s", file.path(output_dirs$basic_analysis, "seurat_clustering_result.rds")))
