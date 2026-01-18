# ==============================================================================
# 单细胞轨迹分析 - 专用脚本
# 运行轨迹分析（monocle2）部分，使用已有的聚类结果
# ==============================================================================
# 移除当前版本
# 卸载当前版本igraph
#remove.packages("monocle")

# 加载必需的R包
library(Seurat)          # 单细胞分析核心包
library(dplyr)           # 数据处理
library(monocle)         # 轨迹分析（monocle2）
library(ggplot2)         # 画图
library(pheatmap)        # 热图绘制

message("✓ 所有必需的R包已成功加载！")

# ==============================================================================
# 参数设置和工作目录
# ==============================================================================
genes_to_analyze <- c("ALDH2", "ARL4C","LY9", "TYMP","FCER1G", "LYZ")
# 设置工作目录
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\33.细胞轨迹的分析"
setwd(workDir)
message("工作目录设置为：", getwd())


# ==============================================================================
# 数据读取
# ==============================================================================

message("\n开始读取数据...")

# 读取表达矩阵
message("  读取单细胞表达数据: Extracted_Cells.rds")
single_cell_data <- readRDS("Extracted_Cells.rds")

# 读取聚类信息
message("  读取聚类信息: umapCluster.txt")
cluster_info <- read.table("umapCluster.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")

# umapCluster.txt 格式：第一列是细胞ID，第二列是聚类号
# 创建聚类数据框
seurat_clusters <- data.frame(
  cell_id = cluster_info[, 1],
  seurat_clusters = cluster_info[, 2]
)
rownames(seurat_clusters) <- seurat_clusters$cell_id

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = single_cell_data,
  min.cells = 3,
  min.features = 200
)

# 数据标准化
seurat_obj <- NormalizeData(
  object = seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# 匹配聚类信息到Seurat对象
# 确保细胞顺序一致
cell_names <- colnames(seurat_obj)
matched_clusters <- seurat_clusters[cell_names, "seurat_clusters"]

# 添加聚类信息到元数据
seurat_obj@meta.data$seurat_clusters <- matched_clusters

# 提取分组信息（对照组vs实验组）
Type <- gsub("(.*?)\\..*", "\\1", colnames(seurat_obj))
names(Type) <- colnames(seurat_obj)
seurat_obj <- AddMetaData(
  object = seurat_obj,
  metadata = Type,
  col.name = "orig.ident"
)

message(sprintf("  数据包含 %d 个基因和 %d 个细胞", nrow(seurat_obj), ncol(seurat_obj)))
message(sprintf("  检测到的分组: %s", paste(unique(Type), collapse = ", ")))
message(sprintf("  检测到的聚类: %s", paste(sort(unique(seurat_obj$seurat_clusters)), collapse = ", ")))

# ==============================================================================
# 轨迹分析 - Monocle2
# ==============================================================================

message("\n开始轨迹分析...")

# ----------- 1. 准备monocle2所需的数据 -----------
message("  准备monocle2分析数据...")

# 提取标准化表达矩阵
monocle_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data"))

# 准备细胞元数据
monocle_sample <- seurat_obj@meta.data

# 准备基因注释数据
monocle_geneAnn <- data.frame(
  gene_short_name = rownames(monocle_matrix),
  row.names = rownames(monocle_matrix)
)

# 创建细胞类型注释（基于聚类）
# 简单的聚类标记标识为细胞类型
monocle_sample$cell_type2 <- paste0("Cluster_", monocle_sample$seurat_clusters)

# 获取高变基因用于轨迹分析
message("  识别高变基因...")
seurat_obj <- FindVariableFeatures(
  object = seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

trajectory_genes <- VariableFeatures(seurat_obj)
message(sprintf("  使用 %d 个高变基因进行轨迹分析", length(trajectory_genes)))

# ----------- 2. 创建monocle对象 -----------
message("  创建monocle CDS对象...")

# 将数据转换为稀疏矩阵
data_sparse <- as(as.matrix(monocle_matrix), 'sparseMatrix')

# 创建注释数据框
pd <- new("AnnotatedDataFrame", data = monocle_sample)
fd <- new("AnnotatedDataFrame", data = monocle_geneAnn)

# 创建CellDataSet对象
cds <- newCellDataSet(data_sparse, phenoData = pd, featureData = fd)

# 修改聚类列名和格式
names(pData(cds))[names(pData(cds)) == "seurat_clusters"] <- "Cluster"
pData(cds)[, "Cluster"] <- paste0("C", pData(cds)[, "Cluster"])

message(sprintf("  CDS对象创建成功，包含 %d 个细胞和 %d 个基因", ncol(cds), nrow(cds)))

# ----------- 3. 轨迹分析 -----------
message("  开始轨迹分析计算...")

# 估计尺度因子和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 设置用于排序的基因
trajectory_genes_valid <- intersect(trajectory_genes, rownames(cds))
message(sprintf("  有效的轨迹基因数量: %d", length(trajectory_genes_valid)))

# 优化基因数量
if(length(trajectory_genes_valid) > 2000) {
  message("  基因数量过多，筛选top 2000个基因...")
  trajectory_genes_valid <- trajectory_genes_valid[1:2000]
}

cds <- setOrderingFilter(cds, trajectory_genes_valid)

# 降维和细胞排序
message("  进行DDRTree降维...")
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')

message("  进行细胞排序...")
cds <- orderCells(cds)

message("  轨迹分析计算完成！")

# ----------- 4. 生成轨迹图 -----------
message("  生成轨迹分析图表...")

# 图1：按State着色的轨迹图
message("    生成State轨迹图...")
pdf("trajectory.State.pdf", width = 6.5, height = 6)
p1 <- plot_cell_trajectory(cds, color_by = "State") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Cell Trajectory by State")
print(p1)
dev.off()

# 图2：按Pseudotime着色的轨迹图
message("    生成Pseudotime轨迹图...")
pdf("trajectory.Pseudotime.pdf", width = 6.5, height = 6)
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Cell Trajectory by Pseudotime")
print(p2)
dev.off()

# 图3：按Cluster着色的轨迹图
message("    生成Cluster轨迹图...")
pdf("trajectory.cluster.pdf", width = 6.5, height = 6)
p3 <- plot_cell_trajectory(cds, color_by = "Cluster") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Cell Trajectory by Cluster") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
print(p3)
dev.off()

# 图4：按分组着色的轨迹图（对照vs实验）
message("    生成分组轨迹图...")
pdf("trajectory.Group.pdf", width = 6.5, height = 6)
p4 <- plot_cell_trajectory(cds, color_by = "orig.ident") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Cell Trajectory by Group (Control vs Treatment)")
print(p4)
dev.off()

# ==============================================================================
# 拟时序热图分析 - 绘制指定基因的拟时序表达热图
# ==============================================================================

message("\n开始生成拟时序热图...")

# 使用脚本开头定义的genes_to_analyze: ALDH2, EPB41L3, CTSO
message(sprintf("  指定的基因: %s", paste(genes_to_analyze, collapse = ", ")))

# 生成拟时序热图（使用参考文件ck.txt中的方法）
message("  正在生成拟时序热图...")

pdf("pseudotime_heatmap.pdf", width = 10, height = 8)
plot_pseudotime_heatmap(cds[genes_to_analyze,],
                        num_clusters = length(unique(pData(cds)$Cluster)),
                        show_rownames = TRUE,
                        cores = 4,
                        return_heatmap = TRUE,
                        hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

message("  拟时序热图已保存: pseudotime_heatmap.pdf")

# 生成基因表达动态曲线图（散点+拟合曲线）
message("  正在生成基因表达动态曲线图...")

pdf("Pseudotime_Gene_Expression_Curves.pdf", width = 6 * length(genes_to_analyze), height = 6)

# 使用monocle的plot_genes_in_pseudotime函数
p_curves <- plot_genes_in_pseudotime(cds[genes_to_analyze,],
                                    color_by = "Cluster",
                                    min_expr = 0.1,
                                    ncol = length(genes_to_analyze)) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Gene Expression Dynamics Along Cell Development Trajectory",
    x = "Pseudotime (Cell Development Progress)",
    y = "Expression Level"
  ) +
  scale_color_manual(values = c("C1" = "#FF6B6B", "C2" = "#4ECDC4", "C3" = "#45B7D1",
                                "C4" = "#96CEB4", "C5" = "#FFEAA7", "C6" = "#DDA0DD"))

print(p_curves)
dev.off()

message("  基因表达动态曲线图已保存: Pseudotime_Gene_Expression_Curves.pdf")

# 手动绘制自定义样式的基因表达曲线图
message("  正在生成自定义样式的基因表达曲线图...")

# 准备数据
pseudotime_data <- pData(cds)
pseudotime_data$Pseudotime <- pseudotime_data$Pseudotime

# 提取基因表达数据
expr_data <- as.matrix(exprs(cds[genes_to_analyze,]))

# 创建绘图数据框
plot_data <- data.frame(
  Cell = rep(colnames(expr_data), length(genes_to_analyze)),
  Pseudotime = rep(pseudotime_data$Pseudotime, length(genes_to_analyze)),
  Gene = rep(genes_to_analyze, each = ncol(expr_data)),
  Expression = as.vector(expr_data),
  Cluster = rep(pseudotime_data$Cluster, length(genes_to_analyze))
)

# 移除NA值
plot_data <- plot_data[!is.na(plot_data$Pseudotime) & !is.na(plot_data$Expression),]

pdf("Custom_Gene_Expression_Curves.pdf", width = 6 * length(genes_to_analyze), height = 6)

# 使用ggplot2绘制自定义图形
p_custom <- ggplot(plot_data, aes(x = Pseudotime, y = Expression, color = Gene)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, alpha = 0.2, size = 1.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = length(genes_to_analyze)) +
  scale_color_manual(values = c("ALDH2" = "#FF6B6B", "EPB41L3" = "#4ECDC4")) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Gene Expression Dynamics Along Cell Development Trajectory",
    x = "Pseudotime (Cell Development Progress)",
    y = "Expression Level",
    color = "Gene"
  )

print(p_custom)
dev.off()

message("  自定义基因表达曲线图已保存: Custom_Gene_Expression_Curves.pdf")

message("\n拟时序热图分析完成！")
message("生成的文件包括：")
message("  1. pseudotime_heatmap.pdf - 指定基因的拟时序热图（增大字体）")
message("  2. Pseudotime_Gene_Expression_Curves.pdf - 基因表达动态曲线图（横向排列）")
message("  3. Custom_Gene_Expression_Curves.pdf - 自定义样式曲线图（横向排列）")