#先运行前642行代码
# 核心单细胞分析包
library(Seurat)        # 单细胞数据分析核心包
library(limma)         # 差异表达分析
library(SingleR)       # 单细胞类型注释
library(celldex)       # 细胞类型注释数据集

# 数据处理包
library(dplyr)         # 数据处理和管道操作
library(magrittr)      # 管道操作符
library(tidyr)         # 数据整理
library(stringr)       # 字符串处理

# 可视化包
library(ggplot2)       # 基础绘图
library(ggpubr)        # 统计绘图和发表级图形
library(RColorBrewer)  # 调色板
library(viridis)       # 现代配色方案
library(scales)        # 图形标度
library(patchwork)     # 图形组合

# 统计分析包
library(rstatix)       # 统计检验

# 交互和辅助包
library(DT)            # 交互式表格
library(progress)      # 进度条显示

# 文件输出包
library(openxlsx)      # Excel文件读写

# 基因聚类富集分析包
library(ClusterGVis)   # 基因聚类可视化
library(org.Hs.eg.db)  # 人类基因注释数据库
library(clusterProfiler) # 富集分析
library(ComplexHeatmap) # 复杂热图
library(ggsci)         # 科学期刊配色方案

message("所有必需的R包已成功加载完成！")

# ==============================================================================
# 02. 参数设置和工作目录
# ==============================================================================

# 设置工作目录 - 请根据实际路径修改
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\09.单细胞分析-聚类注释 - 第二个疾病"

# 创建结果输出目录结构
output_dirs <- list(
  basic_analysis = "01.Basic_Analysis_Results",
  umap_plots = "02.UMAP_Clustering_Plots",
  marker_genes = "03.Marker_Gene_Analysis",
  heatmaps = "04.Heatmap_Visualization",
  cell_annotation = "05.Cell_Type_Annotation",
  differential_genes = "06.Differential_Gene_Analysis",
  statistics_reports = "07.Statistical_Reports",
  proportion_analysis = "08.Cell_Proportion_Analysis"
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
  pcSelect = 14,                      # PCA主成分数选择
  logFCfilter = 1,                  # logFC筛选阈值
  adjPvalFilter = 0.05,               # 调整后p值筛选阈值
  min_cells = 3,                      # 基因至少在多少个细胞中表达
  min_features = 200,                 # 细胞至少表达多少个基因
  cluster_resolution = 0.4,           # 聚类分辨率
  n_variable_features = 2000          # 高变基因数量
)

# ==============================================================================
# 03. 数据读取和Seurat对象创建
# ==============================================================================

message("正在读取单细胞数据...")

# 读取数据文件
single_cell_data <- readRDS("filtered_single_cell_data.rds")

# 数据质量检查
if (ncol(single_cell_data) == 0) {
  stop("数据文件为空，请检查文件内容！")
}

message(sprintf("数据包含 %d 个基因和 %d 个细胞",
                nrow(single_cell_data), ncol(single_cell_data)))

# 创建Seurat对象
Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
  counts = single_cell_data,
  min.cells = analysis_params$min_cells,
  min.features = analysis_params$min_features
)

message("Seurat对象创建完成")

# ==============================================================================
# 04. 数据预处理
# ==============================================================================

message("开始数据预处理...")

# 优化的预处理流程
preprocess_seurat <- function(seurat_obj, params) {
  # 检查是否有progress包，如果有则创建进度条
  if ("progress" %in% loadedNamespaces()) {
    pb <- progress::progress_bar$new(
      total = 7,
      format = "  数据处理进度 [:bar] :percent | 耗时: :elapsed | 预计剩余: :eta"
    )
    use_progress <- TRUE
  } else {
    use_progress <- FALSE
    message("开始数据预处理...")
  }

  # 数据归一化
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  if (use_progress) pb$tick() else message("1/7 数据归一化完成")

  # 寻找高变基因
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj,
    selection.method = "vst",
    nfeatures = params$n_variable_features
  )
  if (use_progress) pb$tick() else message("2/7 高变基因识别完成")

  # 数据标准化
  seurat_obj <- ScaleData(seurat_obj)
  if (use_progress) pb$tick() else message("3/7 数据标准化完成")

  # PCA分析
  seurat_obj <- RunPCA(
    object = seurat_obj,
    npcs = 20,
    features = VariableFeatures(object = seurat_obj)
  )
  if (use_progress) pb$tick() else message("4/7 PCA分析完成")

  # 邻居查找
  seurat_obj <- FindNeighbors(
    object = seurat_obj,
    dims = 1:params$pcSelect
  )
  if (use_progress) pb$tick() else message("5/7 邻居查找完成")

  # 聚类分析
  seurat_obj <- FindClusters(
    object = seurat_obj,
    resolution = params$cluster_resolution
  )
  if (use_progress) pb$tick() else message("6/7 聚类分析完成")

  # UMAP降维
  seurat_obj <- RunUMAP(
    object = seurat_obj,
    dims = 1:params$pcSelect
  )
  if (use_progress) pb$tick() else message("7/7 UMAP降维完成")

  message("数据预处理完成")
  return(seurat_obj)
}

Peripheral_Blood_Mononuclear_Cells <- preprocess_seurat(
  Peripheral_Blood_Mononuclear_Cells,
  analysis_params
)

# ==============================================================================
# 05. 基础可视化
# ==============================================================================

message("生成基础可视化图形...")

# 优化的可视化函数
create_umap_plot <- function(seurat_obj, filename, title = "UMAP Clustering Results") {
  pdf(file = filename, width = 8, height = 6)
  print(
    DimPlot(
      object = seurat_obj,
      reduction = "umap",
      pt.size = 1.5,
      label = TRUE,
      label.size = 4
    ) +
      ggtitle(title) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.3)
      )
  )
  dev.off()
  message(sprintf("保存UMAP图: %s", filename))
}

create_umap_plot(Peripheral_Blood_Mononuclear_Cells, file.path(output_dirs$umap_plots, "UMAP_Clusters.pdf"))

# 生成分组UMAP图（对照组vs疾病组）
message("生成分组UMAP图...")

# 确保有Type分组信息
if (!("Type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))) {
  Type <- gsub("(.*?)\\..*", "\\1", colnames(Peripheral_Blood_Mononuclear_Cells))
  names(Type) <- colnames(Peripheral_Blood_Mononuclear_Cells)
  Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
    object = Peripheral_Blood_Mononuclear_Cells,
    metadata = Type,
    col.name = "Type"
  )
  message("已添加Type分组信息")
}

# 创建分面UMAP图（对照组和疾病组并排展示）
pdf(file = file.path(output_dirs$umap_plots, "UMAP_Clusters_Split_by_Group.pdf"), width = 14, height = 6)
print(
  DimPlot(
    Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1.2,
    label = TRUE,
    label.size = 3.5,
    split.by = "Type",
    ncol = 2
  ) +
    ggtitle("UMAP Clustering: Control vs Disease") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.3),
      strip.text = element_text(face = "bold", size = 12)
    )
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

# ==============================================================================
# 06. Marker基因分析
# ==============================================================================

message("进行Marker基因分析...")

# 寻找所有聚类的marker基因
find_all_markers <- function(seurat_obj, params) {
  markers <- FindAllMarkers(
    object = seurat_obj,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = params$logFCfilter
  )

  # 筛选显著的marker基因
  sig_markers <- markers[
    (abs(markers$avg_log2FC) > params$logFCfilter &
       markers$p_val_adj < params$adjPvalFilter),
  ]

  # 保存结果
  write.table(
    sig_markers,
    file = file.path(output_dirs$marker_genes, "clusterMarkers.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message(sprintf("找到 %d 个显著marker基因", nrow(sig_markers)))
  return(list(all_markers = markers, sig_markers = sig_markers))
}

marker_results <- find_all_markers(Peripheral_Blood_Mononuclear_Cells, analysis_params)

# ==============================================================================
# 07. 热图可视化
# ==============================================================================

message("生成热图...")

# 优化的热图生成
create_heatmap <- function(seurat_obj, markers_df, filename = file.path(output_dirs$heatmaps, "umapHeatmap.pdf")) {
  top10 <- markers_df %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

  pdf(file = filename, width = 12, height = 9)
  print(
    DoHeatmap(
      object = seurat_obj,
      features = top10$gene
    ) +
      NoLegend()
  )
  dev.off()
  message(sprintf("保存热图: %s", filename))
}

create_heatmap(Peripheral_Blood_Mononuclear_Cells, marker_results$all_markers)

# ==============================================================================
# 08. 细胞类型标记基因分析和气泡图生成（Manual_Annotation_Analysis）
# ==============================================================================

message("进行细胞类型标记基因分析...")

# 创建标记基因分析结果文件夹
annotation_dir <- file.path(output_dirs$cell_annotation, "Manual_Annotation_Analysis")
if (!dir.exists(annotation_dir)) {
  dir.create(annotation_dir, recursive = TRUE)
}

# 读取细胞类型标记基因注释文件
marker_file <- "Cell Annotation.txt"
if (!file.exists(marker_file)) {
  message(paste("警告：标记基因文件", marker_file, "不存在，跳过标记基因分析"))
} else {
  message(paste("读取标记基因文件：", marker_file))

  # 读取标记基因文件
  refMarker <- read.table(marker_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

  # 解析基因列表
  genes <- list()
  for(i in 1:nrow(refMarker)){
    celltype_name <- refMarker[i, 1]
    # 跳过标题行
    if(celltype_name == "Cell Type (English)" || celltype_name == "Cell Type") next

    gene_string <- refMarker[i, 2]
    if(!is.na(gene_string) && gene_string != "Marker Genes") {
      # 移除引号并分割基因
      gene_string <- gsub('"', '', gene_string)
      genes[[celltype_name]] <- unlist(strsplit(gene_string, "\\,|\\;|\\|"))
      genes[[celltype_name]] <- trimws(genes[[celltype_name]])
      genes[[celltype_name]] <- genes[[celltype_name]][genes[[celltype_name]] != ""]
    }
  }

  # 移除空的基因列表
  genes <- genes[sapply(genes, length) > 0]

  message("解析的细胞类型和基因：")
  for(ct in names(genes)) {
    message(sprintf("  %s: %d个基因 (%s)", ct, length(genes[[ct]]),
                    paste(head(genes[[ct]], 3), collapse = ", ")))
  }

  if(length(genes) == 0) {
    message("警告：没有找到有效的基因列表，请检查文件格式")
  } else {

    # 验证基因是否在数据中存在
    all_marker_genes <- unique(unlist(genes))
    valid_marker_genes <- all_marker_genes[all_marker_genes %in% rownames(Peripheral_Blood_Mononuclear_Cells)]
    missing_genes <- setdiff(all_marker_genes, valid_marker_genes)

    message(sprintf("总标记基因: %d, 在数据中找到: %d, 缺失: %d",
                    length(all_marker_genes), length(valid_marker_genes), length(missing_genes)))

    if(length(missing_genes) > 0) {
      write.csv(data.frame(Missing_Genes = missing_genes),
                file = file.path(annotation_dir, "Missing_Marker_Genes.csv"),
                row.names = FALSE)
    }

    # 过滤基因列表，只保留存在的基因
    genes_filtered <- lapply(genes, function(x) x[x %in% rownames(Peripheral_Blood_Mononuclear_Cells)])
    genes_filtered <- genes_filtered[sapply(genes_filtered, length) > 0]

    if(length(genes_filtered) == 0) {
      message("警告：没有找到有效的标记基因，跳过分析")
    } else {

      # 绘制标记基因的气泡图
      message("绘制标记基因气泡图...")

      total_genes_for_plot <- sum(sapply(genes_filtered, length))
      if(total_genes_for_plot > 0) {

        # 限制每个细胞类型最多显示前8个基因
        genes_for_plot <- lapply(genes_filtered, function(x) head(x, 8))

        # 绘制气泡图 - 添加错误处理以兼容新版本ggplot2
        tryCatch({
          # 先尝试创建DotPlot对象，如果成功再打开PDF
          p_dot <- DotPlot(Peripheral_Blood_Mononuclear_Cells,
                           features = genes_for_plot,
                           group.by = "seurat_clusters",
                           dot.scale = 8,
                           cols = c("lightgrey", "red")) +
            theme_classic(base_size = 12) +
            theme(
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 9),
              axis.text.y = element_text(face = "bold", size = 10),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              legend.title = element_text(face = "bold"),
              strip.text = element_text(face = "bold", size = 10),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(color = "black", size = 0.5),
              axis.ticks = element_line(color = "black", size = 0.3)
            ) +
            labs(
              title = "Marker Genes Expression Across Clusters (For Manual Annotation)",
              x = "Genes",
              y = "Clusters"
            ) +
            coord_flip()

          # 如果DotPlot创建成功，再打开PDF并绘制
          pdf(file.path(annotation_dir, "Marker_Gene_DotPlot_for_Manual_Annotation.pdf"), width = 30, height = 12)
          print(p_dot)
          dev.off()
          message("已保存: Marker_Gene_DotPlot_for_Manual_Annotation.pdf")

        }, error = function(e) {
          # 关闭可能打开的PDF设备
          if (dev.cur() != 1) dev.off()

          # 删除可能生成的空PDF文件
          empty_pdf <- file.path(annotation_dir, "Marker_Gene_DotPlot_for_Manual_Annotation.pdf")
          if (file.exists(empty_pdf)) {
            file.remove(empty_pdf)
            message("已删除空的PDF文件")
          }

          if(any(grepl("facets.*defunct", e$message, ignore.case = TRUE))) {
            message("警告: DotPlot函数遇到ggplot2版本兼容性问题 (facets参数已废弃)")
            message("建议更新Seurat包到最新版本或降级ggplot2版本")
            message("跳过该气泡图的生成，继续执行后续分析...")
          } else {
            message(sprintf("生成气泡图时出现其他错误: %s", e$message))
          }
        })

        # 生成汇总气泡图（不显示具体细胞类型名字）
        message("生成汇总标记基因气泡图...")

        # 收集所有标记基因（不分细胞类型）
        all_marker_genes <- unique(unlist(genes_for_plot))
        # 限制基因数量，避免图形过于拥挤
        if(length(all_marker_genes) > 30) {
          all_marker_genes <- head(all_marker_genes, 30)
          message(sprintf("基因数量过多，限制显示前30个基因"))
        }

        # 生成汇总气泡图 - 添加错误处理
        tryCatch({
          # 先创建DotPlot对象
          p_comprehensive <- DotPlot(Peripheral_Blood_Mononuclear_Cells,
                                   features = all_marker_genes,
                                   group.by = "seurat_clusters",
                                   dot.scale = 6,
                                   cols = c("#E8E8E8", "#4A90E2", "#2E5BBA", "#1A365D")) +
          theme_classic(base_size = 11) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 9),
            axis.text.y = element_text(face = "bold", size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 20)),
            plot.subtitle = element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 15), color = "gray50"),
            legend.title = element_text(face = "bold", size = 11),
            legend.text = element_text(size = 9),
            legend.position = "bottom",
            legend.box = "horizontal",
            panel.grid.major = element_line(color = "gray92", size = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 0.8),
            axis.line = element_line(color = "black", size = 0.5),
            axis.ticks = element_line(color = "black", size = 0.3),
            plot.margin = ggplot2::margin(20, 20, 20, 20)
          ) +
          labs(
            title = "Comprehensive Marker Genes Expression Analysis",
            subtitle = paste0("Expression patterns of ", length(all_marker_genes), " marker genes across all clusters"),
            x = "Marker Genes",
            y = "Clusters",
            caption = "Dot size represents percentage of expressing cells; Color intensity represents average expression level"
          ) +
          guides(
            color = guide_colorbar(
              title = "Average\nExpression",
              title.position = "top",
              title.hjust = 0.5,
              barwidth = 12,
              barheight = 1
            ),
            size = guide_legend(
              title = "Percent\nExpressed",
              title.position = "top",
              title.hjust = 0.5,
              override.aes = list(color = "black")
            )
          )

        # 如果DotPlot创建成功，再打开PDF并绘制
        pdf(file.path(annotation_dir, "Comprehensive_Marker_Genes_DotPlot.pdf"), width = 16, height = 10)
        print(p_comprehensive)
        dev.off()
        message("已保存: Comprehensive_Marker_Genes_DotPlot.pdf")

        }, error = function(e) {
          if (dev.cur() != 1) dev.off()

          # 删除可能生成的空PDF文件
          empty_pdf <- file.path(annotation_dir, "Comprehensive_Marker_Genes_DotPlot.pdf")
          if (file.exists(empty_pdf)) {
            file.remove(empty_pdf)
            message("已删除空的Comprehensive DotPlot PDF文件")
          }

          if(any(grepl("facets.*defunct", e$message, ignore.case = TRUE))) {
            message("警告: Comprehensive DotPlot遇到ggplot2版本兼容性问题，跳过生成")
          } else {
            message(sprintf("生成Comprehensive DotPlot时出现错误: %s", e$message))
          }
        })

        # 保存使用的基因列表
        write.table(all_marker_genes,
                   file.path(annotation_dir, "Comprehensive_DotPlot_Genes_List.txt"),
                   quote = FALSE, row.names = FALSE, col.names = FALSE)
        message("已保存基因列表: Comprehensive_DotPlot_Genes_List.txt")
      }

      # 计算每个细胞类型标记基因的聚类富集分数
      message("计算细胞类型-聚类富集分数...")

      enrichment_results <- data.frame()
      cluster_levels <- sort(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))

      for(celltype in names(genes_filtered)) {
        marker_genes <- genes_filtered[[celltype]]
        message(sprintf("  分析细胞类型: %s (%d个基因)", celltype, length(marker_genes)))

        for(cluster in cluster_levels) {
          cluster_cells <- colnames(Peripheral_Blood_Mononuclear_Cells)[Peripheral_Blood_Mononuclear_Cells$seurat_clusters == cluster]

          if(length(cluster_cells) == 0) next

          tryCatch({
            if(length(marker_genes) == 1) {
              gene_expr <- FetchData(Peripheral_Blood_Mononuclear_Cells, vars = marker_genes, cells = cluster_cells)[,1]
            } else {
              gene_expr_matrix <- FetchData(Peripheral_Blood_Mononuclear_Cells, vars = marker_genes, cells = cluster_cells)
              gene_expr <- rowMeans(gene_expr_matrix, na.rm = TRUE)
            }

            cluster_score <- data.frame(
              CellType = celltype,
              Cluster = as.character(cluster),
              Mean_Expression = mean(gene_expr, na.rm = TRUE),
              Median_Expression = median(gene_expr, na.rm = TRUE),
              Positive_Cells = sum(gene_expr > 0, na.rm = TRUE),
              Total_Cells = length(gene_expr),
              Positive_Percentage = round(sum(gene_expr > 0, na.rm = TRUE) / length(gene_expr) * 100, 2),
              Gene_Count = length(marker_genes),
              stringsAsFactors = FALSE
            )

            enrichment_results <- rbind(enrichment_results, cluster_score)
          }, error = function(e) {
            message(sprintf("    聚类 %s 分析出错: %s", cluster, e$message))
          })
        }
      }

      if(nrow(enrichment_results) > 0) {
        # 保存富集结果
        write.csv(enrichment_results,
                  file = file.path(annotation_dir, "CellType_Cluster_Enrichment_Scores.csv"),
                  row.names = FALSE)

        # 生成推荐的聚类-细胞类型对应表
        message("生成聚类-细胞类型推荐对应表...")

        cluster_celltype_prediction <- enrichment_results %>%
          group_by(Cluster) %>%
          slice_max(Mean_Expression, n = 1, with_ties = FALSE) %>%
          select(Cluster, Predicted_CellType = CellType, Max_Score = Mean_Expression,
                 Positive_Percentage, Gene_Count) %>%
          ungroup() %>%
          arrange(as.numeric(Cluster))

        write.csv(cluster_celltype_prediction,
                  file = file.path(annotation_dir, "Cluster_CellType_Prediction.csv"),
                  row.names = FALSE)

        message("\n=== 基于标记基因的聚类注释推荐 ===")
        for(i in 1:nrow(cluster_celltype_prediction)) {
          message(sprintf("聚类 %s: %s (得分: %.3f, 阳性比例: %.1f%%)",
                          cluster_celltype_prediction$Cluster[i],
                          cluster_celltype_prediction$Predicted_CellType[i],
                          cluster_celltype_prediction$Max_Score[i],
                          cluster_celltype_prediction$Positive_Percentage[i]))
        }
      }
    }
  }
}

# ==============================================================================
# 09. 细胞类型手动注释（基于标记基因气泡图）
# ==============================================================================

根据标记基因分析结果进行的手动注释
# 注：每个聚类的注释来自clusterProfiler富集分析，括号内为得分和阳性标记基因比例

cluster_names <- c(
  "0" = "T cells",
  "1" = "Monocytes",
  "2" = "Monocytes",
  "3" = "Monocytes",
  "4" = "Monocytes",
  "5" = "Dendritic Cells",
  "6" = "Monocytes",
  "7" = "Monocytes",
  "8" = "Monocytes",
  "9" = "Monocytes",
  "10" = "Monocytes",
  "11" = "Monocytes",
  "12" = "Monocytes",
  "13" = "Monocytes",
  "14" = "Monocytes",
  "15" = "Monocytes",
  "16" = "Monocytes",
  "17" = "Monocytes",
  "18" = "Dendritic Cells",
  "19" = "Monocytes",
  "20" = "Monocytes"
)

# 获取当前的聚类信息
current_clusters <- sort(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))
message(sprintf("当前聚类: %s", paste(current_clusters, collapse = ", ")))

# 只为存在的聚类分配注释
valid_annotations <- cluster_names[names(cluster_names) %in% as.character(current_clusters)]

# 应用手动注释
current_idents <- Idents(Peripheral_Blood_Mononuclear_Cells)
new_idents <- as.character(current_idents)

for(cluster_id in names(valid_annotations)) {
  new_idents[current_idents == cluster_id] <- valid_annotations[cluster_id]
}

# 更新细胞身份
Idents(Peripheral_Blood_Mononuclear_Cells) <- factor(new_idents)

# 添加细胞类型到元数据
Peripheral_Blood_Mononuclear_Cells$cell_type <- Idents(Peripheral_Blood_Mononuclear_Cells)

# 保存聚类-细胞类型对应表
annotation_dir <- file.path(output_dirs$cell_annotation, "Manual_Annotation_Analysis")
if (!dir.exists(annotation_dir)) {
  dir.create(annotation_dir, recursive = TRUE)
}

cluster_annotation_table <- data.frame(
  Cluster = as.character(current_clusters),
  Manual_CellType = cluster_names[as.character(current_clusters)],
  stringsAsFactors = FALSE
)

write.csv(cluster_annotation_table,
          file = file.path(annotation_dir, "Manual_Cluster_CellType_Annotation.csv"),
          row.names = FALSE)

message("手动注释完成！")
message(sprintf("注释文件已保存到: %s", annotation_dir))

# 生成细胞ID到细胞类型映射表
message("生成细胞ID到细胞类型映射表...")

# 创建详细的细胞映射表
cell_id_to_celltype_mapping <- data.frame(
  Cell_ID = colnames(Peripheral_Blood_Mononuclear_Cells),
  Cluster_ID = as.character(Peripheral_Blood_Mononuclear_Cells$seurat_clusters),
  Cell_Type = as.character(Idents(Peripheral_Blood_Mononuclear_Cells)),
  stringsAsFactors = FALSE
)

# 添加额外信息
cell_id_to_celltype_mapping$nFeature_RNA <- Peripheral_Blood_Mononuclear_Cells$nFeature_RNA
cell_id_to_celltype_mapping$nCount_RNA <- Peripheral_Blood_Mononuclear_Cells$nCount_RNA

# 检查线粒体百分比列名并添加
if ("percent.mt" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  cell_id_to_celltype_mapping$percent.mito <- Peripheral_Blood_Mononuclear_Cells$percent.mt
} else if ("percent.mito" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)) {
  cell_id_to_celltype_mapping$percent.mito <- Peripheral_Blood_Mononuclear_Cells$percent.mito
} else {
  message("警告：未找到线粒体百分比列，跳过此项")
  cell_id_to_celltype_mapping$percent.mito <- NA
}

# 按聚类ID排序
cell_id_to_celltype_mapping <- cell_id_to_celltype_mapping[order(as.numeric(cell_id_to_celltype_mapping$Cluster_ID)), ]

# 保存到指定位置
write.csv(cell_id_to_celltype_mapping,
          file = file.path(workDir, "Cell_ID_to_CellType_Mapping.csv"),
          row.names = FALSE)

# 生成细胞类型统计摘要
celltype_summary <- cell_id_to_celltype_mapping %>%
  group_by(Cell_Type, Cluster_ID) %>%
  summarise(
    Cell_Count = n(),
    Mean_Features = round(mean(nFeature_RNA, na.rm = TRUE), 2),
    Mean_Counts = round(mean(nCount_RNA, na.rm = TRUE), 2),
    Mean_Mito_Percent = round(mean(percent.mito, na.rm = TRUE), 2),
    .groups = 'drop'
  ) %>%
  arrange(Cluster_ID)

# 保存细胞类型统计摘要
write.csv(celltype_summary,
          file = file.path(workDir, "CellType_Summary_Statistics.csv"),
          row.names = FALSE)

# 输出映射表信息
total_cells <- nrow(cell_id_to_celltype_mapping)
unique_celltypes <- unique(cell_id_to_celltype_mapping$Cell_Type)
unique_clusters <- unique(cell_id_to_celltype_mapping$Cluster_ID)

message(sprintf("细胞ID到细胞类型映射表已生成"))
message(sprintf("- 总细胞数: %d", total_cells))
message(sprintf("- 细胞类型数: %d (%s)", length(unique_celltypes), paste(unique_celltypes, collapse = ", ")))
message(sprintf("- 聚类数: %d", length(unique_clusters)))
message(sprintf("- 映射表保存位置: %s", file.path(workDir, "Cell_ID_to_CellType_Mapping.csv")))
message(sprintf("- 统计摘要保存位置: %s", file.path(workDir, "CellType_Summary_Statistics.csv")))

# 显示每个细胞类型的细胞数量
message("\n=== 各细胞类型细胞数量统计 ===")
celltype_counts <- table(cell_id_to_celltype_mapping$Cell_Type)
for(celltype in names(celltype_counts)) {
  message(sprintf("%s: %d 个细胞", celltype, celltype_counts[celltype]))
}

# ==============================================================================
# 10. 注释后可视化
# ==============================================================================

message("生成注释后的可视化图...")

# 细胞类型注释UMAP图
pdf(file = file.path(output_dirs$cell_annotation, "cellTypeAnnotation.pdf"), width = 10, height = 8)
print(
  DimPlot(
    Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1.5,
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Manual Cell Type Annotations") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.3)
    )
)
dev.off()

# 添加分组信息
Type <- gsub("(.*?)\\..*", "\\1", colnames(Peripheral_Blood_Mononuclear_Cells))
names(Type) <- colnames(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
  object = Peripheral_Blood_Mononuclear_Cells,
  metadata = Type,
  col.name = "Type"
)

# 分组可视化
pdf(file = file.path(output_dirs$cell_annotation, "groupComparison.pdf"), width = 14, height = 6)
print(
  DimPlot(
    Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1,
    label = TRUE,
    split.by = "Type"
  ) +
    ggtitle("Cell Types: Control vs Disease") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.3)
    )
)
dev.off()

# ==============================================================================
# 11. 差异表达分析
# ==============================================================================

message("进行组间差异表达分析...")

# 优化的差异表达分析函数
perform_differential_analysis <- function(seurat_obj, cell_ann, type_info, params) {
  # 添加组合标签
  groups <- paste0(type_info, "_", cell_ann)
  names(groups) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = groups,
    col.name = "group"
  )

  # 动态获取分组名称
  unique_types <- unique(type_info)
  if (length(unique_types) != 2) {
    message("警告：数据中的分组数量不是2个，跳过差异表达分析")
    return(seurat_obj)
  }

  control_group <- unique_types[1]  # 第一个作为对照组
  treatment_group <- unique_types[2]  # 第二个作为处理组

  message(sprintf("差异表达分析: %s (对照) vs %s (处理)", control_group, treatment_group))

  # 为每种细胞类型进行差异分析
  unique_cell_types <- unique(cell_ann)

  for (cellName in unique_cell_types) {
    con_name <- paste0(control_group, "_", cellName)
    treat_name <- paste0(treatment_group, "_", cellName)

    # 检查标签是否存在
    if (!(con_name %in% seurat_obj$group) | !(treat_name %in% seurat_obj$group)) {
      message(sprintf("跳过 %s - 缺少对照组或处理组", cellName))
      next
    }

    # 获取细胞
    con_cells <- WhichCells(seurat_obj, expression = group == con_name)
    treat_cells <- WhichCells(seurat_obj, expression = group == treat_name)

    if (length(con_cells) >= 3 & length(treat_cells) >= 3) {
      message(sprintf("分析 %s - 对照组: %d 细胞, 处理组: %d 细胞",
                      cellName, length(con_cells), length(treat_cells)))

      tryCatch({
        pbmc_markers <- FindMarkers(
          seurat_obj,
          ident.1 = treat_cells,
          ident.2 = con_cells,
          group.by = 'group',
          logfc.threshold = 0.1
        )

        sig_markers_group <- pbmc_markers[
          (abs(pbmc_markers$avg_log2FC) > params$logFCfilter &
             pbmc_markers$p_val_adj < params$adjPvalFilter),
        ]

        if (nrow(sig_markers_group) > 0) {
          sig_markers_group <- cbind(Gene = rownames(sig_markers_group), sig_markers_group)

          # 清理文件名中的特殊字符
          safe_cell_name <- gsub("[^A-Za-z0-9_]", "_", cellName)

          write.table(
            sig_markers_group,
            file = file.path(output_dirs$differential_genes, paste0(safe_cell_name, "_diffGenes.txt")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
          )

          message(sprintf("保存差异基因结果: %s_diffGenes.txt", safe_cell_name))
        } else {
          message(sprintf("未找到 %s 的显著差异基因", cellName))
        }
      }, error = function(e) {
        message(sprintf("分析 %s 时出错: %s", cellName, e$message))
      })
    } else {
      message(sprintf("跳过 %s - 细胞数量不足 (对照组:%d, 处理组:%d)",
                      cellName, length(con_cells), length(treat_cells)))
    }
  }

  return(seurat_obj)
}

# 获取手动注释的细胞类型信息
cellAnn <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
names(cellAnn) <- colnames(Peripheral_Blood_Mononuclear_Cells)

Peripheral_Blood_Mononuclear_Cells <- perform_differential_analysis(
  Peripheral_Blood_Mononuclear_Cells,
  cellAnn,
  Type,
  analysis_params
)

# ==============================================================================
# 12. 基因聚类富集分析和可视化 (ClusterGVis)
# ==============================================================================

message("开始基因聚类富集分析和可视化...")

# 创建基因聚类富集分析目录
gene_cluster_dir <- file.path(output_dirs$statistics_reports, "Gene_Cluster_Enrichment")
if (!dir.exists(gene_cluster_dir)) {
  dir.create(gene_cluster_dir, recursive = TRUE)
}

# 从已有的marker基因分析中获取top基因
if (exists("marker_results") && !is.null(marker_results$all_markers)) {
  message("使用已有的marker基因进行聚类富集分析...")

  # 获取marker基因的top基因
  top_markers <- marker_results$all_markers %>%
    group_by(cluster) %>%
    arrange(desc(abs(avg_log2FC)), p_val_adj) %>%
    slice_head(n = 20)  # 每个cluster取top 20个基因

  # 获取标准化的表达数据
  normalized_data <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, slot = "data")

  # 筛选marker基因的表达矩阵
  marker_genes <- unique(top_markers$gene)
  marker_genes <- intersect(marker_genes, rownames(normalized_data))
  message(sprintf("总共包含 %d 个标记基因用于聚类分析", length(marker_genes)))

  if (length(marker_genes) > 10) {
    # 提取标记基因的表达矩阵
    gene_expression_matrix <- as.matrix(normalized_data[marker_genes, ])

    # 计算每个细胞类型的平均表达
    cell_types <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
    unique_cell_types <- unique(cell_types)

    # 构建平均表达矩阵
    avg_expression <- matrix(0, nrow = length(marker_genes), ncol = length(unique_cell_types))
    rownames(avg_expression) <- marker_genes
    colnames(avg_expression) <- unique_cell_types

    for (ct in unique_cell_types) {
      cells_of_type <- which(cell_types == ct)
      if (length(cells_of_type) > 1) {
        avg_expression[, ct] <- rowMeans(gene_expression_matrix[, cells_of_type])
      } else if (length(cells_of_type) == 1) {
        avg_expression[, ct] <- gene_expression_matrix[, cells_of_type]
      }
    }

    message(sprintf("已构建平均表达矩阵，维度: %d genes x %d cell types",
                    nrow(avg_expression), ncol(avg_expression)))

    # 使用ClusterGVis进行基因聚类分析
    message("正在进行基因聚类分析...")

    # 确保数据格式正确
    avg_expression <- as.data.frame(avg_expression)

    # 设置聚类数量（可以根据需要调整）
    cluster_number <- min(8, max(4, round(length(marker_genes) / 50)))

    # 进行基因聚类（使用kmeans方法）
    gene_clusters <- clusterData(obj = avg_expression,
                               cluster.method = "kmeans",
                               cluster.num = cluster_number)

    message("基因聚类完成，使用kmeans方法")

    # 进行GO富集分析
    message("正在进行GO富集分析...")

    gene_enrich <- enrichCluster(object = gene_clusters,
                               OrgDb = org.Hs.eg.db,
                               type = "BP",  # 生物过程
                               organism = "hsa",
                               pvalueCutoff = 0.05,
                               topn = 5,
                               seed = 5201314)

    message("GO富集分析完成")

    # 选择要标记的重要基因
    predefined_genes <- c("CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD14", "FCGR3A",
                         "LYZ", "MS4A1", "GNLY", "PRF1", "CST3", "FCER1A")

    # 筛选在我们数据中存在的基因
    mark_genes <- intersect(predefined_genes, rownames(avg_expression))

    # 如果预定义基因数量不足，从marker基因中补充
    if(length(mark_genes) < 15) {
      additional_genes <- sample(rownames(avg_expression)[!rownames(avg_expression) %in% mark_genes],
                                min(15 - length(mark_genes), nrow(avg_expression) - length(mark_genes)))
      mark_genes <- c(mark_genes, additional_genes)
    }

    message(sprintf("将标记 %d 个基因", length(mark_genes)))

    # 生成包含富集分析的完整图形
    message("正在生成基因聚类富集注释热图...")

    # 包含富集分析的完整图形
    pdf(file.path(gene_cluster_dir, "Gene_Cluster_Enrichment_Heatmap_Complete.pdf"),
        width = 16, height = 12, onefile = FALSE)

    # 计算实际需要的颜色数量
    num_terms <- nrow(gene_enrich)
    go_colors <- rep(ggsci::pal_d3()(8), length.out = num_terms)

    visCluster(object = gene_clusters,
              plot.type = "both",
              column_names_rot = 45,
              show_row_dend = FALSE,
              markGenes = mark_genes,
              markGenes.side = "left",
              genes.gp = c('italic', fontsize = 10, col = "black"),
              annoTerm.data = gene_enrich,
              line.side = "left",
              go.col = go_colors,
              go.size = "pval",
              add.bar = TRUE,
              textbar.pos = c(0.85, 0.15))
    dev.off()

    # 仅包含热图版本
    pdf(file.path(gene_cluster_dir, "Gene_Cluster_Heatmap.pdf"),
        width = 12, height = 10, onefile = FALSE)
    visCluster(object = gene_clusters,
              plot.type = "heatmap",
              column_names_rot = 45,
              show_row_dend = FALSE,
              markGenes = mark_genes,
              markGenes.side = "left",
              genes.gp = c('italic', fontsize = 10, col = "black"))
    dev.off()

    # 保存聚类结果和富集结果
    write.table(gene_clusters$wide.res,
                file = file.path(gene_cluster_dir, "Gene_Cluster_Results.csv"),
                sep = ",", quote = FALSE, row.names = TRUE)

    write.table(gene_enrich,
                file = file.path(gene_cluster_dir, "Gene_Enrichment_Results.csv"),
                sep = ",", quote = FALSE, row.names = FALSE)

    # 保存用于聚类的表达矩阵
    write.table(avg_expression,
                file = file.path(gene_cluster_dir, "Average_Expression_Matrix.csv"),
                sep = ",", quote = FALSE, row.names = TRUE)

    # 保存标记基因列表
    write.table(data.frame(Marked_Genes = mark_genes),
                file = file.path(gene_cluster_dir, "Marked_Genes_List.csv"),
                sep = ",", quote = FALSE, row.names = FALSE)

    message("基因聚类富集分析完成！")
    message("输出文件：")
    message("- Gene_Cluster_Enrichment_Heatmap_Complete.pdf: 完整聚类富集注释热图")
    message("- Gene_Cluster_Heatmap.pdf: 基因聚类热图")
    message("- Gene_Cluster_Results.csv: 基因聚类结果")
    message("- Average_Expression_Matrix.csv: 平均表达矩阵")

  } else {
    message("标记基因数量不足，跳过基因聚类富集分析")
  }
} else {
  message("未找到marker基因分析结果，跳过基因聚类富集分析")
}

message("基因聚类富集分析模块完成")

# ==============================================================================
# 13. 对照组与实验组细胞比例分析
# ==============================================================================

message("开始对照组与实验组细胞比例分析...")

# 细胞比例分析函数
analyze_celltype_proportions <- function(seurat_obj) {

  message("计算各组细胞类型比例...")

  # 获取细胞类型和分组信息
  cell_types <- Idents(seurat_obj)
  groups <- seurat_obj@meta.data$Type

  # 创建交叉表
  cross_table <- table(cell_types, groups)

  # 计算各种比例
  # 1. 每组内各细胞类型的比例
  prop_within_group <- prop.table(cross_table, margin = 2) * 100

  # 2. 每个细胞类型在各组中的分布比例
  prop_within_celltype <- prop.table(cross_table, margin = 1) * 100

  # 3. 总体比例
  prop_total <- prop.table(cross_table) * 100

  # 转换为数据框格式
  # 组内比例数据框
  prop_within_group_df <- as.data.frame(prop_within_group)
  colnames(prop_within_group_df) <- c("CellType", "Group", "Percentage_Within_Group")

  # 细胞类型内比例数据框
  prop_within_celltype_df <- as.data.frame(prop_within_celltype)
  colnames(prop_within_celltype_df) <- c("CellType", "Group", "Percentage_Within_CellType")

  # 总体比例数据框
  prop_total_df <- as.data.frame(prop_total)
  colnames(prop_total_df) <- c("CellType", "Group", "Percentage_Total")

  # 细胞计数数据框
  counts_df <- as.data.frame(cross_table)
  colnames(counts_df) <- c("CellType", "Group", "Cell_Count")

  # 合并所有信息
  comprehensive_stats <- merge(counts_df, prop_within_group_df, by = c("CellType", "Group"))
  comprehensive_stats <- merge(comprehensive_stats, prop_within_celltype_df, by = c("CellType", "Group"))
  comprehensive_stats <- merge(comprehensive_stats, prop_total_df, by = c("CellType", "Group"))

  # 添加总计信息
  total_control <- sum(groups == "Control")
  total_treat <- sum(groups == "Treat")
  total_cells <- length(groups)

  comprehensive_stats$Total_Cells_In_Group <- ifelse(
    comprehensive_stats$Group == "Control",
    total_control,
    total_treat
  )

  comprehensive_stats$Total_Cells_Overall <- total_cells

  # 重新排列列顺序
  comprehensive_stats <- comprehensive_stats[, c(
    "CellType", "Group", "Cell_Count", "Total_Cells_In_Group",
    "Percentage_Within_Group", "Percentage_Within_CellType",
    "Percentage_Total", "Total_Cells_Overall"
  )]

  # 按细胞类型和组排序
  comprehensive_stats <- comprehensive_stats[order(comprehensive_stats$CellType, comprehensive_stats$Group), ]

  return(list(
    cross_table = cross_table,
    prop_within_group = prop_within_group,
    prop_within_celltype = prop_within_celltype,
    comprehensive_stats = comprehensive_stats
  ))
}

# 执行比例分析
proportion_results <- analyze_celltype_proportions(Peripheral_Blood_Mononuclear_Cells)

# 保存比例分析详细统计结果
write.csv(
  proportion_results$comprehensive_stats,
  file.path(output_dirs$proportion_analysis, "Comprehensive_CellType_Proportions.csv"),
  row.names = FALSE
)

# 生成细胞类型比例堆积条形图
create_celltype_proportion_barplot <- function(seurat_obj) {

  # 获取细胞类型和分组信息
  cell_types <- as.character(Idents(seurat_obj))
  groups <- seurat_obj@meta.data$Type

  message("调试信息:")
  message(sprintf("  细胞类型数: %d", length(cell_types)))
  message(sprintf("  分组数: %d", length(groups)))
  message(sprintf("  细胞总数: %d", ncol(seurat_obj)))

  # 创建数据框
  proportion_data <- data.frame(
    CellType = cell_types,
    Group = groups,
    stringsAsFactors = FALSE
  )

  message(sprintf("  比例数据框创建成功，行数: %d", nrow(proportion_data)))

  # 计算每组内各细胞类型的比例
  proportion_summary <- proportion_data %>%
    group_by(Group, CellType) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(Group) %>%
    mutate(
      Total = sum(Count),
      Proportion = Count / Total
    ) %>%
    ungroup()

  message(sprintf("  比例汇总数据框创建成功，行数: %d", nrow(proportion_summary)))

  # 自动生成细胞类型颜色
  unique_celltypes <- unique(proportion_summary$CellType)
  celltype_colors <- setNames(
    colorRampPalette(brewer.pal(min(12, length(unique_celltypes)), "Set3"))(length(unique_celltypes)),
    unique_celltypes
  )

  # 创建堆积条形图
  p <- ggplot(proportion_summary, aes(x = Group, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +

    # 设置颜色
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +

    # 设置Y轴
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      expand = c(0, 0)
    ) +

    # 主题设置
    theme_classic(base_size = 12) +
    theme(
      # 坐标轴
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),

      # 坐标轴线条
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),

      # 图例设置
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, "cm"),
      legend.key = element_rect(color = "black", linewidth = 0.3),

      # 面板设置
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # 边距
      plot.margin = ggplot2::margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
    ) +

    # 标签
    labs(y = "Ratio")

  return(p)
}

# 生成并保存堆积条形图
celltype_proportion_plot <- create_celltype_proportion_barplot(Peripheral_Blood_Mononuclear_Cells)

# 保存PDF版本
pdf(file.path(output_dirs$proportion_analysis, "CellType_Proportion_Stacked_Barplot.pdf"),
    width = 8, height = 6)
print(celltype_proportion_plot)
dev.off()

message("细胞类型比例堆积条形图已保存到: 08.Cell_Proportion_Analysis/")

# ==============================================================================
# 14. 最终结果汇总和保存
# ==============================================================================

message("生成分析统计报告...")

# 优化的统计报告函数
generate_statistics_report <- function(seurat_obj) {
  # 细胞类型统计
  cell_type_counts <- table(Idents(seurat_obj),
                            seurat_obj$Type)

  # 动态获取分组名称
  group_names <- colnames(cell_type_counts)

  # 创建基础统计表
  cell_type_stats <- data.frame(
    CellType = rownames(cell_type_counts),
    stringsAsFactors = FALSE
  )

  # 动态添加每个分组的计数
  for (group in group_names) {
    cell_type_stats[[paste0(group, "_Count")]] <- cell_type_counts[, group]
  }

  # 添加总计
  cell_type_stats$Total_Count <- rowSums(cell_type_counts)

  # 动态计算每个分组的百分比
  for (group in group_names) {
    count_col <- paste0(group, "_Count")
    pct_col <- paste0(group, "_Percentage")
    cell_type_stats[[pct_col]] <- round(
      cell_type_stats[[count_col]] / sum(cell_type_stats[[count_col]]) * 100, 2
    )
  }

  # 保存统计结果
  write.csv(cell_type_stats, file.path(output_dirs$statistics_reports, "CellType_Statistics.csv"), row.names = FALSE)

  # 打印分析摘要
  message("\n=== 分析完成摘要 ===")
  message(sprintf("总细胞数: %d", ncol(seurat_obj)))
  message(sprintf("总基因数: %d", nrow(seurat_obj)))
  message(sprintf("聚类数: %d", length(levels(seurat_obj$seurat_clusters))))
  message(sprintf("细胞类型数: %d", length(unique(seurat_obj$cell_type))))
  message(sprintf("实验分组: %s", paste(group_names, collapse = ", ")))

  # 打印各组细胞数量
  message("\n=== 各组细胞数量 ===")
  total_counts <- colSums(cell_type_counts)
  for (group in group_names) {
    message(sprintf("%s组: %d 个细胞", group, total_counts[group]))
  }

  # 打印细胞类型分布
  message("\n=== 细胞类型分布 ===")
  for (i in 1:nrow(cell_type_stats)) {
    celltype <- cell_type_stats$CellType[i]
    total <- cell_type_stats$Total_Count[i]
    group_info <- c()
    for (group in group_names) {
      count_col <- paste0(group, "_Count")
      pct_col <- paste0(group, "_Percentage")
      count <- cell_type_stats[[count_col]][i]
      pct <- cell_type_stats[[pct_col]][i]
      group_info <- c(group_info, sprintf("%s: %d (%.1f%%)", group, count, pct))
    }
    message(sprintf("%s (总计: %d) - %s", celltype, total, paste(group_info, collapse = ", ")))
  }

  return(cell_type_stats)
}

stats_report <- generate_statistics_report(Peripheral_Blood_Mononuclear_Cells)

# 保存Seurat对象
save(Peripheral_Blood_Mononuclear_Cells, cellAnn, file = file.path(output_dirs$basic_analysis, "SingleCell_Manual_Annotation.RData"))

# ==============================================================================
# 08. 分面火山图绘制（使用 jjVolcano - 参考官方方法）
# ==============================================================================

message("\n开始生成分面火山图（jjVolcano）...\n")

# 检查是否安装scRNAtoolVis包
if (!require("scRNAtoolVis", quietly = TRUE)) {
  message("正在安装scRNAtoolVis包...")
  tryCatch({
    install.packages("scRNAtoolVis", repos = "http://cran.us.r-project.org")
    library(scRNAtoolVis)
    message("✓ scRNAtoolVis包安装成功")
  }, error = function(e) {
    message(sprintf("✗ scRNAtoolVis安装失败: %s", e$message))
  })
}

# 获取分组信息
Type <- gsub("(.*?)\\..*", "\\1", colnames(Peripheral_Blood_Mononuclear_Cells))
names(Type) <- colnames(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
  object = Peripheral_Blood_Mononuclear_Cells,
  metadata = Type,
  col.name = "group"
)

message(sprintf("检测到的分组: %s", paste(unique(Type), collapse = ", ")))

# 获取细胞类型信息
cell_types_from_clusters <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
names(cell_types_from_clusters) <- colnames(Peripheral_Blood_Mononuclear_Cells)

# 创建组合标识符（分组_细胞类型）
combined_groups <- paste0(Type, "_", cell_types_from_clusters)
names(combined_groups) <- colnames(Peripheral_Blood_Mononuclear_Cells)
Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
  object = Peripheral_Blood_Mononuclear_Cells,
  metadata = combined_groups,
  col.name = "combined_group"
)

# 获取唯一的细胞类型
unique_cell_types <- unique(cell_types_from_clusters)
message(sprintf("检测到的细胞类型: %s", paste(unique_cell_types, collapse = ", ")))

# 提取分组名称
group_names <- unique(Type)

# 调试信息
message(sprintf("元数据检查:"))
message(sprintf("  - combined_group 列已添加: %s", "combined_group" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data)))
message(sprintf("  - unique combined_group 值数: %d", length(unique(Peripheral_Blood_Mononuclear_Cells$combined_group))))
message(sprintf("  - 示例 combined_group 值: %s", paste(head(unique(Peripheral_Blood_Mononuclear_Cells$combined_group), 3), collapse = ", ")))
if (length(group_names) < 2) {
  message("警告: 分组数少于2个，无法进行火山图分析")
} else {
  message(sprintf("分析分组: %s", paste(group_names, collapse = " vs ")))

  # 准备火山图数据
  volcano_data_list <- list()

  message("\n正在进行差异基因分析...\n")

  for (cellName in unique_cell_types) {
    # 清理细胞类型名称
    cellName_clean <- gsub("[^A-Za-z0-9_]", "_", as.character(cellName))

    conName <- paste0(group_names[1], "_", cellName)
    treatName <- paste0(group_names[2], "_", cellName)

    # 检查该细胞类型是否在两组中都存在
    con_exists <- any(Peripheral_Blood_Mononuclear_Cells$combined_group == conName)
    treat_exists <- any(Peripheral_Blood_Mononuclear_Cells$combined_group == treatName)

    if (!con_exists | !treat_exists) {
      message(sprintf("  ✗ 跳过 %s: 缺失分组(对照:%s, 处理:%s)", cellName, con_exists, treat_exists))
      next
    }

    # 获取该细胞类型在两组中的细胞
    conCells <- which(Peripheral_Blood_Mononuclear_Cells$combined_group == conName)
    treatCells <- which(Peripheral_Blood_Mononuclear_Cells$combined_group == treatName)

    # 转换为细胞名称
    if (length(conCells) > 0) {
      conCells <- colnames(Peripheral_Blood_Mononuclear_Cells)[conCells]
    }
    if (length(treatCells) > 0) {
      treatCells <- colnames(Peripheral_Blood_Mononuclear_Cells)[treatCells]
    }

    # 确保每组至少有3个细胞
    if (length(conCells) >= 3 & length(treatCells) >= 3) {
      message(sprintf("  分析: %s (对照: %d 细胞, 处理: %d 细胞)", cellName, length(conCells), length(treatCells)))

      # 进行差异表达分析
      tryCatch({
        pbmc.markers <- FindMarkers(
          Peripheral_Blood_Mononuclear_Cells,
          ident.1 = treatCells,
          ident.2 = conCells,
          logfc.threshold = 0,
          min.pct = 0.1,
          verbose = FALSE
        )

        # 添加基因名和细胞类型信息
        pbmc.markers$gene <- rownames(pbmc.markers)
        pbmc.markers$cluster <- cellName

        volcano_data_list[[as.character(cellName)]] <- pbmc.markers

        message(sprintf("    ✓ 找到 %d 个基因", nrow(pbmc.markers)))

      }, error = function(e) {
        message(sprintf("    ✗ 分析出错: %s", e$message))
      })

    } else {
      message(sprintf("  ✗ 跳过 %s: 细胞数不足 (对照: %d, 处理: %d)", cellName, length(conCells), length(treatCells)))
    }
  }

  message(sprintf("\n共 %d 个细胞类型有分析数据", length(volcano_data_list)))

  # 如果有数据，生成火山图
  if (length(volcano_data_list) > 0) {

    # 合并所有数据
    volcano_data <- do.call(rbind, volcano_data_list)
    message(sprintf("合并数据: %d 个基因", nrow(volcano_data)))

    # 过滤数据：移除以数字结尾和以LINC开头的基因
    volcano_data <- volcano_data %>%
      filter(
        !str_detect(gene, "\\.[0-9]$"),     # 过滤掉以 .数字 结尾的
        !str_detect(gene, "^LINC")          # 过滤掉以 LINC 开头的
      )

    message(sprintf("过滤后数据: %d 个基因", nrow(volcano_data)))

    # 获取唯一的细胞类型
    unique_clusters <- unique(volcano_data$cluster)
    message(sprintf("唯一细胞类型: %s", paste(unique_clusters, collapse = ", ")))

    # 生成颜色配置
    nColors <- length(unique_clusters)
    myColors <- colorRampPalette(brewer.pal(min(12, nColors), "Set3"))(nColors)
    names(myColors) <- unique_clusters

    # 确保cluster列是因子类型
    volcano_data$cluster <- factor(volcano_data$cluster, levels = unique_clusters)

    message("\n正在绘制jjVolcano火山图...\n")

    # 创建输出目录
    volcano_output_dir <- file.path(output_dirs$statistics_reports, "Volcano_Plot_Analysis")
    if (!dir.exists(volcano_output_dir)) {
      dir.create(volcano_output_dir, recursive = TRUE)
    }

    # 绘制火山图
    tryCatch({

      p <- jjVolcano(
        diffData = volcano_data,
        log2FC.cutoff = 1,
        pvalue.cutoff = 0.05,
        tile.col = myColors,
        size = 4.5,
        topGeneN = 5,
        legend.position = "top"
      ) +
        labs(
          title = "Volcano Plot of Differentially Expressed Genes",
          subtitle = sprintf("%s vs %s", group_names[2], group_names[1]),
          x = "Cell Type",
          y = "-log10(p-value)"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = ggplot2::margin(b = 10)),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50", margin = ggplot2::margin(b = 15)),
          legend.position = "top",
          legend.title = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 13, face = "bold")
        )

      message("✓ 火山图创建成功，正在保存...\n")

      # 保存为PDF
      pdf_path <- file.path(volcano_output_dir, "faceted_volcano_plot.pdf")
      ggsave(
        p,
        filename = pdf_path,
        width = 14,
        height = 9,
        dpi = 600
      )

      message(sprintf("✓ PDF已保存: faceted_volcano_plot.pdf\n"))

      # 保存分析数据
      csv_path <- file.path(volcano_output_dir, "volcano_data_for_plot.csv")
      write.csv(
        volcano_data,
        file = csv_path,
        row.names = FALSE
      )

      message(sprintf("✓ 分析数据已保存: volcano_data_for_plot.csv\n"))

    }, error = function(e) {
      message(sprintf("✗ 绘制火山图时出错: %s\n", e$message))
      print(e)
    })

  } else {
    message("\n✗ 警告: 没有足够的差异基因分析数据来绘制火山图！\n")
  }
}

message("✓ 分面火山图绘制完成！\n")

message("\n所有分析已完成！结果文件已保存到工作目录。")

# ==============================================================================
# 15. 高级细胞类型标记基因气泡图（Advanced Cell Type Marker Gene DotPlot）
# ==============================================================================

message("\n开始生成高级细胞类型标记基因气泡图...\n")

# 定义各细胞类型的经典标记基因（高级版本）
cell_type_markers <- list(
  "T cells" = c("CD3D", "CD3E", "CD4", "CD8A", "IL7R", "TRAC"),
  "NK cells" = c("GNLY", "NKG7", "FGFBP2", "GZMB", "GZMA", "PRF1"),
  "B cells" = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "IGHD"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9", "CST3", "FCGR1A"),
  "Dendritic Cells" = c("FCER1A", "CD1C", "DLA", "CLEC9A", "IDO1", "THBD"),
  "Plasma cells" = c("MZB1", "IGHG1", "IGHA1", "IGHM", "POU2AF1", "XBP1")
)

# 获取数据中实际存在的基因
valid_markers <- list()
for (celltype in names(cell_type_markers)) {
  valid_genes <- intersect(cell_type_markers[[celltype]], rownames(Peripheral_Blood_Mononuclear_Cells))
  if (length(valid_genes) > 0) {
    valid_markers[[celltype]] <- valid_genes
  }
}

message(sprintf("找到 %d 个细胞类型的有效标记基因", length(valid_markers)))

# 如果有有效的标记基因，绘制高级气泡图
if (length(valid_markers) > 0) {

  # 创建高级气泡图函数
  create_advanced_dotplot <- function(seurat_obj, markers_list, output_dir) {

    # 将标记基因列表转换为向量（保持顺序）
    all_markers <- unique(unlist(markers_list))

    message(sprintf("准备绘制包含 %d 个基因的气泡图...", length(all_markers)))

    # 获取细胞类型信息
    cell_types <- as.character(Idents(seurat_obj))

    # 创建DotPlot
    p <- DotPlot(
      seurat_obj,
      features = all_markers,
      group.by = "cell_type",
      dot.scale = 8,
      scale = TRUE,
      cols = c("lightgrey", "#002868", "#BF0A30")  # 蓝色到红色的渐变
    ) +

      # 主题美化
      theme_minimal() +
      theme(
        # 标题设置
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 18,
          margin = ggplot2::margin(b = 15)
        ),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 12,
          color = "gray40",
          margin = ggplot2::margin(b = 20)
        ),

        # 坐标轴
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 11,
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 11,
          face = "bold"
        ),

        # 网格线
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),

        # 图例
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.6, "cm"),

        # 边距
        plot.margin = ggplot2::margin(20, 20, 20, 20, unit = "pt")
      ) +

      # 标签
      labs(
        title = "Cell Type Marker Gene Expression Profile",
        subtitle = "DotPlot showing marker gene expression across annotated cell types",
        x = "Marker Genes",
        y = "Cell Types",
        size = "Percent Expressed",
        color = "Average Expression"
      ) +

      # 调整图例
      guides(
        color = guide_colorbar(
          title = "Avg Expr",
          title.position = "top",
          barwidth = 10,
          barheight = 1
        ),
        size = guide_legend(
          title = "% Expr",
          title.position = "top",
          override.aes = list(color = "#002868")
        )
      )

    return(p)
  }

  # 创建增强版气泡图（包含分组信息）
  message("\n正在生成分组对比高级气泡图...\n")

  tryCatch({

    # 计算每个细胞类型-基因组合在不同分组中的表达统计
    all_markers <- unique(unlist(valid_markers))

    # 获取数据
    expr_data <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, slot = "data")

    # 创建一个增强的数据框用于绘图
    plot_data <- data.frame()

    for (gene in all_markers) {
      if (gene %in% rownames(expr_data)) {
        for (celltype in unique(Idents(Peripheral_Blood_Mononuclear_Cells))) {
          celltype_cells <- which(Idents(Peripheral_Blood_Mononuclear_Cells) == celltype)

          if (length(celltype_cells) > 0) {
            gene_expr <- expr_data[gene, celltype_cells]

            plot_data <- rbind(plot_data, data.frame(
              gene = gene,
              celltype = celltype,
              avg_expr = mean(gene_expr),
              pct_expr = sum(gene_expr > 0) / length(gene_expr) * 100,
              cell_count = length(celltype_cells)
            ))
          }
        }
      }
    }

    # 生成增强气泡图
    p_enhanced <- ggplot(plot_data, aes(x = celltype, y = gene)) +
      geom_point(aes(size = pct_expr, color = avg_expr), alpha = 0.8) +

      # 颜色配置
      scale_color_gradient(
        low = "#F0F0F0",
        high = "#C41E3A",
        name = "Avg Expression"
      ) +
      scale_size_continuous(
        range = c(2, 10),
        name = "% Expressing Cells"
      ) +

      # 主题优化
      theme_minimal() +
      theme(
        # 标题
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 18,
          margin = ggplot2::margin(b = 10)
        ),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 12,
          color = "gray50",
          margin = ggplot2::margin(b = 15)
        ),

        # 坐标轴
        axis.title = element_text(size = 13, face = "bold", color = "black"),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 11,
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 9,
          face = "bold"
        ),

        # 网格
        panel.grid.major.x = element_line(color = "gray95", linewidth = 0.2),
        panel.grid.major.y = element_line(color = "gray95", linewidth = 0.2),
        panel.grid.minor = element_blank(),

        # 背景
        panel.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),

        # 图例
        legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.5, "cm"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),

        # 边距
        plot.margin = ggplot2::margin(20, 20, 20, 20, unit = "pt")
      ) +

      labs(
        title = "Cell Type & Marker Gene Correlation Analysis",
        subtitle = "Expression patterns revealing cell type identity",
        x = "Cell Types",
        y = "Marker Genes"
      )

    # 保存增强版气泡图
    pdf_path2 <- file.path(output_dirs$cell_annotation, "Enhanced_CellType_Gene_Correlation_DotPlot.pdf")
    pdf(pdf_path2, width = 14, height = 9)
    print(p_enhanced)
    dev.off()

    message(sprintf("✓ 增强气泡图已保存: Enhanced_CellType_Gene_Correlation_DotPlot.pdf"))

  }, error = function(e) {
    message(sprintf("✗ 生成增强气泡图时出错: %s", e$message))
  })

} else {
  message("警告: 未找到有效的标记基因，跳过高级气泡图生成")
}

message("\n✓ 高级气泡图生成完成！\n")

