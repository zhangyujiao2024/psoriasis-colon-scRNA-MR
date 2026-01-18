# ==============================================================================
# 单细胞分析 - 基因在注释细胞中的表达差异分析
# ==============================================================================

# 加载必需的R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(rstatix)
library(gridExtra)
library(grid)  # 用于textGrob和gpar函数

message("✓ 所有必需的R包已成功加载完成！")

# ==============================================================================
# 参数设置
# ==============================================================================

# 设置工作目录
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\30.研究的基因在注释细胞的表达情况 -第2种疾病"
setwd(workDir)
message("工作目录设置为：", getwd())

# 设置分析的基因
genes_to_analyze <- c("ALDH2", "ARL4C","LY9", "TYMP","FCER1G", "LYZ")
message(sprintf("待分析基因: %s", paste(genes_to_analyze, collapse = ", ")))
# 创建输出目录
output_dir <- "Gene_Expression_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(sprintf("创建输出目录: %s", output_dir))
}

# ==============================================================================
# 读取数据
# ==============================================================================

message("\n读取单细胞数据和细胞类型注释...")

# 读取Seurat对象
seurat_obj <- readRDS("filtered_single_cell_data.rds")

# 创建Seurat对象（如果还没有）
if (!inherits(seurat_obj, "Seurat")) {
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, min.cells = 3, min.features = 200)
  # 进行基本的数据归一化
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
}

message(sprintf("数据包含 %d 个基因和 %d 个细胞", nrow(seurat_obj), ncol(seurat_obj)))

# 读取细胞类型注释文件
cell_mapping <- read.csv("Cell_ID_to_CellType_Mapping.csv", stringsAsFactors = FALSE)
message(sprintf("读取细胞注释文件: %d 个细胞", nrow(cell_mapping)))

# ==============================================================================
# 添加细胞类型注释到Seurat对象
# ==============================================================================

message("\n将细胞类型注释添加到Seurat对象...")

# 智能细胞ID匹配 - 修复版本
message("开始智能匹配细胞ID...")

# 显示原始Cell_ID格式示例
message("Cell_ID映射文件格式示例:")
message(paste("  ", head(cell_mapping$Cell_ID, 3), collapse = "\n"))

# 显示Seurat对象中细胞名称格式示例
message("Seurat对象中细胞名称格式示例:")
message(paste("  ", head(colnames(seurat_obj), 3), collapse = "\n"))

# 方法1：直接匹配（原始格式）
common_cells_direct <- intersect(cell_mapping$Cell_ID, colnames(seurat_obj))
message(sprintf("方法1 - 直接匹配: %d 个细胞", length(common_cells_direct)))

# 方法2：统一分隔符后匹配
cell_mapping$Cell_ID_dot <- gsub("-", ".", cell_mapping$Cell_ID)
common_cells_dot <- intersect(cell_mapping$Cell_ID_dot, colnames(seurat_obj))
message(sprintf("方法2 - 统一分隔符: %d 个细胞", length(common_cells_dot)))

# 方法3：提取细胞条码部分（移除前缀）
# 从 "Control.GSM5102900_AAACCCAGTTGAGAGC-1" 提取 "AAACCCAGTTGAGAGC-1" 或 "AAACCCAGTTGAGAGC.1"
cell_mapping$Cell_ID_barcode <- gsub(".*_([ATCG]+[-\\.]\\d+)$", "\\1", cell_mapping$Cell_ID)
common_cells_barcode <- intersect(cell_mapping$Cell_ID_barcode, colnames(seurat_obj))
message(sprintf("方法3 - 细胞条码匹配: %d 个细胞", length(common_cells_barcode)))

# 方法4：细胞条码+分隔符转换
cell_mapping$Cell_ID_barcode_dot <- gsub("-", ".", cell_mapping$Cell_ID_barcode)
common_cells_barcode_dot <- intersect(cell_mapping$Cell_ID_barcode_dot, colnames(seurat_obj))
message(sprintf("方法4 - 细胞条码+分隔符转换: %d 个细胞", length(common_cells_barcode_dot)))

# 方法5：尝试更灵活的匹配模式
# 提取 GSM部分 + 细胞条码
cell_mapping$Cell_ID_gsm <- gsub(".*\\.(GSM\\d+_[ATCG]+[-\\.]\\d+)$", "\\1", cell_mapping$Cell_ID)
cell_mapping$Cell_ID_gsm_dot <- gsub("-", ".", cell_mapping$Cell_ID_gsm)
common_cells_gsm <- intersect(cell_mapping$Cell_ID_gsm, colnames(seurat_obj))
common_cells_gsm_dot <- intersect(cell_mapping$Cell_ID_gsm_dot, colnames(seurat_obj))
message(sprintf("方法5a - GSM+条码匹配: %d 个细胞", length(common_cells_gsm)))
message(sprintf("方法5b - GSM+条码+分隔符转换: %d 个细胞", length(common_cells_gsm_dot)))

# 选择最佳匹配方法
match_results <- list(
  direct = common_cells_direct,
  dot = common_cells_dot,
  barcode = common_cells_barcode,
  barcode_dot = common_cells_barcode_dot,
  gsm = common_cells_gsm,
  gsm_dot = common_cells_gsm_dot
)

match_counts <- sapply(match_results, length)
best_method <- names(match_counts)[which.max(match_counts)]
best_count <- max(match_counts)

message(sprintf("🎯 选择最佳匹配方法: %s (%d 个细胞)", best_method, best_count))

# 应用最佳匹配方法
if (best_count == 0) {
  stop("❌ 所有匹配方法都失败了，请检查Cell_ID格式！

可能的解决方案:
1. 检查Cell_ID_to_CellType_Mapping.csv文件格式是否正确
2. 检查Seurat对象是否为正确的数据文件
3. 确认Cell_ID格式与Seurat对象中的细胞名称格式是否一致")
}

common_cells <- match_results[[best_method]]

# 更新cell_mapping的Cell_ID列为匹配成功的格式
if (best_method == "direct") {
  # 保持原格式
} else if (best_method == "dot") {
  cell_mapping$Cell_ID <- cell_mapping$Cell_ID_dot
} else if (best_method == "barcode") {
  cell_mapping$Cell_ID <- cell_mapping$Cell_ID_barcode
} else if (best_method == "barcode_dot") {
  cell_mapping$Cell_ID <- cell_mapping$Cell_ID_barcode_dot
} else if (best_method == "gsm") {
  cell_mapping$Cell_ID <- cell_mapping$Cell_ID_gsm
} else if (best_method == "gsm_dot") {
  cell_mapping$Cell_ID <- cell_mapping$Cell_ID_gsm_dot
}

message(sprintf("✅ 成功匹配到 %d 个细胞 (总计: %d, 匹配率: %.1f%%)",
                length(common_cells), ncol(seurat_obj),
                length(common_cells)/ncol(seurat_obj)*100))

# 只保留匹配的细胞
seurat_obj <- seurat_obj[, common_cells]
cell_mapping <- cell_mapping[cell_mapping$Cell_ID %in% common_cells, ]

# 按照Seurat对象中的细胞顺序排列注释
cell_mapping <- cell_mapping[match(colnames(seurat_obj), cell_mapping$Cell_ID), ]

# 添加细胞类型到Seurat对象
seurat_obj$cell_type <- cell_mapping$Cell_Type

# 提取分组信息（Control vs T2DM）
Type <- gsub("(.*?)\\..*", "\\1", colnames(seurat_obj))
seurat_obj$Group <- Type

message(sprintf("检测到的分组: %s", paste(unique(Type), collapse = ", ")))
message(sprintf("检测到的细胞类型: %s", paste(unique(seurat_obj$cell_type), collapse = ", ")))

# ==============================================================================
# 基因表达数据提取和统计分析
# ==============================================================================

message("\n提取基因表达数据并进行统计分析...")

# 创建存储结果的列表
expression_data_list <- list()
statistical_results <- list()

for (gene in genes_to_analyze) {

  message(sprintf("\n分析基因: %s", gene))

  # 检查基因是否存在
  if (!gene %in% rownames(seurat_obj)) {
    message(sprintf("  警告: 基因 %s 不在数据中，跳过", gene))
    next
  }

  # 提取基因表达数据
  gene_expr <- GetAssayData(seurat_obj, slot = "data")[gene, ]

  # 创建数据框
  expr_df <- data.frame(
    Cell_ID = colnames(seurat_obj),
    Expression = as.numeric(gene_expr),
    Group = seurat_obj$Group,
    Cell_Type = seurat_obj$cell_type,
    stringsAsFactors = FALSE
  )

  # 移除缺失值
  expr_df <- expr_df[!is.na(expr_df$Cell_Type) & expr_df$Cell_Type != "", ]

  message(sprintf("  有效细胞数: %d", nrow(expr_df)))

  # 保存表达数据
  expression_data_list[[gene]] <- expr_df

  # 统计分析：对每个细胞类型进行组间比较
  stat_results <- data.frame()

  for (celltype in unique(expr_df$Cell_Type)) {

    # 该细胞类型的数据
    celltype_data <- expr_df[expr_df$Cell_Type == celltype, ]

    # 检查两组是否都有数据
    groups_present <- unique(celltype_data$Group)
    if (length(groups_present) < 2) {
      message(sprintf("  %s: 只有一个分组，跳过统计", celltype))
      next
    }

    # 获取两组的细胞数
    group_counts <- table(celltype_data$Group)
    if (any(group_counts < 3)) {
      message(sprintf("  %s: 细胞数不足 (%s), 跳过",
                      celltype, paste(names(group_counts), group_counts, sep=":", collapse=", ")))
      next
    }

    # Wilcoxon秩和检验（适合非正态分布）
    test_result <- wilcox.test(
      Expression ~ Group,
      data = celltype_data,
      exact = FALSE
    )

    # 计算统计量
    control_data <- celltype_data[celltype_data$Group == "Control", "Expression"]
    treat_data <- celltype_data[celltype_data$Group == "Treat", "Expression"]

    stat_results <- rbind(stat_results, data.frame(
      Gene = gene,
      Cell_Type = celltype,
      Control_Mean = mean(control_data, na.rm = TRUE),
      Treat_Mean = mean(treat_data, na.rm = TRUE),
      Control_Median = median(control_data, na.rm = TRUE),
      Treat_Median = median(treat_data, na.rm = TRUE),
      Control_N = length(control_data),
      Treat_N = length(treat_data),
      P_Value = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }

  # 多重检验校正
  if (nrow(stat_results) > 0) {
    stat_results$P_Adjusted <- p.adjust(stat_results$P_Value, method = "BH")

    # 添加显著性标记
    stat_results$Significance <- ifelse(stat_results$P_Adjusted < 0.001, "***",
                                        ifelse(stat_results$P_Adjusted < 0.01, "**",
                                               ifelse(stat_results$P_Adjusted < 0.05, "*", "ns")))

    message(sprintf("  完成统计分析: %d 个细胞类型", nrow(stat_results)))

    # 保存统计结果
    statistical_results[[gene]] <- stat_results

    # 保存为CSV
    write.csv(stat_results,
              file.path(output_dir, paste0(gene, "_statistical_results.csv")),
              row.names = FALSE)
  }
}

# ==============================================================================
# 生成小提琴图 - 与参考图相同的风格
# ==============================================================================

message("\n生成基因表达小提琴图...")

for (gene in genes_to_analyze) {

  if (!gene %in% names(expression_data_list)) {
    next
  }

  expr_df <- expression_data_list[[gene]]
  stat_df <- statistical_results[[gene]]

  if (is.null(stat_df) || nrow(stat_df) == 0) {
    message(sprintf("  跳过 %s: 无统计结果", gene))
    next
  }

  message(sprintf("\n绘制 %s 的表达图...", gene))

  # 自动检测实际的分组值
  unique_groups <- unique(expr_df$Group)
  unique_groups <- unique_groups[!is.na(unique_groups) & unique_groups != ""]
  message(sprintf("  检测到分组: %s", paste(unique_groups, collapse = " vs ")))

  # 智能配色 - 自动适配实际分组
  if (length(unique_groups) == 2) {
    # 经典红蓝配色
    color_palette <- c("#E74C3C", "#3498DB")
    names(color_palette) <- unique_groups
  } else {
    # 多分组时使用更多颜色
    color_palette <- RColorBrewer::brewer.pal(max(3, length(unique_groups)), "Set1")[1:length(unique_groups)]
    names(color_palette) <- unique_groups
  }

  # 获取有统计结果的细胞类型
  cell_types_with_stats <- stat_df$Cell_Type
  expr_df_filtered <- expr_df[expr_df$Cell_Type %in% cell_types_with_stats, ]

  # 设置因子水平
  expr_df_filtered$Cell_Type <- factor(expr_df_filtered$Cell_Type, levels = cell_types_with_stats)
  expr_df_filtered$Group <- factor(expr_df_filtered$Group, levels = unique_groups)

  # ==============================================================================
  # 第一部分：分别保存每个细胞类型的小提琴图（保持原来的输出格式）
  # ==============================================================================

  # 创建小提琴图列表
  plot_list <- list()

  for (i in 1:length(cell_types_with_stats)) {
    celltype <- cell_types_with_stats[i]

    # 该细胞类型的数据
    plot_data <- expr_df_filtered[expr_df_filtered$Cell_Type == celltype, ]

    if (nrow(plot_data) == 0) next

    # 获取统计显著性
    stat_row <- stat_df[stat_df$Cell_Type == celltype, ]
    if (nrow(stat_row) == 0) next

    sig <- stat_row$Significance[1]
    p_val <- stat_row$P_Adjusted[1]

    # 计算Y轴最大值（用于放置显著性标记）
    y_max <- max(plot_data$Expression, na.rm = TRUE)
    y_position <- y_max * 1.15

    # 创建小提琴图
    p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +

      # 小提琴图
      geom_violin(trim = TRUE, alpha = 0.7, color = "black", linewidth = 0.5) +

      # 箱线图叠加
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA,
                   color = "black", linewidth = 0.5) +

      # 智能配色
      scale_fill_manual(values = color_palette) +

      # 添加显著性标记
      {if (!is.na(sig) && sig != "ns") {
        annotate("text", x = mean(1:length(unique_groups)), y = y_position,
                 label = sig, size = 5, fontface = "bold")
      }} +

      # 主题设置
      theme_classic() +
      theme(
        # 标题
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),

        # 坐标轴
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),

        # 坐标轴线
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),

        # 图例
        legend.position = "none",

        # 背景
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),

        # 边距
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
      ) +

      # 标签
      labs(
        title = celltype,
        y = "Expression"
      ) +

      # Y轴扩展
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

    plot_list[[i]] <- p
  }

  # 组合所有小提琴图 - 保持原来的格式
  n_plots <- length(plot_list)
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)

  # 创建整体标题
  title_grob <- grid::textGrob(
    paste(gene, "\nExpression in different cell types"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  # 保存PDF - 原来的格式
  pdf_path <- file.path(output_dir, paste0(gene, "_violin_plot.pdf"))
  pdf(pdf_path, width = n_cols * 3, height = n_rows * 2.5 + 1.5)

  grid.arrange(
    title_grob,
    grobs = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    top = grid::textGrob("")
  )

  dev.off()

  message(sprintf("  ✅ 保存PDF（多面板）: %s", basename(pdf_path)))

  # ==============================================================================
  # 第二部分：生成组合单图（所有细胞类型在一个图中）- Nature/Reference风格
  # ==============================================================================

  message(sprintf("  生成组合单图版本（所有细胞类型）..."))

  # 创建单个组合图 - 所有细胞类型在一个图中
  n_celltypes <- length(cell_types_with_stats)

  if (n_celltypes <= 10) {  # 如果细胞类型不太多，可以放在一个图中

    # 为了在一个图中显示多个小提琴图，需要分面处理
    # 创建用于分面的图
    p_combined <- ggplot(expr_df_filtered, aes(x = Group, y = Expression, fill = Group)) +

      # 小提琴图
      geom_violin(trim = TRUE, alpha = 0.7, color = "black", linewidth = 0.4) +

      # 箱线图叠加
      geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA,
                   color = "black", linewidth = 0.4) +

      # 配色
      scale_fill_manual(values = color_palette, name = "Group") +

      # 分面 - 按细胞类型分（每个facet独立Y轴）
      facet_wrap(~Cell_Type, nrow = 1, ncol = ceiling(n_celltypes/1), scales = "free_y") +

      # 主题设置
      theme_classic() +
      theme(
        # 标题
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)),

        # 坐标轴标题
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),

        # 坐标轴文本
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),

        # 分面标签 - 增大字体
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_rect(color = "black", fill = "white", size = 0.5),

        # 坐标轴线
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),

        # 图例 - 增大字体
        legend.position = "top",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),

        # 背景
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.background = element_rect(fill = "white"),

        # 边距
        plot.margin = margin(15, 15, 15, 15)
      ) +

      # 标签
      labs(
        title = paste(gene, "- Expression across cell types"),
        x = "Group",
        y = "Expression Level"
      ) +

      # Y轴扩展
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

    # 添加统计显著性标记
    # 创建一个包含统计信息的数据框
    sig_data <- data.frame()
    for (i in 1:nrow(stat_df)) {
      celltype <- stat_df$Cell_Type[i]
      sig <- stat_df$Significance[i]

      # 计算该细胞类型数据的Y轴位置
      celltype_data <- expr_df_filtered[expr_df_filtered$Cell_Type == celltype, ]
      if (nrow(celltype_data) > 0) {
        y_max <- max(celltype_data$Expression, na.rm = TRUE)
        y_pos <- y_max * 1.1

        sig_data <- rbind(sig_data, data.frame(
          Cell_Type = celltype,
          x_pos = mean(1:length(unique_groups)),
          y_pos = y_pos,
          sig = sig,
          stringsAsFactors = FALSE
        ))
      }
    }

    # 添加显著性标记
    if (nrow(sig_data) > 0) {
      p_combined <- p_combined +
        geom_text(data = sig_data, aes(x = x_pos, y = y_pos, label = sig),
                  size = 6, fontface = "bold", color = "black",
                  inherit.aes = FALSE)
    }

    # 保存组合单图 - PDF
    combined_pdf_path <- file.path(output_dir, paste0(gene, "_combined_violin_plot.pdf"))
    ggsave(combined_pdf_path, plot = p_combined, width = 14, height = 5, dpi = 300, device = "pdf")
    message(sprintf("  ✅ 保存组合单图PDF: %s", basename(combined_pdf_path)))

  } else {
    message(sprintf("  提示：细胞类型太多（%d个），跳过组合单图", n_celltypes))
  }
}

# ==============================================================================
# 生成汇总表
# ==============================================================================

message("\n生成汇总统计表...")

# 合并所有基因的统计结果
if (length(statistical_results) > 0) {
  all_stats <- do.call(rbind, statistical_results)

  # 保存完整统计表
  write.csv(all_stats,
            file.path(output_dir, "All_Genes_Statistical_Summary.csv"),
            row.names = FALSE)

  # 生成简化摘要表
  summary_table <- all_stats %>%
    select(Gene, Cell_Type, Control_Mean, Treat_Mean,
           Control_N, Treat_N, P_Adjusted, Significance) %>%
    arrange(Gene, Cell_Type)

  write.csv(summary_table,
            file.path(output_dir, "Statistical_Summary_Simplified.csv"),
            row.names = FALSE)

  message(sprintf("保存汇总表: %d 行", nrow(all_stats)))
}

# ==============================================================================
# 输出分析结果摘要
# ==============================================================================

message("\n" , paste(rep("=", 60), collapse = ""))
message("分析完成！")
message(paste(rep("=", 60), collapse = ""))

message("\n生成的文件:")
message(sprintf("  输出目录: %s", output_dir))

for (gene in genes_to_analyze) {
  if (gene %in% names(statistical_results)) {
    n_celltypes <- nrow(statistical_results[[gene]])
    n_sig <- sum(statistical_results[[gene]]$P_Adjusted < 0.05, na.rm = TRUE)
    message(sprintf("  - %s: %d 个细胞类型, %d 个显著差异", gene, n_celltypes, n_sig))
  }
}

message("\n所有图表已保存到:", file.path(workDir, output_dir))
message(paste(rep("=", 60), collapse = ""))
