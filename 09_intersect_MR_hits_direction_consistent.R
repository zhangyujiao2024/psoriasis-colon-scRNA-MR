# ==========================================
# MR随机化结果取方向一致的交集
# ==========================================
# 脚本说明：
# 1. 读取目录下所有CSV文件的MR分析结果
# 2. 根据OR值判断效应方向（OR>1为正向，OR<1为负向）
# 3. 取各文件中方向一致的基因交集
# 4. 输出交集结果和Venn图
# ==========================================

# --- 加载必要包 ---
library(ggvenn)
library(ggplot2)

# --- 设置工作目录 ---
setwd("D:\\生信分析狮\\22.3种细胞")

# --- 创建输出文件夹 ---
outputDir <- "output_folder"
if (!dir.exists(outputDir)) dir.create(outputDir)

# --- 获取所有CSV文件 ---
files_csv <- list.files(pattern = "\\.csv$", full.names = FALSE)
files_csv <- files_csv[!grepl("^output_|^Intersection", files_csv)]  # 排除输出文件

cat("==========================================\n")
cat("MR Results Intersection Analysis\n")
cat("==========================================\n")
cat(sprintf("Found %d CSV files to process:\n", length(files_csv)))
for (f in files_csv) {
  cat(sprintf("  - %s\n", f))
}
cat("\n")

# --- 读取所有MR结果文件，按方向分类基因 ---
geneList_positive <- list()  # OR > 1 (正向效应)
geneList_negative <- list()  # OR < 1 (负向效应)
geneData_positive <- list()  # 存储基因和OR值
geneData_negative <- list()  # 存储基因和OR值

for (inputFile in files_csv) {
  cat(sprintf("Processing: %s\n", inputFile))

  # 读取CSV文件
  rt <- read.csv(inputFile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  # 检查是否有exposure列（基因名）和or列（OR值）
  if (!"exposure" %in% colnames(rt)) {
    cat(sprintf("  Warning: No 'exposure' column found, skipping...\n"))
    next
  }

  if (!"or" %in% colnames(rt)) {
    cat(sprintf("  Warning: No 'or' column found, skipping...\n"))
    next
  }

  # 提取基因名和OR值
  geneNames <- rt$exposure
  orValues <- rt$or

  # 移除NA值
  valid_idx <- !is.na(geneNames) & !is.na(orValues) & geneNames != ""
  geneNames <- geneNames[valid_idx]
  orValues <- orValues[valid_idx]

  # 根据OR值分类
  # OR > 1: 正向效应（增加风险）
  positive_idx <- orValues > 1
  positiveGenes <- geneNames[positive_idx]
  positiveOR <- orValues[positive_idx]

  # OR < 1: 负向效应（降低风险）
  negative_idx <- orValues < 1
  negativeGenes <- geneNames[negative_idx]
  negativeOR <- orValues[negative_idx]

  # 文件名作为集合名称
  setName <- tools::file_path_sans_ext(basename(inputFile))

  # 存储正向效应基因
  if (length(positiveGenes) > 0) {
    geneList_positive[[setName]] <- unique(positiveGenes)
    geneData_positive[[setName]] <- data.frame(
      Gene = positiveGenes,
      OR = positiveOR,
      stringsAsFactors = FALSE
    )
  }

  # 存储负向效应基因
  if (length(negativeGenes) > 0) {
    geneList_negative[[setName]] <- unique(negativeGenes)
    geneData_negative[[setName]] <- data.frame(
      Gene = negativeGenes,
      OR = negativeOR,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("  Positive effect genes (OR>1): %d\n", length(unique(positiveGenes))))
  cat(sprintf("  Negative effect genes (OR<1): %d\n", length(unique(negativeGenes))))
}

cat("\n==========================================\n")
cat("Computing Intersections\n")
cat("==========================================\n")

# --- 计算正向效应基因的交集 ---
if (length(geneList_positive) > 1) {
  intersectionGenes_positive <- Reduce(intersect, geneList_positive)
} else {
  intersectionGenes_positive <- c()
}

# --- 计算负向效应基因的交集 ---
if (length(geneList_negative) > 1) {
  intersectionGenes_negative <- Reduce(intersect, geneList_negative)
} else {
  intersectionGenes_negative <- c()
}

intersectionCount_positive <- length(intersectionGenes_positive)
intersectionCount_negative <- length(intersectionGenes_negative)

cat(sprintf("Common positive effect genes (OR>1): %d\n", intersectionCount_positive))
cat(sprintf("Common negative effect genes (OR<1): %d\n", intersectionCount_negative))

# --- 绘制正向效应基因Venn图 ---
if (length(geneList_positive) > 1) {
  pdf(file = file.path(outputDir, "venn_PositiveEffect.pdf"), width = 6, height = 6)
  venn_plot_positive <- ggvenn(
    geneList_positive,
    show_percentage = TRUE,
    stroke_color = "white",
    stroke_size = 0.5,
    fill_color = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")[seq_along(geneList_positive)],
    set_name_color = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8")[seq_along(geneList_positive)],
    set_name_size = 2.5,
    text_size = 2
  ) + ggtitle("Positive Effect Genes (OR > 1)")
  print(venn_plot_positive)
  dev.off()
  cat("Saved: venn_PositiveEffect.pdf\n")
}

# --- 绘制负向效应基因Venn图 ---
if (length(geneList_negative) > 1) {
  pdf(file = file.path(outputDir, "venn_NegativeEffect.pdf"), width = 6, height = 6)
  venn_plot_negative <- ggvenn(
    geneList_negative,
    show_percentage = TRUE,
    stroke_color = "white",
    stroke_size = 0.5,
    fill_color = c("#6C5CE7", "#00B894", "#FDCB6E", "#E17055", "#74B9FF")[seq_along(geneList_negative)],
    set_name_color = c("#6C5CE7", "#00B894", "#FDCB6E", "#E17055", "#74B9FF")[seq_along(geneList_negative)],
    set_name_size = 2.5,
    text_size = 2
  ) + ggtitle("Negative Effect Genes (OR < 1)")
  print(venn_plot_negative)
  dev.off()
  cat("Saved: venn_NegativeEffect.pdf\n")
}

# --- 输出交集基因的详细信息 ---
if (length(intersectionGenes_positive) > 0 || length(intersectionGenes_negative) > 0) {

  # 处理正向效应交集基因
  if (length(intersectionGenes_positive) > 0) {
    intersection_data_positive <- data.frame(
      Gene = intersectionGenes_positive,
      Direction = "Positive (OR>1)",
      stringsAsFactors = FALSE
    )

    # 添加每个数据集的OR值
    for (setName in names(geneData_positive)) {
      setData <- geneData_positive[[setName]]
      or_values <- c()
      for (gene in intersectionGenes_positive) {
        idx <- which(setData$Gene == gene)
        if (length(idx) > 0) {
          or_values <- c(or_values, setData$OR[idx[1]])
        } else {
          or_values <- c(or_values, NA)
        }
      }
      intersection_data_positive[[paste0(setName, "_OR")]] <- or_values
    }
  } else {
    intersection_data_positive <- data.frame()
  }

  # 处理负向效应交集基因
  if (length(intersectionGenes_negative) > 0) {
    intersection_data_negative <- data.frame(
      Gene = intersectionGenes_negative,
      Direction = "Negative (OR<1)",
      stringsAsFactors = FALSE
    )

    # 添加每个数据集的OR值
    for (setName in names(geneData_negative)) {
      setData <- geneData_negative[[setName]]
      or_values <- c()
      for (gene in intersectionGenes_negative) {
        idx <- which(setData$Gene == gene)
        if (length(idx) > 0) {
          or_values <- c(or_values, setData$OR[idx[1]])
        } else {
          or_values <- c(or_values, NA)
        }
      }
      intersection_data_negative[[paste0(setName, "_OR")]] <- or_values
    }
  } else {
    intersection_data_negative <- data.frame()
  }

  # 合并正向和负向数据
  if (nrow(intersection_data_positive) > 0 && nrow(intersection_data_negative) > 0) {
    intersection_data_combined <- rbind(intersection_data_positive, intersection_data_negative)
  } else if (nrow(intersection_data_positive) > 0) {
    intersection_data_combined <- intersection_data_positive
  } else {
    intersection_data_combined <- intersection_data_negative
  }

  # 输出合并的CSV文件
  if (nrow(intersection_data_combined) > 0) {
    write.csv(intersection_data_combined,
              file = file.path(outputDir, "Intersection_MR_hits_directional.csv"),
              row.names = FALSE)
    cat("Saved: Intersection_MR_hits_directional.csv\n")
  }
}

# --- 输出正向效应交集基因列表 ---
if (length(intersectionGenes_positive) > 0) {
  write.table(intersectionGenes_positive,
              file = file.path(outputDir, "IntersectionGenes_PositiveEffect.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  cat("Saved: IntersectionGenes_PositiveEffect.txt\n")
}

# --- 输出负向效应交集基因列表 ---
if (length(intersectionGenes_negative) > 0) {
  write.table(intersectionGenes_negative,
              file = file.path(outputDir, "IntersectionGenes_NegativeEffect.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  cat("Saved: IntersectionGenes_NegativeEffect.txt\n")
}

# --- 绘制柱状图 ---
if (length(geneList_positive) > 1) {
  setSizes_positive <- sapply(geneList_positive, length)
  pdf(file = file.path(outputDir, "barplot_PositiveEffect.pdf"), width = 6, height = 4)
  print(ggplot(data.frame(Set = names(setSizes_positive), Size = setSizes_positive),
               aes(x = Set, y = Size, fill = Set)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylab("Gene Count") +
    xlab("") +
    ggtitle("Positive Effect Genes (OR > 1)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.position = "none"))
  dev.off()
  cat("Saved: barplot_PositiveEffect.pdf\n")
}

if (length(geneList_negative) > 1) {
  setSizes_negative <- sapply(geneList_negative, length)
  pdf(file = file.path(outputDir, "barplot_NegativeEffect.pdf"), width = 6, height = 4)
  print(ggplot(data.frame(Set = names(setSizes_negative), Size = setSizes_negative),
               aes(x = Set, y = Size, fill = Set)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylab("Gene Count") +
    xlab("") +
    ggtitle("Negative Effect Genes (OR < 1)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.position = "none"))
  dev.off()
  cat("Saved: barplot_NegativeEffect.pdf\n")
}

cat("\n==========================================\n")
cat("Analysis Complete!\n")
cat("==========================================\n")
cat(sprintf("Common positive effect genes (OR>1): %d\n", intersectionCount_positive))
cat(sprintf("Common negative effect genes (OR<1): %d\n", intersectionCount_negative))
cat("==========================================\n")

