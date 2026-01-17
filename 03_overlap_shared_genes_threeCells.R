# --- 加载必要包 ---
library(ggvenn)
library(gridExtra)
library(ggplot2)
library(eulerr)

# --- 设置工作目录 ---
setwd("D:\\生信分析狮\\10.三种细胞")

# --- 创建输出文件夹 ---
outputDir <- "output_folder"
if (!dir.exists(outputDir)) dir.create(outputDir)

# --- 获取所有CSV和TXT文件（不包含compound.txt、IntersectionGenes.txt等输出文件）---
files_txt <- list.files(pattern = "\\.txt$")
files_csv <- list.files(pattern = "\\.csv$")
files <- c(files_txt, files_csv)
files <- files[!files %in% c("compound.txt", "IntersectionGenes.txt", "IntersectionGenes_UpRegulated.txt", "IntersectionGenes_DownRegulated.txt")]

# --- 读取所有文件中的基因信息，分别整理上调和下调基因 ---
geneList_up <- list()
geneList_down <- list()
geneData_up <- list()  # 存储基因和logFC
geneData_down <- list()  # 存储基因和logFC

for (inputFile in files) {
  # 判断文件类型并设置分隔符
  if (tolower(tools::file_ext(inputFile)) == "csv") {
    sep_char <- ","
    file_type <- "CSV"
  } else {
    sep_char <- "\t"
    file_type <- "TXT"
  }

  # 先读取第一行检查是否含有表头"gene"
  first_line <- readLines(inputFile, n = 1)
  has_header <- tolower(first_line) == "gene" || grepl("^gene$", first_line, ignore.case = TRUE) || grepl("Gene", first_line)

  # 根据是否有表头选择跳过第一行
  if (has_header) {
    rt <- read.table(inputFile, header = TRUE, sep = sep_char, check.names = FALSE, stringsAsFactors = FALSE)
    cat(sprintf("  %s (%s): Contains header, data loaded\n", inputFile, file_type))
  } else {
    rt <- read.table(inputFile, header = FALSE, sep = sep_char, check.names = FALSE, stringsAsFactors = FALSE)
    cat(sprintf("  %s (%s): No header, all rows included\n", inputFile, file_type))
  }

  # 检查是否有avg_log2FC列来区分上调和下调
  if ("avg_log2FC" %in% colnames(rt)) {
    # 有avg_log2FC列，分别提取上调和下调基因
    geneNames <- rt[, 1]
    log2fc <- rt[, "avg_log2FC"]

    # 上调基因 (avg_log2FC > 0)
    upGenes <- geneNames[log2fc > 0]
    upLogFC <- log2fc[log2fc > 0]
    downGenes <- geneNames[log2fc < 0]
    downLogFC <- log2fc[log2fc < 0]

    # 清理基因名
    upGenes <- trimws(upGenes)
    upGenes <- upGenes[upGenes != ""]
    downGenes <- trimws(downGenes)
    downGenes <- downGenes[downGenes != ""]

    setName <- tools::file_path_sans_ext(basename(inputFile))

    # 存储上调和下调基因
    geneList_up[[setName]] <- unique(upGenes)
    geneList_down[[setName]] <- unique(downGenes)

    # 存储基因和logFC的对应关系
    geneData_up[[setName]] <- data.frame(Gene = upGenes, logFC = upLogFC, stringsAsFactors = FALSE)
    geneData_down[[setName]] <- data.frame(Gene = downGenes, logFC = downLogFC, stringsAsFactors = FALSE)

    cat(sprintf("    Upregulated genes: %d\n", length(unique(upGenes))))
    cat(sprintf("    Downregulated genes: %d\n", length(unique(downGenes))))
  } else {
    # 没有avg_log2FC列，所有基因作为一个集合
    geneNames <- unlist(strsplit(as.vector(rt[, 1]), " "))
    geneNames <- trimws(geneNames)
    geneNames <- geneNames[geneNames != ""]
    setName <- tools::file_path_sans_ext(basename(inputFile))

    geneList_up[[setName]] <- unique(geneNames)
    geneData_up[[setName]] <- data.frame(Gene = geneNames, logFC = NA, stringsAsFactors = FALSE)
    cat(sprintf("    Total genes: %d\n", length(unique(geneNames))))
  }
}

# --- 计算共同上调和共同下调基因的交集 ---
if (length(geneList_up) > 0) {
  intersectionGenes_up <- Reduce(intersect, geneList_up)
} else {
  intersectionGenes_up <- c()
}

if (length(geneList_down) > 0) {
  intersectionGenes_down <- Reduce(intersect, geneList_down)
} else {
  intersectionGenes_down <- c()
}

intersectionCount_up <- length(intersectionGenes_up)
intersectionCount_down <- length(intersectionGenes_down)

cat(sprintf("\nCo-upregulated genes: %d\n", intersectionCount_up))
cat(sprintf("Co-downregulated genes: %d\n", intersectionCount_down))

# --- 绘制共同上调基因Venn图 ---
if (length(geneList_up) > 1) {
  pdf(file = file.path(outputDir, "venn_UpRegulated.pdf"), width = 6, height = 6)
  venn_plot_up <- ggvenn(
    geneList_up,
    show_percentage = TRUE,
    stroke_color = "white",
    stroke_size = 0.5,
    fill_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList_up)],
    set_name_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList_up)],
    set_name_size = 2.5,
    text_size = 2
  )
  print(venn_plot_up)
  dev.off()
  cat("Saved: venn_UpRegulated.pdf\n")
}

# --- 绘制共同下调基因Venn图 ---
if (length(geneList_down) > 1) {
  pdf(file = file.path(outputDir, "venn_DownRegulated.pdf"), width = 6, height = 6)
  venn_plot_down <- ggvenn(
    geneList_down,
    show_percentage = TRUE,
    stroke_color = "white",
    stroke_size = 0.5,
    fill_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList_down)],
    set_name_color = c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(geneList_down)],
    set_name_size = 2.5,
    text_size = 2
  )
  print(venn_plot_down)
  dev.off()
  cat("Saved: venn_DownRegulated.pdf\n")
}

# --- 输出共同上调和下调基因的统一CSV ---
if (length(intersectionGenes_up) > 0 || length(intersectionGenes_down) > 0) {
  # 创建上调基因数据框
  if (length(intersectionGenes_up) > 0) {
    intersection_data_up <- data.frame(Gene = intersectionGenes_up,
                                       Regulation = "Upregulated",
                                       stringsAsFactors = FALSE)

    # 从每个样本的数据中提取交集基因的logFC值
    for (sampleName in names(geneData_up)) {
      sampleData <- geneData_up[[sampleName]]
      logfc_values <- c()
      for (gene in intersectionGenes_up) {
        idx <- which(sampleData$Gene == gene)
        if (length(idx) > 0) {
          logfc_values <- c(logfc_values, sampleData$logFC[idx[1]])
        } else {
          logfc_values <- c(logfc_values, NA)
        }
      }
      intersection_data_up[[paste0(sampleName, "_logFC")]] <- logfc_values
    }
  } else {
    intersection_data_up <- data.frame()
  }

  # 创建下调基因数据框
  if (length(intersectionGenes_down) > 0) {
    intersection_data_down <- data.frame(Gene = intersectionGenes_down,
                                         Regulation = "Downregulated",
                                         stringsAsFactors = FALSE)

    # 从每个样本的数据中提取交集基因的logFC值
    for (sampleName in names(geneData_down)) {
      sampleData <- geneData_down[[sampleName]]
      logfc_values <- c()
      for (gene in intersectionGenes_down) {
        idx <- which(sampleData$Gene == gene)
        if (length(idx) > 0) {
          logfc_values <- c(logfc_values, sampleData$logFC[idx[1]])
        } else {
          logfc_values <- c(logfc_values, NA)
        }
      }
      intersection_data_down[[paste0(sampleName, "_logFC")]] <- logfc_values
    }
  } else {
    intersection_data_down <- data.frame()
  }

  # 合并上调和下调数据
  if (nrow(intersection_data_up) > 0 && nrow(intersection_data_down) > 0) {
    intersection_data_combined <- rbind(intersection_data_up, intersection_data_down)
  } else if (nrow(intersection_data_up) > 0) {
    intersection_data_combined <- intersection_data_up
  } else {
    intersection_data_combined <- intersection_data_down
  }

  # 输出合并的CSV文件
  if (nrow(intersection_data_combined) > 0) {
    write.csv(intersection_data_combined,
              file = file.path(outputDir, "IntersectionGenes_Combined.csv"),
              row.names = FALSE, quote = FALSE)
    cat("Saved: IntersectionGenes_Combined.csv\n")
  }
}

# --- 输出共同上调基因TXT ---
if (length(intersectionGenes_up) > 0) {
  write.table(intersectionGenes_up, file = file.path(outputDir, "IntersectionGenes_UpRegulated.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  cat("Saved: IntersectionGenes_UpRegulated.txt\n")
}

# --- 输出共同下调基因TXT ---
if (length(intersectionGenes_down) > 0) {
  write.table(intersectionGenes_down, file = file.path(outputDir, "IntersectionGenes_DownRegulated.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  cat("Saved: IntersectionGenes_DownRegulated.txt\n")
}

# --- 绘制共同上调基因柱状图 ---
if (length(geneList_up) > 1) {
  setSizes_up <- sapply(geneList_up, length)
  pdf(file = file.path(outputDir, "barplot_UpRegulated.pdf"), width=6, height=4)
  print(ggplot(data.frame(Set=names(setSizes_up), Size=setSizes_up), aes(x=Set, y=Size, fill=Set)) +
    geom_bar(stat="identity") +
    theme_minimal() + ylab("Gene Count") + xlab("") + ggtitle("Co-upregulated Genes") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          axis.title = element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8)))
  dev.off()
  cat("Saved: barplot_UpRegulated.pdf\n")
}

# --- 绘制共同下调基因柱状图 ---
if (length(geneList_down) > 1) {
  setSizes_down <- sapply(geneList_down, length)
  pdf(file = file.path(outputDir, "barplot_DownRegulated.pdf"), width=6, height=4)
  print(ggplot(data.frame(Set=names(setSizes_down), Size=setSizes_down), aes(x=Set, y=Size, fill=Set)) +
    geom_bar(stat="identity") +
    theme_minimal() + ylab("Gene Count") + xlab("") + ggtitle("Co-downregulated Genes") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          axis.title = element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8)))
  dev.off()
  cat("Saved: barplot_DownRegulated.pdf\n")
}

cat("\nAnalysis complete!\n")
cat(sprintf("Co-upregulated genes: %d\n", intersectionCount_up))
cat(sprintf("Co-downregulated genes: %d\n", intersectionCount_down))

