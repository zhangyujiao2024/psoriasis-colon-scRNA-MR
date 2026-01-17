# ========================================================================
# 森林图绘制 - 为每个CSV文件单独绘制森林图
# ========================================================================
# 引用包
library(grid)
library(readr)
library(forestploter)

# 设置工作目录
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\21.森林图"
setwd(workDir)

cat("工作目录已设置为:", getwd(), "\n")

# ========================================================================
# 自动读取目录下的CSV文件
# ========================================================================
files <- dir()                           # 获取目录下所有文件
files <- grep("csv$", files, value = TRUE)  # 提取csv结尾的文件

cat(sprintf("找到 %d 个CSV文件:\n", length(files)))
print(files)

# ========================================================================
# 循环处理每个CSV文件，分别绘制森林图
# ========================================================================
for(csv_file in files){
  cat("\n========================================\n")
  cat("正在处理:", csv_file, "\n")

  # 读取当前CSV文件
  data <- read.csv(csv_file, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

  cat(sprintf("读取了 %d 行数据\n", nrow(data)))

  # 检查数据是否为空
  if(nrow(data) == 0) {
    warning(paste("文件", csv_file, "没有数据，跳过"))
    next
  }

  # 显示数据的列名
  cat("数据列名:", paste(colnames(data), collapse = ", "), "\n")

  # 检查是否存在Inverse variance weighted方法
  if(!"Inverse variance weighted" %in% data$method) {
    cat("可用的方法包括：", unique(data$method), "\n")
    warning(paste("文件", csv_file, "中没有找到'Inverse variance weighted'方法，跳过"))
    next
  }

  # 只保留Inverse variance weighted方法
  data <- data[data$method == "Inverse variance weighted", ]

  cat(sprintf("筛选后保留 %d 行IVW方法的数据\n", nrow(data)))

  # 检查筛选后是否还有数据
  if(nrow(data) == 0) {
    warning(paste("文件", csv_file, "筛选后没有数据，跳过"))
    next
  }

  # ========================================================================
  # 数据整理和格式化
  # ========================================================================
  # 重新计算lineVec
  lineVec <- cumsum(c(1, table(data[, c('exposure','outcome')])))

  # 对数据进行整理
  data$' ' <- paste(rep(" ", 10), collapse = " ")
  data$'OR(95% CI)' <- ifelse(is.na(data$or), "", sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))
  data$pval <- ifelse(data$pval < 0.001, "<0.001", sprintf("%.3f", data$pval))
  data$exposure <- ifelse(is.na(data$exposure), "", data$exposure)
  data$outcome <- ifelse(is.na(data$outcome), "", data$outcome)
  data$nsnp <- ifelse(is.na(data$nsnp), "", data$nsnp)

  # 将方法名称从"Inverse variance weighted"改为"IVW"
  data$method <- "IVW"

  # 处理重复数据的安全方法
  data2 <- data[, c('exposure','outcome')]
  duplicate_rows <- duplicated(data2)
  if(sum(duplicate_rows) > 0) {
    data[duplicate_rows, "exposure"] <- ""
    data[duplicate_rows, "outcome"] <- ""
  }

  # ========================================================================
  # 准备图形参数
  # ========================================================================
  tm <- forest_theme(
    base_size = 15,   # 图形整体的大小
    # 可信区间的形状、线条类型、宽度、颜色、两端竖线高度
    ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2,
    # 参考线条的形状、宽度、颜色
    refline_lty = "dashed", refline_lwd = 1, refline_col = "grey20",
    # x轴刻度字体的大小
    xaxis_cex = 0.8,
    # 脚注大小、颜色
    footnote_cex = 0.6, footnote_col = "blue"
  )

  # ========================================================================
  # 绘制森林图
  # ========================================================================
  cat("正在绘制森林图...\n")

  plot <- forestploter::forest(
    data[, c("exposure","outcome","nsnp","method","pval"," ","OR(95% CI)")],
    est = data$or,
    lower = data$or_lci95,
    upper = data$or_uci95,
    ci_column = 6,     # 可信区间所在的列
    ref_line = 1,      # 参考线条的位置
    xlim = c(0, 3),    # X轴的范围
    theme = tm         # 图形的参数
  )

  # ========================================================================
  # 修改图形样式
  # ========================================================================
  # 修改图形中可信区间的颜色（现在只有一种方法，可以使用固定颜色）
  boxcolor <- "#E64B35"  # 使用固定颜色，因为只有一种方法
  for(i in 1:nrow(data)){
    plot <- edit_plot(plot, col = 6, row = i, which = "ci", gp = gpar(fill = boxcolor, fontsize = 25))
  }

  # 设置pvalue的字体（显著性p<0.05的加粗）
  pos_bold_pval <- which(as.numeric(gsub('<',"", data$pval)) < 0.05)
  if(length(pos_bold_pval) > 0){
    for(i in pos_bold_pval){
      plot <- edit_plot(plot, col = 5, row = i, which = "text", gp = gpar(fontface = "bold"))
    }
  }

  # 在图形中增加线段
  plot <- add_border(plot, part = "header", row = 1, where = "top", gp = gpar(lwd = 2))
  plot <- add_border(plot, part = "header", row = lineVec, gp = gpar(lwd = 1))

  # 设置字体大小, 并且将文字居中
  plot <- edit_plot(plot, col = 1:ncol(data), row = 1:nrow(data), which = "text", gp = gpar(fontsize = 12))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"), part = "header",
                    x = unit(0.5, "npc"))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"),
                    x = unit(0.5, "npc"))

  # ========================================================================
  # 输出图形 - 使用CSV文件名命名
  # ========================================================================
  # 从CSV文件名中提取基础名称（去掉.csv扩展名）
  base_name <- tools::file_path_sans_ext(csv_file)
  output_file <- paste0(base_name, "_forest_plot.pdf")

  pdf(output_file, width = 12, height = 20)
  print(plot)
  dev.off()

  cat("森林图已保存为:", output_file, "\n")
}

cat("\n========================================\n")
cat("所有森林图绘制完成!\n")
cat(sprintf("共处理了 %d 个CSV文件\n", length(files)))
