# 设置工作目录
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\31.提取细胞亚型"  # 工作目录
setwd(workDir)  # 切换到指定工作目录



# 读取 Seurat 对象（文件在当前目录）
seurat_obj <- readRDS("filtered_single_cell_data.rds")

# 读取细胞ID映射文件
cell_id_mapping <- read.csv("Cell_ID_to_CellType_Mapping.csv", header = TRUE)

# 提取细胞名称（第一列是细胞ID）
cell_ids_to_extract <- cell_id_mapping[, 1]

# 确保这些细胞ID在 Seurat 对象中
cells_in_seurat <- colnames(seurat_obj)
cells_to_extract <- intersect(cell_ids_to_extract, cells_in_seurat)

# 提取指定的细胞（使用 seurat_obj[, cells_to_extract] 语法）
seurat_obj_extracted <- seurat_obj[, cells_to_extract]

# 保存提取后的 Seurat 对象
saveRDS(seurat_obj_extracted, file = paste0(workDir, "/Extracted_Cells.rds"))

# 提示完成
message("已成功提取细胞并保存为 RDS 文件！")
