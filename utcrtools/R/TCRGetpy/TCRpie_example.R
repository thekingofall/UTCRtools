analyze_TCR <- function(data_path, output_dir = "TCR_analysis", output_prefix = "TCR") {
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 完整输出路径
  output_path <- file.path(output_dir, output_prefix)
  
  # 读取数据
  TRall <- fread(data_path)
  TRB <- TRall[,1:2]
  colnames(TRB) <- c("cloneId","Count")
  
  # 1. 密度分布图
  p1 <- data.frame(count = as.numeric(names(table(TRB$Count))), 
                   freq = as.numeric(table(TRB$Count))) %>%
    ggplot(aes(x = log10(count), weight = freq)) +
    geom_density(fill = "#a6e3e9", alpha = 0.6, color = "#364f6b",size=1) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Log10 Density Distribution of Count",
      x = "log10(Count)",
      y = "Density"
    ) +
    # geom_vline(aes(xintercept = weighted.mean(log10(count), freq)),
    #            linetype = "dashed", color = "#E74C3C", alpha = 0.7,size=1) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0))
  
  ggsave(paste0(output_path, "_density.png"), 
         plot = p1, 
         width = 6, 
         height = 4, 
         dpi = 300)
  
  # 2. 累积分布图
  df <- data.frame(count = as.numeric(names(table(TRB$Count))), 
                   freq = as.numeric(table(TRB$Count)))
  
  p2 <- ggplot(df, aes(x = log10(count), weight = freq)) +
    stat_ecdf(geom = "line", size = 1.2, color = "#2E86C1") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    ) +
    labs(
      title = "Cumulative Distribution of Count",
      x = "log10(Count)",
      y = "Cumulative Probability"
    ) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2)) +
    geom_hline(yintercept = c(0.2, 0.6, 0.8), 
               linetype = "dashed", 
               color = "gray70",
               alpha = 0.5)
  
  ggsave(paste0(output_path, "_cumulative.png"), 
         plot = p2, 
         width = 6, 
         height = 4, 
         dpi = 300)
  
  # 3. 饼图函数
  make_pie <- function(TCR, name) {
    count = as.data.frame(table(TCR$Count))
    count$Var1 <- as.numeric(as.character(count$Var1))
    count$sum <- count$Var1 * count$Freq
    count <- count[order(count$sum, decreasing = T),]
    count$percent <- signif((count$sum/sum(count$sum))*100, 2)
    
    count$TCRname <- paste0("clone", count$Var1)
    count$TCRname[1:5] <- paste0("Top", 1:5, ": clone", count$Var1[1:5])
    
    count$group <- count$TCRname
    count$group[6:nrow(count)] <- case_when(
      count$Var1[6:nrow(count)] < 50 ~ "Others (<50)",
      count$Var1[6:nrow(count)] >= 50 & count$Var1[6:nrow(count)] < 1000 ~ "Others (50-1000)",
      count$Var1[6:nrow(count)] >= 1000 & count$Var1[6:nrow(count)] < 10000 ~ "Others (1000-10000)",
      count$Var1[6:nrow(count)] >= 10000 & count$Var1[6:nrow(count)] < 50000 ~ "Others (10000-50000)",
      count$Var1[6:nrow(count)] >= 50000 ~ "Others (>50000)"
    )
    
    top5 <- count[1:5,]
    others <- count[6:nrow(count),]
    others_grouped <- aggregate(sum ~ group, data = others, sum)
    others_grouped$percent <- signif((others_grouped$sum/sum(count$sum))*100, 2)
    
    final_data <- rbind(
      data.frame(group = top5$TCRname, sum = top5$sum, percent = top5$percent),
      others_grouped
    )
    
    write.csv(final_data, paste0(output_path, "_", name, "_pure.csv"))
    
    labels <- paste0(final_data$group, " (", final_data$percent, "%)")
    
    png(paste0(output_path, "_", name, "_pie.png"), 
        width=600*4, height=3*600, res=72*3)
    pie(final_data$sum, 
        labels = labels,
        col = colorRampPalette(rev(brewer.pal(10,'Spectral')))(nrow(final_data)),
        border = "white",
        angle = 10,
        main = paste0(name, "\n", "CDR3 total ", sum(count$sum), 
                      "\n", "CDR3 types ", nrow(TCR)))
    dev.off()
  }
  
  make_pie(TRB, name=output_prefix)
  
  plot_VJ_heatmap <- function(data, output_file = "VJ_heatmap.png", width = 12, height = 10) {
    library(pheatmap)
    library(dplyr)
    library(stringr)
    
    # 提取V基因名和J基因名
    data$V <- gsub("\\*.*$", "", data$allVHitsWithScore)
    data$J <- gsub("\\*.*$", "", data$allJHitsWithScore)
    
    # 创建频率矩阵并计算百分比
    freq_matrix <- table(data$V, data$J) %>% 
      as.data.frame.matrix() %>%
      as.matrix()
    
    # 计算总体百分比
    freq_matrix_pct <- freq_matrix/sum(freq_matrix) * 100
    
    # 创建自定义排序函数
    custom_sort <- function(x) {
      # 将名称拆分为基因和数字部分
      gene_nums <- str_match(x, "(TRBV\\d+)-?(\\d*)")
      # 转换为数值进行排序
      gene_base <- as.numeric(str_extract(gene_nums[,2], "\\d+"))
      gene_suffix <- ifelse(gene_nums[,3] == "", 0, as.numeric(gene_nums[,3]))
      # 按基因基础号码和后缀号码排序
      x[order(gene_base, gene_suffix)]
    }
    
    # 对行名进行排序
    row_order <- custom_sort(rownames(freq_matrix_pct))
    freq_matrix_pct <- freq_matrix_pct[row_order,]
    
    # 绘制热图
    pheatmap(freq_matrix_pct,
             color = colorRampPalette(c("navy", "white", "red"))(100),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             angle_col = 45,
             main = "TRV-TRJ Usage (%)",
             filename = output_file,
             width = width,
             height = height)
    
    # 返回频率矩阵（可选）
    return(freq_matrix_pct)
  }
  # 4. VJ热图
  plot_VJ_heatmap(TRall, 
                  output_file = paste0(output_path, "_VJ_heatmap.png"),
                  width = 10, 
                  height = 10)
  
  # 返回输出目录路径
  return(output_dir)
}

## 单文件
analyze_TCR("data/TCR_seq/20240824.TRA.txt", 
            output_dir = "results/20240824.TRA",
            output_prefix = "20240824.TRA")

## 多文件
# 获取所有需要分析的文件
files <- list.files("data/TCR_seq", pattern = "\\.txt$", full.names = TRUE)

# 为每个文件创建对应的输出目录
for(file in files) {
  # 从文件名提取日期和类型
  filename <- basename(file)
  date_type <- tools::file_path_sans_ext(filename)  # 去除.txt后缀
  
  # 创建对应的输出目录
  output_dir <- file.path("TCR_analysis", date_type)
  
  # 分析数据
  message("Processing file: ", filename)
  analyze_TCR(file, 
              output_dir = output_dir,
              output_prefix = date_type)
}