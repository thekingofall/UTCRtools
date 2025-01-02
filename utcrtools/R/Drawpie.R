library(RColorBrewer)

# 获取命令行参数
args = commandArgs(trailingOnly=TRUE)

# 检查参数数量
if (length(args) != 2) {
    print("Usage:<filename><Give the png a name>")
} else {
    # 读取输入参数
    filename = args[1]
    pngname = args[2]
    
    # 读取数据文件
    TCR_condition2 <- read.table(filename, sep = "\t", header = TRUE)
    
    # 分离TRA和TRB数据
    TRA = TCR_condition2[which(substring(TCR_condition2$V.segments, 1, 3) == "TRA"), ]
    TRB = TCR_condition2[which(substring(TCR_condition2$V.segments, 1, 3) == "TRB"), ]
    
    # 定义绘制饼图函数
    Runpie <- function(TCR, name) {
        # 计算克隆计数
        count = as.data.frame(table(TCR$Count))
        count$Var1 <- as.numeric(as.character(count$Var1))
        
        # 计算总和和百分比
        count$sum <- count$Var1 * count$Freq
        count <- count[order(count$sum, decreasing = TRUE), ]
        count$percent <- signif((count$sum / sum(count$sum)) * 100, 2)
        count$TCRname <- paste0("Clone_", count$Var1)
        count$label <- paste(count$TCRname, "(", count$percent, "%)")
        
        # 保存结果
        write.csv(count, paste0(pngname, "_pure_", name, ".csv"))
        
        # 绘制饼图
        png(paste0(pngname, "_", name, ".png"), width = 600 * 4, height = 3 * 600, res = 72 * 3)
        pie(count$sum, 
            labels = count$Var1[1:5],
            col = colorRampPalette(rev(brewer.pal(10, 'Spectral')))(34),
            border = "white",
            angle = 10,
            main = paste0(name, "\n", "  CDR3 total ", sum(count$sum), "\n", "  CDR3 types ", nrow(TCR)))
        legend("left", 
               legend = count$label[1:5], 
               cex = 1.0,
               fill = colorRampPalette(rev(brewer.pal(10, 'Spectral')))(34)[1:4],
               border = "white",
               bty = "n")
        dev.off()
    }
    
    # 根据数据类型选择绘制方式
    if (all(substring(TCR_condition2$V.segments, 1, 3) == "TRB")) {
        Runpie(TRB, name = "TCRβ")
    } else if (all(substring(TCR_condition2$V.segments, 1, 3) == "TRA")) {
        Runpie(TRA, name = "TCRα")
    } else {
        Runpie(TRA, name = "TCRα")
        Runpie(TRB, name = "TCRβ")
    }
}