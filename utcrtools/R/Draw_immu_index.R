# Remove all objects from environment
#rm(list = ls())

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

# Check if correct number of arguments provided
if (length(args) != 2) {
  print("Usage:<filename><Give a name for the folder of result >")
} else {
  # Load required libraries
  library(immunarch)
  library(gridExtra)
  library(ggsci)
  library(RColorBrewer)
  library(ggplot2)
  
  # Get input arguments
  filename = args[1]
  pngname = args[2]
  
  # Load data and create output directory
  immdata <- repLoad(filename)
  filetest <- "Results_for_summary"
  file1 <- paste0(pngname, "_", filetest)
  dir.create(file1)
  
  # Process sample names
  names(immdata$data) <- substring(names(immdata$data), 1, 16)
  
  # Find common CDR3 sequences across samples
  data_same = list()
  for (i in 1:length(immdata$data)) {
    data_same[[i]] <- immdata$data[[i]]$CDR3.aa
  }
  data_same_seq = Reduce(intersect, data_same)
  write.csv(data_same_seq, paste0(file1, "/", "All_Same_seq_in.csv"))
  
  # Create CDR3 length distribution plots
  data_length_plot <- list()
  for (i in 1:length(immdata$data)) {
    cdr_length = data.frame(length = nchar(immdata$data[[i]]$CDR3.aa))
    cdr_length_table = as.data.frame(table(cdr_length))
    cdr_length_table$cdr_length <- as.factor(cdr_length_table$cdr_length)
    
    dlp <- ggplot(cdr_length_table) +
      geom_bar(aes(x = cdr_length, y = Freq, fill = cdr_length), stat = "identity") +
      theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25),
            axis.line = element_line(colour = "black", size = 0.25),
            axis.title = element_text(size = 13, face = "plain", color = "black"),
            axis.text = element_text(size = 12, face = "plain", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none") +
      scale_fill_viridis_d() +
      ylab("CDR3 Freq") +
      ggtitle(names(immdata$data)[i])
    
    data_length_plot[[i]] <- dlp
  }
  
  gd <- do.call(grid.arrange, data_length_plot)
  ggsave(plot = gd, filename = paste0(file1, "/", "All_cdr3_length.pdf"), width = 16, height = 16)
  
  # Calculate and visualize sample overlaps
  repOverdata <- as.data.frame(repOverlap(immdata$data))
  data_overlap = as.data.frame(repOverlap(immdata$data))
  data_overlap[upper.tri(data_overlap)] <- NA
  
  corda <- data_overlap
  corda$y <- rownames(corda)
  da <- melt(data = corda) %>% na.omit()
  da$variable <- factor(da$variable, levels = unique(da$variable))
  da$y <- factor(da$y, levels = unique(da$y))
  
  mycolor2 = colorRampPalette(c("#1E3163", "#00C1D4", "#FFED99", "#FF7600"))(10)
  colnames(da)[3] <- "Overlap_values"
  
  pw <- ggplot(da, aes(x = variable, y = y, fill = Overlap_values)) +
    geom_tile(aes(x = variable, y = y, fill = Overlap_values),
              color = 'white', size = 0.6, alpha = 0.8) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    scale_x_discrete(position = 'top') +
    scale_fill_gradientn(colors = mycolor2) +
    geom_text(aes(label = Overlap_values)) +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("") + ylab("")
  
  all_overlap_plot <- pw + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(pw, filename = paste0(file1, "/", "All_overlaplot.pdf"), width = 10, height = 10)
  write.csv(data_overlap, paste0(file1, "/", "All_overlap_value.csv"))
  
  # Analyze and visualize gene usage
  data_geneusage <- list()
  data_geneusage_csv <- data.frame()
  
  for (i in 1:length(immdata$data)) {
    print(i)
    genedata <- as.data.frame(geneUsage(immdata$data[[i]]))
    genedata$group <- rep(names(immdata$data)[i], nrow(genedata))
    
    pn <- ggplot(genedata, aes(x = Names, y = Clones, fill = Names)) +
      geom_bar(stat = "identity") +
      theme(legend.position = "none") +
      xlab("") +
      theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25),
            axis.title = element_text(size = 13, face = "plain", color = "black"),
            axis.text = element_text(size = 12, face = "plain", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none") +
      scale_fill_viridis_d() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
      ggtitle(names(immdata$data)[i])
    
    data_geneusage[[i]] <- pn
    data_geneusage_csv <- rbind(data_geneusage_csv, genedata)
  }
  
  write.csv(data_geneusage_csv, paste0(file1, "/", "ALl_geneusage.csv"))
  write.csv(repOverdata, paste0(file1, "/", "All_overlap_value.csv"))
  
  plot_geneuase = do.call(grid.arrange, data_geneusage)
  ggsave(plot = plot_geneuase, paste0(file1, "/", "plot1_Summary_plot_geneuase.pdf"), height = 40, width = 40)
  
  # Create exploration plots
  myfile2 <- paste0(file1, "/1.Explore")
  dir.create(myfile2)
  
  re1 = repExplore(immdata$data, "len", .col = "aa") %>% vis() + scale_fill_viridis_d()
  ggsave(plot = re1, paste0(myfile2, "/", "1_len.pdf"))
  
  re2 = repExplore(immdata$data, "volume") %>% vis() + scale_fill_viridis_d()
  ggsave(plot = re2, paste0(myfile2, "/", "2_volume.pdf"))
  
  re3 = repExplore(immdata$data, .method = "count", .col = "aa") %>% vis() + scale_fill_viridis_d() + scale_colour_viridis_d()
  ggsave(plot = re3, paste0(myfile2, "/", "3_volume.pdf"))
  
  re4 = repExplore(immdata$data, "clones") %>% vis() + scale_fill_viridis_d()
  ggsave(plot = re4, paste0(myfile2, "/", "4_clones.pdf"))
  
  re_all = grid.arrange(re1, re2, re3, re4)
  ggsave(plot = re_all, filename = paste0(file1, "/All_explot.pdf"), width = 14, height = 14)
  
  # Calculate and visualize diversity metrics
  method_repversity <- c("div", "hill", "gini.simp", "raref", "d50", "chao1")
  myfile = paste0(file1, "/", "2.results_of_repversity")
  dir.create(myfile)
  data_repversity <- list()
  
  for (i in 1:length(method_repversity)) {
    tem_rep <- method_repversity[i]
    print(tem_rep)
    p = repDiversity(immdata$data, tem_rep) %>% vis() + scale_fill_viridis_d() + scale_colour_viridis_d()
    data_repversity$data[[i]] <- as.data.frame(repDiversity(immdata$data, tem_rep))
    ggsave(plot = p, filename = paste0(myfile, "/", tem_rep, ".pdf"))
    assign(paste0(tem_rep, "_plot"), p)
  }
  
  data_rep_plot <- list(chao1_plot, gini.simp_plot, div_plot, hill_plot, raref_plot, d50_plot)
  dall <- do.call(grid.arrange, data_rep_plot)
  ggsave(plot = dall, paste0(file1, "/", "All_diversity_index.pdf"), height = 20, width = 20)
  
  # Process diversity data
  names(data_repversity$data) <- c("div", "hill", "gini.simp", "raref", "d50", "chao1")
  
  data_hill <- data_repversity$data[["hill"]]
  data_hill <- dcast(data_hill, Sample ~ Q)
  lmn <- length(colnames(data_hill)) - 1
  colnames(data_hill) <- c("Sample", paste0("hill_Q", 1:lmn))
  data_repversity$data[["hill"]] <- data_hill
  
  colnames(data_repversity$data[["div"]]) <- c("Sample", "div_value")
  colnames(data_repversity$data[["gini.simp"]]) <- c("Sample", "gini.simp_value")
  
  data_repversity$data[["d50"]][, 3] <- rownames(data_repversity$data[["d50"]])
  colnames(data_repversity$data[["d50"]]) <- c("d50_Clones", "d50_Percentage", "Sample")
  
  data_repversity$data[["chao1"]][, 5] <- rownames(data_repversity$data[["chao1"]])
  colnames(data_repversity$data[["chao1"]]) <- c(paste0("chao1_", colnames(data_repversity$data[["chao1"]])[1:4]), "Sample")
  
  # Merge diversity metrics
  names = c("hill", "gini.simp", "d50", "chao1")
  data_tem <- as.data.frame(data_repversity$data[["div"]])
  for (i in names) {
    data_tem <- merge(data_tem, data_repversity$data[[i]], by = "Sample")
  }
  
  # Define plotting function
  mycolor = colorRampPalette(rev(brewer.pal(10, "Paired")))(34)
  plotdef <- function(data, xnames, yvalue, linelong = c(1), ylabname = "Values") {
    p1 = ggplot(data = data, aes(x = get(xnames), y = get(yvalue), fill = get(xnames))) +
      geom_bar(stat = "identity") +
      theme(legend.position = "none") +
      xlab("") +
      ylab(ylabname) +
      theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25),
            axis.title = element_text(size = 20, face = "plain", color = "black"),
            axis.text = element_text(size = 12, face = "plain", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none") +
      scale_fill_manual(values = mycolor) +
      theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = .5)) +
      geom_hline(yintercept = linelong, linetype = "dashed", alpha = 0.3)
    
    return(p1)
  }
  
  # Calculate and visualize Gini coefficient
  data_gini <- data.frame(names = rownames(repDiversity(immdata$data)), ginivalues = repDiversity(immdata$data, "gini"))
  pgini <- plotdef(data_gini, xnames = "names", yvalue = "ginivalues", ylabname = "Gini coefficient")
  ggsave(plot = pgini, filename = paste0(myfile, "/", "gini.plot.pdf"))
  write.csv(data_gini, paste0(myfile, "/", "gini.csv"))
  
  # Save results
  write.csv(data_tem, paste0(myfile, "/", "half_repversity_tem.csv"))
  write.csv(data_repversity$data[["raref"]], paste0(myfile, "/raref.csv"))
  
  print("thanks for use--------MLP20211115")
}