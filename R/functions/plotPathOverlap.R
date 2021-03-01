plotPathOverlap <- function(path1, path2, filePath = "", ratio = 0.8) {
  path1 <- as.data.frame(path1)
  path2 <- as.data.frame(path2)
  path_length   <- ncol(path1)
  colnames(path1) <- sprintf("t%03d", 0:(path_length-1)) # Replace t..0 with t000
  colnames(path2) <- sprintf("t%03d", 0:(path_length-1)) # Replace t..0 with t000
  path1$symbols  <- rownames(path1)
  path2$symbols  <- rownames(path2)

  # Rearranges data frame 1
  pathGather1         <- gather(path1, "t", "value", 1:path_length)
  pathGather1$symbols <- factor(pathGather1$symbols, levels = rownames(path1)[nrow(path1):1])
  pathGather1$t       <- factor(pathGather1$t)
  pathGather1$value   <- factor(pathGather1$value, levels = c("0", "1"))

  # Rearranges data frame 2
  pathGather2         <- gather(path2, "t", "value", 1:path_length)
  pathGather2$symbols <- factor(pathGather2$symbols, levels = rownames(path2)[nrow(path2):1])
  pathGather2$t       <- factor(pathGather2$t)
  pathGather2$value   <- factor(pathGather2$value, levels = c("0", "1", "2"))
  pathGather2$value[pathGather2$value == "1"] <- "2"

  # Plot first simulation
  compColors <- c("0" = "transparent", "1" = "red", "2" = "steelblue")
  p          <- ggplot() +
                  geom_tile(data = pathGather1, mapping = aes(t, symbols, fill = value), colour = "white", alpha = 0.7) +
                  geom_tile(data = pathGather2, mapping = aes(t, symbols, fill = value), colour = "white", alpha = 0.5) +
                  scale_fill_manual(values = compColors, labels=c("Neither", "1", "2", "1+2", "NA"), drop=FALSE)
  base_size <- 8
  ratio <- 0.8
  p         <- p + theme_grey(base_size = base_size) + theme(panel.background = element_blank()) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
                coord_fixed(ratio=ratio) + guides(fill = guide_legend(title="Active in path:"))

  # Save plot to file
  if (filePath != "") {
    filePath <- gsub(".pdf$", "", filePath, ignore.case = TRUE) # Remove extension if present
    ggsave(paste0(filePath, ".pdf"), plot = last_plot(), height=(5 + length(levels(pathGather1$symbols)) * 0.25), width=(10 + length(levels(pathGather1$t)) * 0.45), scale = 1, units = "cm")
  }

  return(p)
}