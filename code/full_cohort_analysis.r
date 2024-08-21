# Load required libraries
library(TraMineR)

# create and save down a distribution plot
create_dist_plot <- function(seqact, labels, palette, plot_file = NA) {
  if (!dir.exists("./results")) {
    dir.create("./results")
  }
  # if no file directory provided then save in automatic location
  if (is.na(plot_file)) {
    plot_file <- "./results/State-Distribution-Plot.png"
  }
  
  png(file = plot_file,
      width = 9000, height = 6000, units = "px", res = 900)
  par(mar = c(9, 6, 2, 2), mgp = c(4, 1, 0))
  ylab_text <- sprintf("Relative Frequency (n=%s)", nrow(seqact))
  seqdplot(seqact, with.legend = FALSE,
           main = "State Distribution Plot", border = NA, las = 2, cex.axis = .8,
           xtstep = round(ncol(seqact) / 43, 0), tick.last = TRUE,
           xlab = "Days in proximity to Index Date", ylab = ylab_text)
  
  # Add a vertical line at day 366
  usr <- par("usr")
  segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4],
           col = "black", lty = 2, lwd = 2)
  
  # Add text label at the vertical line
  text(x = 355, y = 0.1, srt = 90, "Index Date")
  
  # Add a legend
  par(xpd = TRUE)
  if(length(alphabet) > 4) {
    ncol_leg <- 3
  } else {
    ncol_leg <- 2
  }
  legend("bottom", inset = c(-0.37, -0.37), ncol = ncol_leg,
         legend = labels,
         fill = palette, bty = "n", x.intersp = 0.2, cex = 0.8)
  dev.off()
}


# create entropy plot
create_entropy_plot <- function(seqact, seqact2 = NULL, plot_file = NA) {
  # if no file directory provided then save in automatic location
  if (is.na(plot_file)) {
    if (!dir.exists("./results")) {
      dir.create("./results")
    }
    plot_file <- "./results/Entropy-Plot.png"
  }
  
  # set up image file details
  png(file = plot_file, width = 9000, height = 6000, units = "px", res = 900)
  
  # create label for y axis
  ylab_text <- sprintf("Entropy of patient trajectories (n=%s)", nrow(seqact))
  
  # create main title based on whether there are two entropy lines or not
  main_title <- if (is.null(seqact2)) {
    "Entropy Plot of Treatment Patterns over time (GLP-1 RA medications grouped)"
  } else {
    "Entropy Plot of Treatment Patterns over time"
  }
  
  # Calculate entropy for the first sequence
  entropy1 <- seqstatd(seqact, weighted = FALSE, with.missing = FALSE, norm = TRUE)$Entropy
  
  # Calculate entropy for the second sequence if provided
  if (!is.null(seqact2)) {
    entropy2 <- seqstatd(seqact2, weighted = FALSE, with.missing = FALSE, norm = TRUE)$Entropy
  }
  
  par(mar = c(7, 6, 2, 2), mgp = c(4, 1, 0), oma = c(0, 0, 1, 0), bty = "n")
  
  # Plot the entropy for the first sequence
  plot(entropy1, type = "l", lwd = 1.5, col = "black", lty = 1,
       main = main_title,
       xlab = "Days in proximity to Index Date", ylab = ylab_text,
       xaxt = 'n', yaxt = 'n', ylim = c(0, 1))
  
  # Customize x and y axes
  xticks <- seq(1, length(entropy1), by = round(length(entropy1) / 43))
  if (tail(xticks, 1) != length(entropy1)) {
    xticks <- c(xticks, length(entropy1))
  }
  axis(1, at = xticks, labels = colnames(seqact)[xticks], las = 2, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
  
  # Add grid lines
  grid(nx = NULL, ny = NULL, col = rgb(0.5, 0.5, 0.5, 0.5), lty = "dotted")
  
  # Add a vertical line at day 366
  usr <- par("usr")
  segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4], col = "red", lty = 2, lwd = 1.5)
  
  # Add text label at the vertical line
  text(x = 355, y = 0.2, srt = 90, "Index Date", col = "red", lwd = 1.5)
  
  # Plot the entropy for the second sequence if provided
  if (!is.null(seqact2)) {
    legend <- if (length(unique(seqact)) > 4) {
      c("Distinct GLP-1 RAs", "GLP-1 RAs grouped")
    } else {
      c("GLP-1 RAs grouped", "Distinct GLP-1 RAs")
    }
    lines(entropy2, type = "l", lwd = 1.5, col = "blue", lty = 3)
    legend("topright", legend = legend,
           col = c("black", "blue"), lwd = 1.5, lty = c(1, 3), bty = "n")
  }
  
  dev.off()
  print(sprintf("Entropy plot saved here: %s", plot_file))
}

