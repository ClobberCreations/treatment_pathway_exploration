### this file contains functions for analysing the models further ###

# function to save down PatIds and cluster groups
save_clusters_id <- function(model, k = 3) {
  # cut into required number of clusters
  clustered <- cutree(model$clust, k = k)
  # create a data fram of Patient IDs and cluster group names
  cluster_df <- data.frame(PatId = names(clustered), Cluster = as.vector(clustered))
  # preview the data frame
  head(cluster_df)
  # filename
  file_name <- sprintf("./results/PatId_%s-clusters_%s.csv", k, model$name)
  # write to csv
  write.csv(cluster_df, file_name, row.names = FALSE)
  print(sprintf("Patient IDs and their cluster saved here: %s", file_name))
}


# factorise clusters with labels
label_clusters <- function(model, k = 3, label = "Cluster") {
  # cut into clusters
  clustered <- cutree(model$clust, k = k)
  # label
  labelled_k <- factor(clustered, labels = paste(label, 1:k))
  # return labeled factors
  return(labelled_k)
}


# create state distribution plot
create_dist_plot <- function(seqact, labels, palette, model_name = "model", clusters = NULL) {
  if (!dir.exists("./results")) {
    dir.create("./results")
  }
  # if no clusters, plot full cohort
  if (is.null(clusters)) {
    plot_file <- sprintf("./results/State-Distribution-Plot_%s.png", model_name)
    png(file = plot_file,
        width = 8000, height = 6000, units = "px", res = 900)
    par(mar = c(9, 6, 2, 2), mgp = c(4, 1, 0))
    ylab_text <- sprintf("Relative Frequency (n=%s)", nrow(seqact))
    seqdplot(seqact, with.legend = FALSE, main = "State Distribution Plot", 
             border = NA, las = 2, cex.axis = .8,
             xtstep = round(ncol(seqact) / 43, 0), tick.last = TRUE,
             xlab = "Days in proximity to Index Date", ylab = ylab_text)
    
    # Add a vertical line at day 366
    usr <- par("usr")
    segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4],
             col = "black", lty = 2, lwd = 2)
    
    # Add text label at the vertical line
    text(x = 355, y = 0.2, srt = 90, "Index Date")
    
    # Add a legend
    par(xpd = TRUE)
    # work out number of columns depending on number of labels
    ncol_leg <- if (length(labels) > 4) 3 else 2
    legend("bottom", inset = c(-0.37, -0.37), ncol = ncol_leg,
           legend = labels,
           fill = palette, bty = "n", x.intersp = 0.2, cex = 0.8)
  } else {
    # get number of clusters
    num_clusters <- length(unique(clusters))
    plot_file <- sprintf("./results/State-Distribution-Plot-by-%s-clusters_%s.png", num_clusters, model_name)
    png(file = plot_file,
        width = 8000, height = 6000, units = "px", res = 900)
    par(mar=c(6,6,2,2), mgp = c(4, 1, 0))
    if (num_clusters %% 2 == 0) {
      if (num_clusters == 2) {
        png(file = plot_file,
            width = 8000, height = 5000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
        layout_mat <- matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2,
                             byrow = TRUE)
        layout(mat = layout_mat, heights=c(3, 1.5))
      } else if (num_clusters == 4) {
        png(file = plot_file,
            width = 8000, height = 9000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
        layout_mat <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, ncol = 2,
                             byrow = TRUE)
        layout(mat = layout_mat, heights=c(4, 4, 1.5), widths=c(8, 8))
      } else if (num_clusters == 6) {
        png(file = plot_file,
            width = 8000, height = 10000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
        layout(matrix(c(1,2,3,4,5,6,7,7), nrow = 4, ncol=2, byrow=TRUE), heights=c(3, 3, 3, 1.5))
      }
      #par(mai=rep(0.5, 0.5))
    } else {
      par(mfrow = c(ceiling(num_clusters/2), 2), mar = c(6, 6, 2, 2), mgp = c(4, 1, 0))
    }
    
    for (cluster in unique(clusters)) {
      cluster_indices <- which(clusters == cluster)
      seqdplot(seqact[cluster_indices, ], with.legend = FALSE,
               main = paste("State Distribution Plot -", cluster),
               border = NA, las = 2, cex.axis = .6, cex.main = 1.2,
               cex.lab = 0.8,
               xtstep = ceiling(ncol(seqact) / 43), tick.last = TRUE,
               xlab = "Days in proximity to Index Date", ylab = NULL)
      
      # Add a vertical line at day 366
      usr <- par("usr")
      segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4],
               col = "black", lty = 2, lwd = 1)
      
      # Add text label at the vertical line
      text(x = 345, y = 0.3, srt = 90, "Index Date", cex = 0.9)
    }
    
    if (num_clusters %% 2 == 0) {
      sizing <- if (num_clusters == 2) 0.8 else 1
      # Plot the legend in the final plot space
      ncol_leg <- if (length(labels) > 4) 3 else 2
      seqlegend(seqact, position = "center", bty="n", x.intersp = 0.2, cex = sizing,
                ncol = ncol_leg)
    } else {
      plot.new()
      par(xpd = TRUE)
      ncol_leg <- if (length(labels) > 4) 3 else 2
      legend("center", inset = 0,
             legend = labels,
             fill = palette, bty = "n", x.intersp = 0.2, cex = 0.8)
    }
    
  }
  print(sprintf("Distribution plot saved here: %s", plot_file))
  #dev.off()
}


# create plot of most representative sequences
create_rep_plot <- function(seqact, model, labels, palette, clusters = NULL) {
  if (!dir.exists("./results")) {
    dir.create("./results")
  }
  
  if (is.null(clusters)) {
    plot_file <- sprintf("./results/Representative Sequences-Plot-all_%s.png", model$name)
    png(file = plot_file,
        width = 8000, height = 7000, units = "px", res = 900)
    par(mar = c(9, 6, 2, 2), mgp = c(4, 2, 0))
    seqrplot(seqact, diss = model$dissim, with.legend = FALSE, main = "Representative Sequences for Treatment Patterns", 
             border = NA, las = 2, cex.axis = .8,
             xtstep = round(ncol(seqact) / 43, 0), tick.last = TRUE,
             xlab = "Days in proximity to Index Date")
    
    if (ncol(seqact) == 731) {
      # Add a vertical line at day 366, reducing the height by one-third
      usr <- par("usr")
      segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4] * 0.5,
               col = "black", lty = 2, lwd = 2)
      
      # Adjust text label position
      text(x = 355, y = 0.5, srt = 90, "Index Date")
    }
    
    
    # Add a legend
    par(xpd = TRUE)
    ncol_leg <- if (length(labels) > 4) 3 else 2
    legend("bottom", inset = c(-0.37, -0.37), ncol = ncol_leg,
           legend = labels,
           fill = palette, bty = "n", x.intersp = 0.2, cex = 0.8)
  } else {
    num_clusters <- length(unique(clusters))
    plot_file <- sprintf("./results/Representative Sequences-Plot-%s-by-%s-clusters.png", model$name, num_clusters)
    png(file = plot_file,
        width = 8000, height = 6000, units = "px", res = 900)
    par(mar=c(6,6,2,2), mgp = c(4, 1, 0))
    if (num_clusters %% 2 == 0) {
      if (num_clusters == 2) {
        png(file = plot_file,
            width = 8000, height = 5000, units = "px", res = 900)
        par(mar=c(6,3,2,2), mgp = c(4, 0.5, 0))
        layout_mat <- matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2,
                             byrow = TRUE)
        layout(mat = layout_mat, heights=c(3, 1.5))
      } else if (num_clusters == 4) {
        png(file = plot_file,
            width = 8000, height = 9000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
        layout_mat <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, ncol = 2,
                             byrow = TRUE)
        layout(mat = layout_mat, heights=c(4, 4, 1.5), widths=c(8, 8))
      } else if (num_clusters == 6) {
        png(file = plot_file,
            width = 8000, height = 10000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
        layout(matrix(c(1,2,3,4,5,6,7,7), nrow = 4, ncol=2, byrow=TRUE), heights=c(3, 3, 3, 1.5))
      }
      #par(mai=rep(0.5, 0.5))
    } else {
      par(mfrow = c(ceiling(num_clusters/2), 2), mar = c(6, 6, 2, 2), mgp = c(4, 1, 0))
    }
    for (cluster in unique(clusters)) {
      cluster_indices <- which(clusters == cluster)
      seqrplot(seqact[cluster_indices, ], diss = model$dissim, with.legend = FALSE,
               main = paste("Representative Sequences -", cluster),
               border = NA, las = 2, cex.axis = .6, cex.main = 0.9,
               cex.lab = 0.8,
               xtstep = ceiling(ncol(seqact) / 43), tick.last = TRUE,
               xlab = "Days in proximity to Index Date")
      
      if (ncol(seqact) == 731) {
        # Add a vertical line at day 366, reducing the height by one-third
        usr <- par("usr")
        segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4] * 0.5,
                 col = "black", lty = 2, lwd = 1)
        
        # Adjust text label position
        text(x = 355, y = 0.5, srt = 90, "Index Date", cex = 0.7)
      }
      
    }
    
    if (num_clusters %% 2 == 0) {
      sizing <- if (num_clusters == 2) 0.8 else 1
      # Plot the legend in the final plot space
      ncol_leg <- if (length(labels) > 4) 3 else 2
      seqlegend(seqact, position = "center", bty="n", x.intersp = 0.2, cex = sizing,
                ncol = ncol_leg)
    } else {
      plot.new()
      par(xpd = TRUE)
      ncol_leg <- if (length(labels) > 4) 3 else 2
      legend("center", inset = 0,
             legend = labels,
             fill = palette, bty = "n", x.intersp = 0.2, cex = 0.8)
    }
  }
  print(sprintf("Representative sequences plot saved here: %s", plot_file))
  #dev.off()
}



# create plot of mean distribution of states
create_mean_time_plot <- function(seqact, model, labels, palette, clusters = NULL) {
  if (!dir.exists("./results")) {
    dir.create("./results")
  }
  
  if (is.null(clusters)) {
    plot_file <- sprintf("./results/mean-time-state-plot-all_%s.png", model$name)
    png(file = plot_file,
        width = 8000, height = 7000, units = "px", res = 900)
    par(mar = c(9, 6, 2, 2), mgp = c(4, 2, 0))
    seqmtplot(seqact, with.legend = FALSE, main = "Mean time spent in every state", 
             border = NA, las = 2, cex.axis = .8)
  } else {
    num_clusters <- length(unique(clusters))
    plot_file <- sprintf("./results/mean-time-state-plot-%s-by-%s-clusters.png", model$name, num_clusters)
    png(file = plot_file,
        width = 8000, height = 6000, units = "px", res = 900)
    par(mar=c(6,6,2,2), mgp = c(4, 1, 0))
    if (num_clusters %% 2 == 0) {
      if (num_clusters == 2) {
        png(file = plot_file,
            width = 8000, height = 3000, units = "px", res = 900)
        par(mar=c(6,3,2,2), mgp = c(4, 1, 0))
      } else if (num_clusters == 4) {
        png(file = plot_file,
            width = 8000, height = 9000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
      } else if (num_clusters == 6) {
        png(file = plot_file,
            width = 8000, height = 10000, units = "px", res = 900)
        par(mar=c(6,5,2,2), mgp = c(4, 1, 0))
      }
      #par(mai=rep(0.5, 0.5))
    } else {
      par(mfrow = c(ceiling(num_clusters/2), 2), mar = c(6, 6, 2, 2), mgp = c(4, 1, 0))
    }
    seqmtplot(seqact, group = clusters, with.legend = FALSE,
               border = NA, las = 2, cex.axis = .6,
              main = "Mean time spent in each state", cex.main = 0.9,
               cex.lab = 0.8)
  }
  print(sprintf("Mean time plot saved here: %s", plot_file))
  #dev.off()
}


# entropy plot function for whole dataset
create_entropy_plot <- function(seqact, seqact2 = NULL) {
  if (!dir.exists("./results")) {
    dir.create("./results")
  }
  plot_file <- "./results/Entropy-Plot.png"
  png(file = plot_file, width = 9000, height = 6000, units = "px", res = 900)
  
  ylab_text <- sprintf("Entropy of patient trajectories (n=%s)", nrow(seqact))
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
    lines(entropy2, type = "l", lwd = 1.5, col = "blue", lty = 3)
    legend("topright", legend = c("GLP-1 RAs grouped", "Distinct GLP-1 RAs"),
           col = c("black", "blue"), lwd = 1.5, lty = c(1, 3), bty = "n")
  }
  print(sprintf("Entropy plot saved here: %s", plot_file))
  #dev.off()
}
