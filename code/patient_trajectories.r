# Load required libraries
library(TraMineR)
library(cluster)
library(WeightedCluster)
library(FactoMineR)
library(ade4)
library(RColorBrewer)
library(questionr)
library(descriptio)
library(dplyr)
library(purrr)
library(ggplot2)
library(seqhandbook)
library(JLutils)
library(fastcluster)
library(R6)
library(entropy)
library(seriation)
library(gplots)
library(dendextend)
library(caret)
library(clValid)
library(splitstackshape)
library(IRanges)


# sequence data
data_file <- "./data/distinct_GLP1_sequences.csv"
data_file_pooled <- "./data/GLP1_sequences.csv"
pat_details_file <- "./data/pat_level_info.csv"

alphabet <- c(0, 1, 2, 3, 4, 5, 6, 7)
palette <- brewer.pal(length(alphabet), 'Set2')
labels <- c("No GLP-1 RA or Other Anti-Obesity Medication",
            "Semaglutide",
            "Liraglutide",
            "Dulaglutide",
            "Tirzepatide",
            "Exenatide",
            "Other Anti-Obesity Medication",
            "GLP-1 RA & Other Anti-Obesity Medication")
scodes <- c("No medication",
            "Semaglutide",
            "Liraglutide",
            "Dulaglutide",
            "Tirzepatide",
            "Exenatide",
            "Other",
            "GLP-1 RA & Other")


alphabet_pooled <- c(0, 1, 2, 3)
palette_pooled <- brewer.pal(length(alphabet_pooled), 'Set2')
labels_pooled <- c("No GLP-1 RA or Other Anti-Obesity Medication",
                   "GLP-1 RA",
                   "Other Anti-Obesity Medication",
                   "GLP-1 RA & Other Anti-Obesity Medication")
scodes_pooled <- c("No medication",
                   "GLP-1 RA",
                   "Other",
                   "GLP-1 RA & Other")

# toggle to false to remove print statements
trace_flag <- TRUE
# toggle to TRUE if need visuals
heatmap_flag <- FALSE

# Load the data
trajact <- read.table(data_file, header = TRUE, sep = ",", row.names = 1)
colnames(trajact) <- paste('Day', -365:365, sep = ' ')
trajact_pooled <- read.table(data_file_pooled, header = TRUE, sep = ",", row.names = 1)
colnames(trajact_pooled) <- paste('Day', -365:365, sep = ' ')
pat_details <- read.table(pat_details_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
#str(pat_details)


create_stratified_cohort <- function(seqact, pat_details_file, sample_size = 0.1) {
  pat_details <- read.table(pat_details_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
  
  # Subset the data to include columns 1 (ID), 3 (Age), and the year part of column 5 (Index Date)
  pat_info <- pat_details %>%
    select(Gender, `Age at Index Date`, `Index Date`) %>%
    mutate(Year = as.numeric(format(as.Date(`Index Date`, format = "%Y-%m-%d"), "%Y"))) %>%
    select(Gender, `Age at Index Date`, Year)
  
  # get indicators
  indics3 <- seqi1epi(seqact)
  
  # Merge by row names
  stratified_vars <- merge(pat_info, indics3, by = "row.names")
  
  # Set the row names to the original row names
  row.names(stratified_vars) <- stratified_vars$Row.names
  stratified_vars <- stratified_vars[, -1]  # Remove the added 'Row.names' column
  
  # Print the combined data frame
  stratified_vars <- as.data.frame(scale(stratified_vars))
  
  # Remove the column 'epiNo medication'
  stratified_vars <- stratified_vars %>% select(-`epiNo medication`)
  
  # Create stratified sample
  set.seed(123)  # For reproducibility
  stratified_sample_indices <- stratified(stratified_vars, c("Gender","Age at Index Date","Year","epiOther","epiGLP-1 RA & Other"),
                                          sample_size, keep.rownames = TRUE)
  
  stratified_sample <- seqact[stratified_sample_indices$rn, ]
  
  return(stratified_sample)
}


# create sequence for info
seqact_pooled_all <- seqdef(trajact_pooled, cpal = palette_pooled,
                            alphabet = alphabet_pooled,
                            states = scodes_pooled,
                            labels = labels_pooled)
seqact_all <- seqdef(trajact, cpal = palette,
                     alphabet = alphabet,
                     states = scodes,
                     labels = labels,
                     xstep = round(ncol(trajact) / 43, 0),
                     tick.last = TRUE)


rand_entropy_all <- 0.33
rand_entropy_rhs <- 0.628

# Define the class
TrajectoryAnalysis <- R6Class(
  "TrajectoryAnalysis",
  public = list(
    trace_flag = NULL,  # option to print additional information at various stages
    name = NULL,  # option to add a name for the model
    trajact = NULL,  # to store trajectories information
    seqact = NULL,  # to store the sequence object of the trajectories
    palette = NULL,  # colour palette to use for visualisations
    alphabet = NULL,  # alphabet of original trajectory file
    labels = NULL,  # labels for sequence object
    scodes = NULL,  # scodes for sequence object
    cost_method = NULL,  # substitution cost method
    time_varies = NULL,  # substitution cost matrix - for DHD
    con_val = NULL,  # cval in creating substitution cost matrix
    sub_costs = NULL,  
    seqdist_method = NULL,  # method to create distance matrix
    seqdist_indel = NULL,  # indel if using for distance matrix (e.g. OM)
    seqdist_sm = NULL,  # substitution cost matrix result
    seqdist_oth_args = NULL,  # other seqdist arguments
    start_time = NULL,  # start time for creating model
    end_time = NULL,  # end time for creating model
    time_taken = NULL,  # time taken when model is run
    cluster_method = NULL,  # cluster method (default is Ward)
    clust = NULL,  # to store hierarchical clustering object once run
    olo = NULL,  # boolean to speficy whether to optimise leaf ordering
    results = NULL,  # to store a data frame of results
    tpow = NULL,  # parameter for OMspell distance matrix
    expcost = NULL,  # parameter for OMspell distance matrix
    range_1 = NULL,  # range 1 to manually adjust
    range_1_val = NULL,  # value to replace all those in range 1
    state_to_adjust = NULL,  # the state to adjust manually (defaults to first state in scodes list)
    state_cost = NULL,  # amount to adjust the state by (multiply by)
    heatmap_flag = NULL,  # flag to indicate whether to save down heatmap of clustering
    seq_norm = NULL,  # for seqdist()
    rand_entropy = NULL,  # option to incorporate random entropy score into the results
    rand_entropy_rhs = NULL,  # option to incorporate excluding pre-index date random entropy into results
    dissim = NULL,  # distance matrix once run
    nu = NULL,  # for TWED dissim
    h = NULL,  # for TWED dissim
    pi_dissim = NULL,  # for making the distance matrix from index only
    glp1_sub_vals = NULL,  # resetting the values for GLP1s when using distinct
    breaks = NULL,  # for CHI2 seqdist
    step = NULL,  # for CHI2 seqdist
    lag = NULL,  # for TRATE and FUTURE substitution cost matrices
    
    initialize = function(seqact, palette, alphabet, labels, scodes, name, trace_flag = FALSE,
                          cost_method = "CONSTANT", time_varies = FALSE, con_val = 2,
                          seqdist_method = "HAM", seqdist_indel = "auto", seqdist_sm = NULL,
                          cluster_method = "ward.D2", expcost = 0.5, tpow = 1.0, seq_norm = "auto",
                          range_1 = NULL, range_1_val = NULL, state_to_adjust = scodes[1],
                          state_cost = NULL, heatmap_flag = TRUE, olo = TRUE, 
                          rand_entropy = NA, rand_entropy_rhs = NA, nu = NA, h = 0.5,
                          pi_dissim = TRUE, glp1_sub_vals = NULL, breaks = NULL, step = 1,
                          lag = 1, ...) {
      self$trajact <- trajact
      self$seqact <- seqact
      self$palette <- palette
      self$alphabet <- alphabet
      self$labels <- labels
      self$scodes <- scodes
      self$trace_flag <- trace_flag
      self$name <- name
      self$cost_method <- cost_method
      self$time_varies <- time_varies
      self$con_val <- con_val
      self$seqdist_method <- seqdist_method
      self$seqdist_indel <- seqdist_indel
      self$seqdist_sm <- seqdist_sm
      self$cluster_method <- cluster_method
      self$seqdist_oth_args <- list(...)
      self$results <- data.frame()
      self$tpow <- tpow
      self$expcost <- expcost
      self$seq_norm <- seq_norm
      self$range_1 <- range_1
      self$range_1_val <- range_1_val
      self$state_to_adjust <- state_to_adjust
      self$state_cost <- state_cost
      self$heatmap_flag <- heatmap_flag
      self$olo <- olo
      self$rand_entropy <- rand_entropy
      self$rand_entropy_rhs <- rand_entropy_rhs
      self$nu <- nu
      self$h <- h
      self$pi_dissim <- pi_dissim
      self$glp1_sub_vals <- glp1_sub_vals
      self$breaks <- breaks
      self$step <- step
      self$lag <- lag
    },
    
    trace = function(info) {
      if (self$trace_flag) {
        print(info)
      }
    },
    
    # Create the substitution cost matrix
    set_costs = function() {
      seq = if(self$pi_dissim) {
        self$seqact[, 366:731]
      } else {
        self$seqact
      }
      if (self$cost_method %in% c("INDELSLOG", "INDELS", "FEATURES", "FUTURE")) {
        self$sub_costs <- seqcost(seq, method = self$cost_method,
                                  time.varying = self$time_varies,
                                  lag = self$lag)
      } else {
        self$sub_costs <- seqcost(seq, method = self$cost_method,
                                  time.varying = self$time_varies,
                                  cval = self$con_val,
                                  lag = self$lag)
        if (is.numeric(self$sub_costs$indel)) {
          if (self$seqdist_indel == "auto") {
            self$seqdist_indel <- self$sub_costs$indel
          }
        }
      }
    },
    
    # Check whether an attribute is not empty (return boolean)
    check_not_empty = function(attr) {
      if ((is.vector(attr) && length(attr) > 0) ||
          (is.list(attr) && length(attr) > 0)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    
    # method to adjust event
    adjust_event_subs = function(mat) {
      mat[self$state_to_adjust,][mat[self$state_to_adjust,]>0] <- mat[self$state_to_adjust,][mat[self$state_to_adjust,]>0] * self$state_cost
      mat[,self$state_to_adjust][mat[,self$state_to_adjust]>0] <- mat[,self$state_to_adjust][mat[,self$state_to_adjust]>0] * self$state_cost
      
      return(mat)
    },
    
    adjust_distinct_glp1s = function(mat) {
      GLP1_states <- c("Semaglutide",
                       "Liraglutide",
                       "Dulaglutide",
                       "Tirzepatide",
                       "Exenatide")
      for (state1 in GLP1_states) {
        for (state2 in GLP1_states) {
          if (state1 != state2) {
            mat[state1, state2] <- self$glp1_sub_vals
            mat[state2, state1] <- self$glp1_sub_vals
          }
        }
      }
      return(mat)
    },
    
    # method to trigger adjusting importance of first event (no meds)
    adjust_events = function() {
      seq = if(self$pi_dissim) {
        self$seqact[, 366:731]
      } else {
        self$seqact
      }
      if (isTRUE(self$time_varies)) {
        for (i in 1:ncol(seq)) {
          mat <- self$sub_costs$sm[,,i]
          if (!is.null(self$state_cost)) {
            self$sub_costs$sm[,,i] <- self$adjust_event_subs(mat)
          }
          
          if (!is.null(self$glp1_sub_vals)) {
            self$sub_costs$sm[,,i] <- self$adjust_distinct_glp1s(mat)
          }
        }
      } else {
        mat <- self$sub_costs$sm
        if(!is.null(self$state_cost)) {
          self$sub_costs$sm <- self$adjust_event_subs(mat)
        }
        
        if (!is.null(self$glp1_sub_vals)) {
          self$sub_costs$sm <- self$adjust_distinct_glp1s(mat)
        }
      }
      if(!is.null(self$state_cost)) {
        self$trace(sprintf("Substitution cost matrix adjusted for %s x %s", self$scodes[1], self$state_cost))
      }
      if(!is.null(self$glp1_sub_vals)) {
        self$trace(sprintf("Substitution cost matrix adjust for between GLP-1 RA states to %s", self$glp1_sub_vals))
      }
      
    },
    
    # create entropy measure
    entropy_measure = function(rand = FALSE, rhs = FALSE, glp1_pooled = FALSE) {
      # Obtain the order from the clustering
      order <- self$clust$order
      
      # Order the sequences
      ordered_seq <- self$seqact[order,]
      
      # Update sequences if glp1_pooled is TRUE
      if (glp1_pooled == TRUE) {
        glp1_values <- c("Semaglutide", "Liraglutide", "Dulaglutide", "Tirzepatide", "Exenatide")
        ordered_seq <- apply(ordered_seq, 2, function(x) ifelse(x %in% glp1_values, "GLP-1 RA", x))
      }
      
      # Extract the matrix of sequences
      seq_matrix <- as.matrix(ordered_seq)
      
      # if random order
      if (rand == TRUE) {
        rand_rows <- sample(nrow(seq_matrix))
        seq_matrix <- seq_matrix[rand_rows,]
      }
      
      # if RHS only
      if (isTRUE(rhs)) {
        seq_matrix <- seq_matrix[, 366:731]
      }
      # Define the grid size
      grid_size <- 3
      
      # Initialize a matrix to store entropy values
      n_rows <- nrow(seq_matrix)
      n_cols <- ncol(seq_matrix)
      entropy_matrix <- matrix(NA, nrow = ceiling(n_rows / grid_size), ncol = ceiling(n_cols / grid_size))
      
      # Compute entropy for each non-overlapping 3x3 grid, handling edge cases
      for (i in seq(1, n_rows, by = grid_size)) {
        for (j in seq(1, n_cols, by = grid_size)) {
          # Determine the bounds of the sub-matrix
          end_row <- min(i + grid_size - 1, n_rows)
          end_col <- min(j + grid_size - 1, n_cols)
          
          # Extract the sub-matrix
          sub_matrix <- seq_matrix[i:end_row, j:end_col]
          sub_matrix_flat <- as.vector(sub_matrix) # Flatten the sub-matrix
          
          # Calculate probability distribution
          prob <- table(sub_matrix_flat) / length(sub_matrix_flat)
          
          # Compute entropy and store in the corresponding position
          row_idx <- ceiling(i / grid_size)
          col_idx <- ceiling(j / grid_size)
          
          if (row_idx <= nrow(entropy_matrix) && col_idx <= ncol(entropy_matrix)) {
            entropy_matrix[row_idx, col_idx] <- entropy::entropy.empirical(prob, unit = "log2")
          }
        }
      }
      
      # Compute the mean entropy score of the matrix
      mean_entropy <- mean(entropy_matrix)
      return(mean_entropy)
    },
    
    dend_heatmap = function(k = 10) {
      seq = if(self$pi_dissim) {
        self$seqact[, 366:731]
      } else {
        self$seqact
      }
      clusters <- cutree(self$clust, k = k)
      order <- self$clust$order
      sorted_clusters <- clusters[order]
      
      cluster_dir <- sprintf("./results/%s", self$name)
      if (!dir.exists(cluster_dir)) {
        dir.create(cluster_dir)
      }
      heatmap_start <- Sys.time()
      
      self$trace("Directory for results created/exists for saving down dendrogram cluster annotated heatmap")
      plot_file <- sprintf("%s/%s_Treatment-Pattern-Sequences-Ordered_dendrogram-%s-clusters.png", cluster_dir, self$name, k)
      png(file = plot_file,
          width = 9000, height = 10000, units = "px", res = 900)
      par(mar = c(5, 1, 1, 1))
      seq_heatmap_alt <- function (seq, tree, sorted_clusters, with.missing=FALSE, ...) {
        if (class(tree)!="dendrogram") tree <- as.dendrogram(tree)
        
        # Order the sequences
        mat <- seq
        for (i in 1:length(self$seqact)){
          mat[mat[,i]=="%",i] <- NA
          mat[,i] <- as.numeric(mat[,i])
        }
        mat <- as.matrix(mat)
        
        cluster_lines = c()
        last_cluster <- sorted_clusters[1]
        cluster_start <- 1
        for (i in 2:length(sorted_clusters)) {
          if (sorted_clusters[i] != last_cluster) {
            cluster_lines[i] <- i
            last_cluster<- sorted_clusters[i]
            cluster_start <- i
          }
        }
        yheading <- sprintf("Cohort: %s patients, ordered", nrow(self$seqact))
        title <- sprintf("Treatment Pattern Sequences and Dendrogram - %s Clusters", k)
        
        #self$trace(cluster_lines)
        heatmap(mat, tree, col=self$palette, scale = 'none', Colv = NA,
                main = title, las = 2, cex.axis = .8,
                ylab = yheading, xlab = NA,
                labRow = NA, margins = c(6, 2),
                add.expr = {abline(v = 366, col = "black", lty = 2, lwd = 1.5);
                  abline(h = cluster_lines, col = rgb(0, 0, 0, 0.6), lty = 2, lwd = 1.5);
                  text(x = 355, y = nrow(self$seqact) * 0.09, srt = 90, "Index Date");
                  last_cluster <- sorted_clusters[1]
                  cluster_start <- 1
                  for (i in 1:length(sorted_clusters)) {
                    if (sorted_clusters[i] != last_cluster) {
                      text(x = ncol(self$seqact) / 5, y = (cluster_start + i - 1) / 2 , labels = paste("Cluster", last_cluster), pos = 2, col = "black")
                      last_cluster<- sorted_clusters[i]
                      cluster_start <- i
                    }
                  }
                  text(x = ncol(self$seqact) / 5, y = (cluster_start + length(sorted_clusters)) / 2 , labels = paste("Cluster", last_cluster), pos = 2, col = "black")
                })
        
        if(length(self$alphabet) > 4) {
          ncol_leg <- 3
        } else {
          ncol_leg <- 2
        }
        
        legend("bottom", title = "Days in proximity to Index Date (first GLP-1 RA prescription)",
               title.cex = 1.05, ncol = ncol_leg, legend = self$scodes, fill = self$palette,
               bty = "n", x.intersp = 1.5, cex = 0.8)
        
      }
      
      seq_heatmap_alt(self$seqact, self$clust, sorted_clusters)
      
      dev.off()
      heatmap_end <- Sys.time()
      heatmap_time_taken <- round(difftime(heatmap_end, heatmap_start, units = "secs"), 3)
      self$trace(sprintf("Dendrogram cluster annotated heatmap saved for model % s: %s seconds", self$name, heatmap_time_taken))
    },
    
    # create visuals with cluster divisions mapped out (no dendrogram)
    cluster_heatmap = function(name= NA, k=3) {
      clusters <- cutree(self$clust, k = k)
      
      order <- self$clust$order
      sorted_clusters <- clusters[order]
      
      cluster_dir = sprintf("./results/%s", self$name)
      if (!dir.exists(cluster_dir)) {
        dir.create(cluster_dir)
      }
      heatmap_start <- Sys.time()
      self$trace("Directory for results created/exists for saving down cluster annotated heatmap")
      plot_file <- sprintf("%s/%s_Treatment-Pattern-Sequences-Ordered_%s-clusters.png", cluster_dir, self$name, k)
      png(file = plot_file,
          width = 9000, height = 10000, units = "px", res = 900)
      par(mar = c(10, 5, 2, 2), mgp = c(4, 1, 0))
      
      title <- if (!(is.na(name))) {
        sprintf("%s Treatment Pattern Sequences by Cluster Order", name)
      } else {
        "Treatment Pattern Sequences by Cluster Order"
      }
      
      seqIplot(self$seqact, with.legend = FALSE, sortv = cutree.order(self$clust, k = nrow(self$seqact)),
               main = title, xlab = "Days in proximity to Index Date",
               las = 2, cex.axis = .8,
               xtstep = round(ncol(self$seqact) / 43, 0), tick.last = TRUE)
      # Add a vertical line and text label
      usr <- par("usr")
      segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4], col = "black", lty = 2, lwd = 2)
      text(x = 355, y = nrow(self$seqact) * 0.09, srt = 90, "Index Date")
      
      last_cluster <- sorted_clusters[1]
      cluster_start <- 1
      for (i in 2:length(sorted_clusters)) {
        if (sorted_clusters[i] != last_cluster) {
          abline(h = i - 0.5, col = rgb(0, 0, 0, 0.7), lty = 2, lwd = 2)
          text(x = ncol(self$seqact) / 5, y = (cluster_start + i - 1) / 2 , labels = paste("Cluster", last_cluster), pos = 2, col = "black")
          last_cluster<- sorted_clusters[i]
          cluster_start <- i
        }
      }
      
      text(x = ncol(self$seqact) / 5, y = (cluster_start + length(sorted_clusters)) / 2 , labels = paste("Cluster", last_cluster), pos = 2, col = "black")
      
      # Add a legend
      par(xpd = TRUE)
      if(length(self$alphabet) > 4) {
        ncol_leg <- 3
      } else {
        ncol_leg <- 2
      }
      legend("bottom", inset = c(-0.18, -0.18), ncol = ncol_leg, legend = self$labels, fill = self$palette,
             bty = "n", x.intersp = 0.2, cex = 0.8)
      dev.off()
      heatmap_end <- Sys.time()
      heatmap_time_taken <- round(difftime(heatmap_end, heatmap_start, units = "secs"), 3)
      self$trace(sprintf("Cluster annotated heatmap saved for model % s: %s seconds", self$name, heatmap_time_taken))
      
    },
    
    # create and save down dendrogram
    create_dendrogram = function() {
      cluster_dir = sprintf("./results/%s", self$name)
      if (!dir.exists(cluster_dir)) {
        dir.create(cluster_dir)
      }
      self$trace(sprintf("Creating dendrogram of %s...", self$name))
      plot_file <- sprintf("%s/%s_dendrogram.png", cluster_dir, self$name)
      png(file = plot_file,
          width = 12000, height = 8000, units = "px", res = 900)
      par(mar = c(2, 2, 2, 2))
      plot(self$clust, labels = FALSE,
           main = "Dendrogram of agglomerative clustering (Ward's method)")
      dev.off()
      self$trace(sprintf("Dendrogram saved: %s", plot_file))
    },
    
    # Run the sequence analysis
    run_analysis = function() {
      seq = if(self$pi_dissim) {
        self$seqact[, 366:731]
      } else {
        self$seqact
      }
      self$start_time <- Sys.time()
      self$trace(sprintf("Starting model %s", self$name))
      
      # Compute the substitution costs
      if (!(self$cost_method %in% c("CHI2", "EUCLID"))) {
        self$set_costs()
        self$trace("Substitution cost matrix built")
        if (self$check_not_empty(self$range_1)) {
          # Get the subarray based on the range
          sub_array <- self$sub_costs$sm[,,self$range_1]
          
          # Modify the subarray where the condition is met
          sub_array[sub_array > 0] <- self$range_1_val
          
          # Assign the modified subarray back to the original array
          self$sub_costs$sm[,,self$range_1] <- sub_array
          
          # Log the trace information
          self$trace(sprintf("Value %s input for time units in range", self$range_1_val))
        }
        # check if manually editing values
        if (!is.null(self$state_cost) |
            !is.null(self$glp1_sub_vals)) {
          self$adjust_events()
          self$trace("Manual state costs used...")
        }
      } else {
        self$sub_costs$sm <- NULL
      }
      
      
      # Compute the distance matrix and time it
      self$trace("Building Distance Matrix")
      dist_start_time <- Sys.time()
      
      self$dissim <- seqdist(seq, method = self$seqdist_method, sm = self$sub_costs$sm, 
                             indel = self$seqdist_indel, expcost = self$expcost, tpow = self$tpow,
                             norm = self$seq_norm, full.matrix=FALSE, nu = self$nu, h = self$h,
                             breaks = self$breaks, step = self$step)
      dist_end_time <- Sys.time()
      dist_time_taken <- round(difftime(dist_end_time, dist_start_time, units = "secs"), 3)
      self$trace(sprintf("Distance Matrix built: %s seconds", dist_time_taken))
      
      # Perform hierarchical clustering
      self$trace("Starting clustering model")
      clust_start_time <- Sys.time()
      self$clust <- as.dist(self$dissim) %>% hclust(method = self$cluster_method)
      
      # Record the end time and compute the total time taken
      clust_end <- Sys.time()
      
      clust_time_taken <- round(difftime(clust_end, clust_start_time, units = "secs"), 3)
      
      self$trace(sprintf("Finished clustering model: %s seconds", clust_time_taken))
      
      
      # optimal leaf ordering
      if (isTRUE(self$olo)) {
        self$trace("Starting optimal leaf ordering...")
        olo_start <- Sys.time()
        self$clust <- reorder(self$clust, as.dist(self$dissim), method = "GW")
        self$end_time <- Sys.time()
        olo_duration <- round(difftime(self$end_time, olo_start, units = "secs"), 3)
        self$trace(sprintf("Optimal leaf ordering completed: %s seconds", olo_duration)) 
      } else {
        olo_duration <- NA
        self$end_time <- Sys.time()
      }
      
      
      self$time_taken <- round(difftime(self$end_time, self$start_time, units = "secs"), 3)
      self$trace(sprintf("Total model time: %s seconds", self$time_taken))
      
      # Plot the sequences with clustering
      if (isTRUE(self$heatmap_flag)) {
        cluster_dir = sprintf("./results/%s", self$name)
        if (!dir.exists(cluster_dir)) {
          dir.create(cluster_dir)
        }
        heatmap_start <- Sys.time()
        self$trace("Directory for results created/exists for saving down heatmap")
        plot_file <- sprintf("%s/%s_Treatment-Pattern-Sequences-Ordered.png", cluster_dir, self$name)
        png(file = plot_file,
            width = 9000, height = 10000, units = "px", res = 900)
        par(mar = c(10, 5, 2, 2))
        
        seqIplot(self$seqact, with.legend = FALSE, sortv = cutree.order(self$clust, k = nrow(seqact)),
                 main = "Treatment Pattern Sequences by Cluster Order", las = 2, cex.axis = .8,
                 xtstep = round(ncol(self$seqact) / 43, 0), tick.last = TRUE)
        # Add a vertical line and text label
        usr <- par("usr")
        segments(x0 = 366, y0 = usr[3], x1 = 366, y1 = usr[4], col = "black", lty = 2, lwd = 2)
        text(x = 355, y = nrow(self$seqact) * 0.03, srt = 90, "Index")
        # Add a legend
        par(xpd = TRUE)
        legend("bottom", inset = c(-0.2, -0.2), ncol = 2, legend = self$labels, fill = self$palette,
               bty = "n", x.intersp = 0.2, cex = 0.8)
        dev.off()
        heatmap_end <- Sys.time()
        heatmap_time_taken <- round(difftime(heatmap_end, heatmap_start, units = "secs"), 3)
        self$trace(sprintf("Heatmap saved for model % s: %s seconds", self$name, heatmap_time_taken))
      }
      
      # get entropy result
      entropy_res <- self$entropy_measure()
      
      # get entropy result RHS only
      entropy_rhs <- self$entropy_measure(rhs = TRUE)
      
      if (nrow(self$seqact) < 5000) {
        # get ASW result
        wardRange <- as.clustrange(self$clust, diss=self$dissim, nclusters = 10)
        asww_value <- summary(wardRange, max.rank = 1)["ASWw",2]
        asww_num_clusters <- summary(wardRange, max.rank = 1)["ASWw",1]
        
        # CHsq
        CHsq_value <- summary(wardRange, max.rank = 1)["CHsq",2]
        CHsq_clusters <- summary(wardRange, max.rank = 1)["CHsq",1]
        
      } else {
        asww_value <- NA
        asww_num_clusters <- NA
        CHsq_value <- NA
        CHsq_clusters <- NA
      }
      
      
      # Handle NULL values by replacing them with appropriate placeholders
      Range1 <- if (is.null(self$range_1)) "NA" else as.character(sprintf("%s:%s", min(self$range_1), max(self$range_1)))
      Range1_val <- if (is.null(self$range_1_val)) NA else self$range_1_val
      
      # Store the results
      self$results <- data.frame(
        ModelName = self$name,
        TimeTaken = as.numeric(self$time_taken),
        DistMatrixTime = as.numeric(dist_time_taken),
        ClusteringTime = as.numeric(clust_time_taken),
        OLOTime = as.numeric(olo_duration),
        ASWw = round(asww_value, 3),
        ASWw_num_clusters = asww_num_clusters,
        CHsq_value = round(CHsq_value, 3),
        CHsq_clusters = CHsq_clusters,
        Entropy = round(entropy_res, 3),
        EntropyRHS = round(entropy_rhs, 3),
        EntropyRandom = self$rand_entropy,
        EntropyRandomRHS = self$rand_entropy_rhs,
        SeqDistMethod = self$seqdist_method,
        ClustMethod = self$cluster_method,
        SubsMethod = self$cost_method,
        ConstantValue = self$con_val,
        Range1 = Range1,
        Range1_val = Range1_val,
		    Indel = self$seqdist_indel,
        NPatients = nrow(self$seqact)
      )
      self$trace(sprintf("Results saved for %s", self$name))
    }
  )
)
