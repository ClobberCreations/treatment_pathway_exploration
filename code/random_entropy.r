# libraries
library(entropy)
library(parallel)
library(doParallel)
library(foreach)


# create entropy measure
rand_entropy_measure <- function(seq, rhs = FALSE) {
    # Extract the matrix of sequences
    seq_matrix <- as.matrix(seq)
    
    rand_rows <- sample(nrow(seq_matrix))
    seq_matrix <- seq_matrix[rand_rows,]
    
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
}