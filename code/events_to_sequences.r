### Code to transform the event log into a wide format sequence table ###

# Load libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(data.table)


sequence_format <- function(dataset, col_names, dataset_header = TRUE,
                            date_cols = NA, combine_GLP1s = TRUE) {
  # Step 1: Read in the CSV file, rename columns and specify datetime fields
  df <- read.csv(dataset, header=dataset_header)
  colnames(df) <- col_names
  
  # Turn date cols into datetime type
  if (length(date_cols) > 0) {
    df[date_cols] <- lapply(df[date_cols], as.Date, format = "%Y-%m-%d")
  }
  
  # Step 2: Remove rows with no events for GLP1s or Other_medications
  df <- df %>%
    filter(!(GLP1s == 0 & Other_medications == 0))
  
  # Step 3: Split the data frame into two new data frames - other and GLP1s
  
  # other to include all rows where another anti-obesity medication is recorded
  other <- df %>%
    filter(Other_medications == 1)
  # GLP1s to include all rows where a GLP1 is recorded
  GLP1s <- df %>%
    filter(GLP1s == 1)
  
  # Step 4: Add column ‘Event’ to both new data frames
  other <- other %>%
    mutate(Event = "Other_Med")
  
  GLP1s <- GLP1s %>%
    mutate(Event = "GLP1")
  
  # Step 5: Add column Duration to other df
  other <- other %>%
    mutate(Duration = Other_medications_days_supply)
  
  # Step 6: Add column Duration to GLP1s and include GLP1 column to highlight
  # which GLP1 medication
  GLP1s <- GLP1s %>%
    mutate(Duration = GLP1_days_supply)
  
  # list of GLP1 medication column names - matches ARA platform for ease
  GLP1_medications <- c("DULAGLUTIDE",                
                        "EXENATIDE  (BYDUREON, BYETTA)",
                        "SEMAGLUTIDE (OZEMPIC)",
                        "LIRAGLUTIDE (VICTOZA)",
                        "TIRZEPATIDE (MOUNJARO)",
                        "LIXISENATIDE (ADLYXIN)",
                        "LIRAGLUTIDE (SAXENDA)",
                        "ALBIGLUTIDE (TANZEUM)",
                        "SEMAGLUTIDE (WEGOVY)")
  
  # Add column for which GLP1
  GLP1s$GLP1 <- colnames(GLP1s[, GLP1_medications])[apply(GLP1s[, GLP1_medications], 1, which.max)]
  
  # Step 7: Load in csv containing the prescription rules
  presc_rules <- read.csv("./data/Prescription Duration Rules Table_data_20240624.csv")
  
  # Step 8: Impute Duration column for GLP1s data frame if NA or 0
  # (other contains no NAs)
  GLP1s <- GLP1s %>%
    mutate(Duration = ifelse(is.na(Duration) | Duration == 0,
                             presc_rules$Duration[match(GLP1, presc_rules$Medications)], Duration))
  
  # update to only show generic name
  GLP1s$Event <- sub(" \\(.*\\)", "", GLP1s$GLP1)
  
  #print(unique(GLP1s$Event))
  
  # Step 9: Extend both other and GLP1s dfs so that there is a row and EventDate
  # for every day in the Duration column
  
  # define function to extend dates
  extend_dates <- function(df) {
    df <- as.data.table(df)
    df <- df[, .(EventDate = seq(EventDate, by = "day", length.out = Duration),
                 PatId = PatId,
                 Event = Event,
                 IndexDate = IndexDate),
             by = 1:nrow(df)]
    df[, nrow := NULL]
  }
  
  # apply function to both data frames
  other <- extend_dates(other)
  GLP1s <- extend_dates(GLP1s)
  
  # Step 10: Remove duplicate rows
  other <- distinct(other)
  GLP1s <- distinct(GLP1s)
  
  # Step 11: Create ‘Day’ column for both other and GLP1s data frames which
  # provides the difference in days between EventDate and IndexDate
  other <- other %>%
    mutate(Day = as.integer(difftime(EventDate, IndexDate, units = "days")))
  
  GLP1s <- GLP1s %>%
    mutate(Day = as.integer(difftime(EventDate, IndexDate, units = "days")))
  
  # Step 12: Merge other and GLP1s data frames on Day and PatId columns
  combined <- full_join(other, GLP1s, by = c("PatId", "Day"))
  
  # Step 13: Create new column called “Activity” in combined data frame
  if (combine_GLP1s) {
    combined <- combined %>%
      mutate(Activity = case_when(
        is.na(Event.x) & !is.na(Event.y) ~ 1,
        !is.na(Event.x) & is.na(Event.y) ~ 2,
        !is.na(Event.x) & !is.na(Event.y) ~ 3,
        TRUE ~ 0
      ))
  } else {
    combined <- combined %>%
      mutate(Activity = case_when(
        is.na(Event.x) & Event.y == "SEMAGLUTIDE" ~ 1,
        is.na(Event.x) & Event.y == "LIRAGLUTIDE" ~ 2,
        is.na(Event.x) & Event.y == "DULAGLUTIDE" ~ 3,
        is.na(Event.x) & Event.y == "TIRZEPATIDE" ~ 4,
        is.na(Event.x) & Event.y == "EXENATIDE " ~ 5,
        !is.na(Event.x) & is.na(Event.y) ~ 6,
        !is.na(Event.x) & !is.na(Event.y) ~ 7,
        TRUE ~ 0
      ))
  }
  
  # Step 14: Create a new empty data frame called ‘seq_df’ which is the required
  # structure for TraMineR library - with a column for -1 year to 1 year
  # post-IndexDate
  unique_pat_ids <- unique(combined$PatId)
  # create empty data frame of zero values
  seq_df <- data.frame(matrix(0, nrow = length(unique_pat_ids), ncol = 731))
  rownames(seq_df) <- unique_pat_ids
  colnames(seq_df) <- -365:365
  
  # Step 15: Use the combined df values to fill in the seq_df
  
  # turn to data.table to run a bit faster than data frame
  combined_dt <- as.data.table(combined)
  seq_df_dt <- as.data.table(seq_df, keep.rownames = "PatId")
  
  # set key for both data tables
  setkey(combined_dt, PatId)
  setkey(seq_df_dt, PatId)
  
  # transform seq_df_dt data table to long format to help transfer values
  seq_df_long <- melt(seq_df_dt, id.vars = "PatId", variable.name = "Day",
                      value.name = "Activity")
  # order by PatId
  seq_df_long <- seq_df_long[order(seq_df_long$PatId),]
  # ensure Day field is interpreted as an integer
  seq_df_long$Day <- as.integer(as.character(seq_df_long$Day))
  # ensure PatId field is interpreted as a character
  combined_dt$PatId <- as.character(combined_dt$PatId)
  # transfer Activity values available in combined_dt to seq_df_long
  seq_df_long <- seq_df_long[combined_dt, Activity := i.Activity, on = .(PatId, Day)]
  
  # create seq_df_final by casting seq_df_long back to wide format
  seq_df_final <- dcast(seq_df_long, PatId ~ Day, value.var = "Activity")
  # convert PatId column to the row names
  seq_df_final <- seq_df_final %>%
    column_to_rownames(var = "PatId")
  # add prefix to the column names
  colnames(seq_df_final)<- paste('Day', colnames(seq_df_final), sep = '_')
  # preview changes
  head(seq_df_final)
  
  # add index date marker
  # seq_df_final$Day_0 <- 9
  
  # file name depending on whether distinct GLP1s or not
  if (combine_GLP1s) {
    seq_file_path <- "./data/GLP1_sequences.csv"
  } else {
    seq_file_path <- "./data/distinct_GLP1_sequences.csv"
  }
  
  # Step 16: Save down as a csv called “GLP1_sequences.csv”
  write.csv(seq_df_final, seq_file_path, row.names = TRUE)
  
  # print to advise completed
  sprintf("Wide format sequence saved down here: %s", seq_file_path)
}
