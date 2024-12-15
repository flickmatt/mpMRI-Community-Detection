# This would create events for the RNA seq data

# Load tidyverse library
library(tidyverse)
library(tibble)
library(gsubfn)
library(dendextend)

#------------------------------------------------------------------------------#
# Perform hierarchical clustering and return the events
#------------------------------------------------------------------------------#
get_img_events <- function(img_df) {
  
  # Dissimilarity matrix based on euclidean distance
  d <- dist(img_df, method = "euclidean")
  # Hierarchical clustering using median linkage
  hc1 <- hclust(d, method = "median" )
  # expression rate dendrogram
  dend <- as.dendrogram(hc1)
  # Return all the nodes in the tree with corresponding samples
  events <- partition_leaves(dend)
  
  return(events)
}

# get imaging
img <- data.frame(c(master[19435:19441], master[19373:19418]))
img$Iavarone_ID <- master$Iavarone_ID

# All ----
setwd("/Users/m254284/Desktop/Neo4j_Projects/Img_events/All")
mri_annotation <- data.frame(master$MRI.contrast.enhancing.annotation, master$Iavarone_ID)
colnames(mri_annotation) <- c("MRI.contrast.enhancing.annotation", "Iavarone_ID")
img_matrix_1 <- img
img_matrix_1 <- merge(img_matrix_1, mri_annotation, by = "Iavarone_ID")
img_matrix_1$MRI.contrast.enhancing.annotation <- NULL

img_matrix_1 <- na.omit(img_matrix_1)
ids <- img_matrix_1$Iavarone_ID
img_matrix_1$Iavarone_ID <- NULL
img_matrix <- data.frame(t(img_matrix_1))
img_Symbol <- rownames(img_matrix)
img_matrix <- data.frame(lapply(img_matrix, as.numeric))
colnames(img_matrix) <- ids
row.names(img_matrix) <- NULL
no_of_imgs <- dim(img_matrix)[1]
img_names <- img_Symbol

# For each img we make the tree
for (img_no in seq(1,no_of_imgs)) {
  
  # Print the img name and the number
  print(paste("Working on img no.", img_no, ":", img_names[img_no]))
  
  # Get the img  data
  ge_df_matrix <- as.matrix(img_matrix[1:length(img_matrix)])
  rownames(ge_df_matrix) <- img_Symbol
  
  # Transpose the above matrix
  transposed_ge_df <- t(ge_df_matrix)
  transposed_ge_df <- as.data.frame(transposed_ge_df)
  
  # Get the data for the img
  img_df <- as.data.frame(transposed_ge_df[img_no])
  # head(img_df)
  img_name = colnames(img_df)
  
  # Get the img events
  img_events <- get_img_events(img_df = img_df)
  
  # Create an empty data frame to append the result of each img tree events
  df <- data.frame(matrix(ncol = 7, nrow = 0))
  
  for (i in seq(1, length(img_events))) {
    
    # Get event characteristics
    event <- i
    # Get event patients
    unlisted_patients <- unlist(img_events[i])
    patients <- paste(unlisted_patients, collapse=",")
    # Get median img 
    patients_ge <- img_df %>% filter(row.names(img_df) %in% unlisted_patients)
    median_ge <- median(patients_ge[[1]])
    # Get the number of patients and if there is only one patient, its a leaf
    no_of_patients <- length(unlisted_patients)
    if (no_of_patients == 1) {
      leaf_status <- 1
    } else {
      leaf_status <- 0
    }
    tissue_region <- "All"
    # Append the above data frame to df
    df <- rbind(df, c(event, patients, median_ge, no_of_patients, img_name, 'img_value', leaf_status, tissue_region))
  }
  
  # Provide column names
  colnames(df) <- c('event', 'samples', 'median_signal', 'no_of_patients', 'img_name', 'molecular_data', 'leaf_status', 'tissue_region')
  
  file_name = paste(img_name, ".csv", sep="")
  write.csv(df, file=file_name, row.names = FALSE)
}


# NE ----
#NE
setwd("/Users/m254284/Desktop/Neo4j_Projects/Img_events/NE")
tissue_region <- "NE"
mri_annotation <- data.frame(master$MRI.contrast.enhancing.annotation, master$Iavarone_ID)
colnames(mri_annotation) <- c("MRI.contrast.enhancing.annotation", "Iavarone_ID")
img_matrix_1 <- img
img_matrix_1 <- merge(img_matrix_1, mri_annotation, by = "Iavarone_ID")
img_matrix_1 <- subset(img_matrix_1, MRI.contrast.enhancing.annotation == tissue_region)
img_matrix_1$MRI.contrast.enhancing.annotation <- NULL

img_matrix_1 <- na.omit(img_matrix_1)
matched_pts <- data.frame(all_nodes_deg$Iavarone_ID)
colnames(matched_pts) <- "Iavarone_ID"
img_matrix_1 <- merge(img_matrix_1, matched_pts, by = 'Iavarone_ID', all.y = FALSE)
ids <- img_matrix_1$Iavarone_ID
img_matrix_1$Iavarone_ID <- NULL
img_matrix <- data.frame(t(img_matrix_1))
img_Symbol <- rownames(img_matrix)
img_matrix <- data.frame(lapply(img_matrix, as.numeric))
colnames(img_matrix) <- ids
row.names(img_matrix) <- NULL
no_of_imgs <- dim(img_matrix)[1]
img_names <- img_Symbol

# For each img we make the tree
for (img_no in seq(1,no_of_imgs)) {
  
  # Print the img name and the number
  print(paste("Working on img no.", img_no, ":", img_names[img_no]))
  
  # Get the img  data
  ge_df_matrix <- as.matrix(img_matrix[1:length(img_matrix)])
  rownames(ge_df_matrix) <- img_Symbol
  
  # Transpose the above matrix
  transposed_ge_df <- t(ge_df_matrix)
  transposed_ge_df <- as.data.frame(transposed_ge_df)
  
  # Get the data for the img
  img_df <- as.data.frame(transposed_ge_df[img_no])
  # head(img_df)
  img_name = colnames(img_df)
  
  # Get the img events
  img_events <- get_img_events(img_df = img_df)
  
  # Create an empty data frame to append the result of each img tree events
  df <- data.frame(matrix(ncol = 7, nrow = 0))
  
  for (i in seq(1, length(img_events))) {
    
    # Get event characteristics
    event <- i
    # Get event patients
    unlisted_patients <- unlist(img_events[i])
    patients <- paste(unlisted_patients, collapse=",")
    # Get median img 
    patients_ge <- img_df %>% filter(row.names(img_df) %in% unlisted_patients)
    median_ge <- median(patients_ge[[1]])
    # Get the number of patients and if there is only one patient, its a leaf
    no_of_patients <- length(unlisted_patients)
    if (no_of_patients == 1) {
      leaf_status <- 1
    } else {
      leaf_status <- 0
    }
    tissue_region <- tissue_region
    # Append the above data frame to df
    df <- rbind(df, c(event, patients, median_ge, no_of_patients, img_name, 'img_value', leaf_status, tissue_region))
  }
  
  # Provide column names
  colnames(df) <- c('event', 'samples', 'median_signal', 'no_of_patients', 'img_name', 'molecular_data', 'leaf_status', 'tissue_region')
  
  file_name = paste(img_name, ".csv", sep="")
  write.csv(df, file=file_name)
}

# fixing the files

# Set the directory where your CSV files are located
directory_path <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/NE"

# List all CSV files in the directory
csv_files <- list.files(directory_path, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file and make the required changes
for (csv_file in csv_files) {
  # Read the CSV file
  df <- read.csv(csv_file, header = TRUE)
  
  # Rename the "patients" column to "samples"
  colnames(df)[colnames(df) == "patients"] <- "samples"
  colnames(df)[colnames(df) == "median_ge"] <- "median_signal"
  
  
  # Remove the first column
  df <- df[, -1]
  
  # Write the modified dataframe back to the CSV file
  write.csv(df, file = csv_file, row.names = FALSE)
}

# Print a message when all files are processed
cat("All CSV files have been updated.\n")

# CE ----
#CE
setwd("/Users/m254284/Desktop/Neo4j_Projects/Img_events/CE")
tissue_region <- "CE"
mri_annotation <- data.frame(master$MRI.contrast.enhancing.annotation, master$Iavarone_ID)
colnames(mri_annotation) <- c("MRI.contrast.enhancing.annotation", "Iavarone_ID")
img_matrix_1 <- img
img_matrix_1 <- merge(img_matrix_1, mri_annotation, by = "Iavarone_ID")
img_matrix_1 <- subset(img_matrix_1, MRI.contrast.enhancing.annotation == tissue_region)
img_matrix_1$MRI.contrast.enhancing.annotation <- NULL
img_matrix_1 <- na.omit(img_matrix_1)

#prepare size of new subset and also ID collection dataframe
count_NE <- table(all_nodes_deg$MRI.contrast.enhancing.annotation.x)["NE"]
sim_sampleID <- data.frame(Index = 1:count_NE)

for (s in 1:count_NE) {set.seed(s)  # Replace 123 with your desired seed value
  #finish matrix prep for dendrogram
  matched_pts <- data.frame(all_nodes_deg$Iavarone_ID)
  colnames(matched_pts) <- "Iavarone_ID"
  img_matrix_1 <- merge(img_matrix_1, matched_pts, by = 'Iavarone_ID', all.y = FALSE)
  ids <- img_matrix_1$Iavarone_ID
  #need to bootstrap here
  sample_size <- count_NE  # Adjust as needed
  # Randomly sample rows based on the seed
  sampled_rows <- sample(nrow(img_matrix_1), size = sample_size, replace = FALSE)
  # Create a new dataframe with the sampled rows, save in simulation sample ID dataframe
  img_matrix_2 <- img_matrix_1[sampled_rows, ]
  col <- paste("Sim", s, sep = "_")
  ids <- img_matrix_2$Iavarone_ID
  sim_sampleID[[col]] <- ids
  sim_sampleID$Index <- NULL
  
  img_matrix_2$Iavarone_ID <- NULL
  img_matrix <- data.frame(t(img_matrix_2))
  img_Symbol <- rownames(img_matrix)
  img_matrix <- data.frame(lapply(img_matrix, as.numeric))
  colnames(img_matrix) <- ids
  row.names(img_matrix) <- NULL
  no_of_imgs <- dim(img_matrix)[1]
  img_names <- img_Symbol
  
  # For each img we make the tree
  for (img_no in seq(1,no_of_imgs)) {
    
    # Print the img name and the number
    print(paste("Working on img no.", img_no, ":", img_names[img_no]))
    
    # Get the img  data
    ge_df_matrix <- as.matrix(img_matrix[1:length(img_matrix)])
    rownames(ge_df_matrix) <- img_Symbol
    
    # Transpose the above matrix
    transposed_ge_df <- t(ge_df_matrix)
    transposed_ge_df <- as.data.frame(transposed_ge_df)
    
    # Get the data for the img
    img_df <- as.data.frame(transposed_ge_df[img_no])
    # head(img_df)
    img_name = colnames(img_df)
    
    # Get the img events
    img_events <- get_img_events(img_df = img_df)
    
    # Create an empty data frame to append the result of each img tree events
    df <- data.frame(matrix(ncol = 7, nrow = 0))
    
    for (i in seq(1, length(img_events))) {
      
      # Get event characteristics
      event <- i
      # Get event patients
      unlisted_patients <- unlist(img_events[i])
      patients <- paste(unlisted_patients, collapse=",")
      # Get median img 
      patients_ge <- img_df %>% filter(row.names(img_df) %in% unlisted_patients)
      median_ge <- median(patients_ge[[1]])
      # Get the number of patients and if there is only one patient, its a leaf
      no_of_patients <- length(unlisted_patients)
      if (no_of_patients == 1) {
        leaf_status <- 1
      } else {
        leaf_status <- 0
      }
      tissue_region <- tissue_region
      sim <- s
      # Append the above data frame to df
      df <- rbind(df, c(event, patients, median_ge, no_of_patients, img_name, 'img_value', leaf_status, tissue_region, sim))
    }
    
    # Provide column names
    colnames(df) <- c('event', 'samples', 'median_signal', 'no_of_patients', 'img_name', 'molecular_data', 'leaf_status', 'tissue_region', 'sim')
    
    file_name = paste(img_name, ".csv", sep="")
    base_directory_path <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/CE/"
    # simulation within directory
    variable_to_add <- paste("Sim", s, sep = "_")  # Example: "Value_1", "Value_2", ...
    # Combine the base path and the variable to create the new directory path
    directory_path <- file.path(base_directory_path, variable_to_add)
    setwd(directory_path)
    write.csv(df, file=file_name)
  }
  
  # fixing the files. dont actually need this
  
  # Set the directory where your CSV files are located
  base_directory_path <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/CE/"
  
  # simulation within directory
  variable_to_add <- paste("Sim", s, sep = "_")  # Example: "Value_1", "Value_2", ...
  
  # Combine the base path and the variable to create the new directory path
  directory_path <- file.path(base_directory_path, variable_to_add)
  
  # List all CSV files in the directory
  csv_files <- list.files(directory_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Loop through each CSV file and make the required changes
  for (csv_file in csv_files) {
    # Read the CSV file
    df <- read.csv(csv_file, header = TRUE)
    
    # Rename the "patients" column to "samples"
    colnames(df)[colnames(df) == "patients"] <- "samples"
    colnames(df)[colnames(df) == "median_ge"] <- "median_signal"
    
    
    # Remove the first column
    df <- df[, -1]
    
    # Write the modified dataframe back to the CSV file
    write.csv(df, file = csv_file, row.names = FALSE)
  }
  # Print a message when all files are processed
  cat("All CSV files have been updated.\n")
}

# gathering the bootstrapped sample IDs ----
#gathering the bootstrapped sample IDs
sim_sampleID <- data.frame(Index = 1:count_NE)
for (s in 1:count_NE) {set.seed(s)  # Replace 123 with your desired seed value
  #finish matrix prep for dendrogram
  matched_pts <- data.frame(all_nodes_deg$Iavarone_ID)
  colnames(matched_pts) <- "Iavarone_ID"
  img_matrix_1 <- merge(img_matrix_1, matched_pts, by = 'Iavarone_ID', all.y = FALSE)
  ids <- img_matrix_1$Iavarone_ID
  #need to bootstrap here
  sample_size <- count_NE  # Adjust as needed
  # Randomly sample rows based on the seed
  sampled_rows <- sample(nrow(img_matrix_1), size = sample_size, replace = FALSE)
  # Create a new dataframe with the sampled rows, save in simulation sample ID dataframe
  img_matrix_2 <- img_matrix_1[sampled_rows, ]
  col <- paste("Sim", s, sep = "_")
  ids <- img_matrix_2$Iavarone_ID
  sim_sampleID[[col]] <- ids 
  sim_sampleID$Index <- NULL}

setwd("/Users/m254284/Desktop/Neo4j_Projects/Img_events/CE")
write.csv(sim_sampleID, "sim_sampleID.csv", row.names = FALSE)




