# IMAGING COMMUNITIES BY CALCULATION TYPE
# INITIAL ITERATIONS PROCESSING ----
# collect sample IDs for further analysis on community/comunity basis
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)

# Set the main directory containing subdirectories with CSV files
main_directory <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Iterations"

# Function to read CSV files and extract relevant data
read_csv_and_extract <- function(file_path) {
  df <- read.csv(file_path)
  df %>%
    group_by(community) %>%
    summarise(
      samples = list(unique(unlist(strsplit(samples, ", ")))),
      Imgs = list(unique(img_name)),
      quartiles = list(quartile)
    ) %>%
    mutate(community = paste0(community, "_", tools::file_path_sans_ext(basename(file_path)))) %>%
    mutate(samples = gsub("^c\\(|\\)$", "", samples)) %>%
    mutate(Imgs = gsub("^c\\(|\\)$", "", Imgs)) %>%
    mutate(quartiles = sapply(quartiles, function(q) paste(q, collapse = ",")))
}

# List all subdirectories in the main directory
subdirs <- list.dirs(main_directory, full.names = TRUE, recursive = FALSE)

# Function to remove duplicates from a string
remove_duplicates <- function(x) {
  unique_values <- unique(unlist(strsplit(x, ",")))
  paste(unique_values, collapse = ",")
}

# Process each subdirectory separately
for (subdir in subdirs) {
  # List all CSV files in the current subdirectory
  csv_files <- list.files(subdir, pattern = "\\.csv$", full.names = TRUE)
  # Read and process all CSV files in the current subdirectory
  data_list <- purrr::map_dfr(csv_files, read_csv_and_extract)
  # Clean data
  data_list <- data.frame(lapply(data_list, function(x) gsub("\"", "", x)))
  # Remove spaces from the "samples" column and remove duplicates
  data_list$samples <- sapply(data_list$samples, function(x) {
    remove_duplicates(gsub(" ", "", x))
  })
  # Define the output file name based on the subdirectory name
  subdir_name <- basename(subdir)
  output_file <- paste0("samples_by_comm_", subdir_name, ".csv")
  # Set the output directory
  setwd('/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration')
  # Save the combined dataframe to a new CSV file
  write.csv(data_list, output_file, row.names = FALSE)
}
# NE AND CE PERCENTAGE BY COMMUNITY ----
# Load required libraries
library(dplyr)
library(purrr)
library(tools)
setwd('/Users/m254284/Desktop/Multiregional_/for_matt')
master <- read.delim('master.tsv', header = TRUE, sep = "\t")

main_directory <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Iterations"
output_directory <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration"

# Function to process each output CSV file separately
process_output_csv <- function(file_path) {
  # Load master dataframe
  master_filtered <- subset(master, MRI.contrast.enhancing.annotation %in% c("NE", "CE"))
  
  # Calculate the global percentage of "NE" and "CE" samples
  total_master_samples <- nrow(master_filtered)
  NE_master <- sum(master_filtered$MRI.contrast.enhancing.annotation == "NE")
  CE_master <- sum(master_filtered$MRI.contrast.enhancing.annotation == "CE")
  global_percentage_NE <- NE_master / total_master_samples
  global_percentage_CE <- CE_master / total_master_samples
  
  # Read the current output CSV file
  df <- read.csv(file_path)
  
  # Create an empty dataframe to store results
  output <- data.frame(row_number = integer(), num_samples = integer(), percentage_NE = numeric(), percentage_CE = numeric(), chi_squared_statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each row of df
  for (i in 1:nrow(df)) {
    # Split sample IDs in the current row
    sample_ids <- unlist(strsplit(df$samples[i], ","))
    
    # Number of samples
    num_samples <- length(sample_ids)
    
    # Subset master dataframe based on sample IDs
    subsetted_master <- master_filtered[master_filtered$Iavarone_ID %in% sample_ids, ]
    
    # Calculate percentage of samples with "NE" or "CE" in the "MRI.contrast.enhancing.annotation" column
    total_samples <- nrow(subsetted_master)
    NE_samples <- sum(subsetted_master$MRI.contrast.enhancing.annotation == "NE")
    CE_samples <- sum(subsetted_master$MRI.contrast.enhancing.annotation == "CE")
    
    # Proceed only if there are samples with "NE" or "CE"
    if (total_samples > 0) {
      percentage_NE <- NE_samples / total_samples
      percentage_CE <- CE_samples / total_samples
      
      # Perform Chi-squared test comparing with global percentages
      expected_NE <- global_percentage_NE * num_samples
      expected_CE <- global_percentage_CE * num_samples
      
      # Normalize expected values to sum up to 1
      total_expected <- expected_NE + expected_CE
      expected <- c(expected_NE / total_expected, expected_CE / total_expected)
      
      observed <- c(NE_samples, CE_samples)
      
      chi_squared_test <- chisq.test(observed, p = expected)
      
      # Store results in output dataframe
      output[i, "row_number"] <- i
      output[i, "num_samples"] <- num_samples
      output[i, "percentage_NE"] <- percentage_NE
      output[i, "percentage_CE"] <- percentage_CE
      output[i, "chi_squared_statistic"] <- chi_squared_test$statistic
      output[i, "p_value"] <- chi_squared_test$p.value
      
      # Append columns from df
      output[i, names(df)] <- df[i, ]
    } else {
      # If no samples with "NE" or "CE", store NA for percentages and p-value
      output[i, c("row_number", "num_samples", "percentage_NE", "percentage_CE", "chi_squared_statistic", "p_value")] <- NA
    }
  }
  
  # Filter out rows with p-value > 0.05
  output <- output[output$p_value <= 0.05,]
  
  # Save the results to a new CSV file
  output_file_name <- paste0("NECE_comm_percentages_", tools::file_path_sans_ext(basename(file_path)), ".csv")
  write.csv(output, file.path(output_directory, output_file_name), row.names = FALSE)
}

# Process each output CSV file generated in the previous step
output_files <- list.files(output_directory, pattern = "^samples_by_comm_.*\\.csv$", full.names = TRUE)

for (output_file in output_files) {
  process_output_csv(output_file)
}

# FILTER HIGH NE COMMUNITIES ----
setwd("/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration")

# Function to calculate expected percentage of "NE" and subset the dataframe
process_comms_file <- function(file_path, master) {
  # Load the 'comms' dataframe
  comms <- read.csv(file_path)
  
  # Calculate the expected percentage of "NE" occurrence in the 'master' dataframe
  cene <- master[master$MRI.contrast.enhancing.annotation == "NE" | master$MRI.contrast.enhancing.annotation == "CE",]
  cene <- cene$MRI.contrast.enhancing.annotation
  cene <- na.omit(cene)
  expected_percentage_NE <- sum(cene == "NE") / length(cene)
  
  # Subset the 'comms' dataframe based on the condition
  comms_NE <- comms[comms$percentage_NE > expected_percentage_NE, ]
  
  # Define the output file name
  output_file_name <- paste0("filtered_", basename(file_path))
  
  # Save the subsetted dataframe to a new CSV file
  write.csv(comms_NE, output_file_name, row.names = FALSE)
}

# List all the NECE comm percentage files
nece_files <- list.files("/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration", pattern = "^NECE_comm_percentages_.*\\.csv$", full.names = TRUE)

# Assuming 'master' dataframe is already loaded in the environment
for (file in nece_files) {
  process_comms_file(file, master)
}

# PATHWAY ENRICHMENT IN VS OUT ----
setwd("/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration")
load("combat_counts.RData")
cc <- data.frame(t(adjusted))
cc$Iavarone_ID <- rownames(cc)
merged_df <- merge(cc, master[, c("Iavarone_ID", "Patient_ID.x", "MRI.contrast.enhancing.annotation")], by = "Iavarone_ID", all.x = TRUE)

# create an in vs. out column
library(dplyr)
library(purrr)

is_in <- function(samples, iavarone_ids, data) {
  sapply(iavarone_ids, function(id) {
    if (id %in% unlist(samples)) {
      "in"
    } else if ("NE" %in% data[data$Iavarone_ID == id, "MRI.contrast.enhancing.annotation"]) {
      "out"
    } else {
      
    }
  })
}

# Load required libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set directories
integration_dir <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration"
pathways_dir <- file.path(integration_dir, "Pathways_regional/NE_All")
method_list_dir <- file.path(integration_dir, "Method1_Lists")

# Ensure the output directories exist
if (!dir.exists(pathways_dir)) {
  dir.create(pathways_dir, recursive = TRUE)
}
if (!dir.exists(method_list_dir)) {
  dir.create(method_list_dir, recursive = TRUE)
}

# List all the filtered CSV files
filtered_files <- list.files(integration_dir, pattern = "^filtered_.*\\.csv$", full.names = TRUE)

# Function to process each filtered CSV file
process_filtered_file <- function(file_path) {
  df_NE <- read.csv(file_path)
  base_file_name <- tools::file_path_sans_ext(basename(file_path))
  
  for (i in 1:nrow(df_NE)) {
    current_row <- df_NE[i, ]
    new_col_name <- paste("in_out_", current_row$community, sep = "")
    merged_df[[new_col_name]] <- factor(
      is_in(
        strsplit(current_row$samples, ","),
        unique(merged_df$Iavarone_ID),
        merged_df
      ),
      levels = c("in", "out")
    )
  }
  
  countData <- as.matrix(t(cc[cc$Iavarone_ID %in% merged_df$Iavarone_ID,]))
  row_to_remove <- which(rownames(countData) == "Iavarone_ID")
  countData <- countData[-row_to_remove, ]
  mode(countData) <- "integer"  # Convert to integer
  in_out_cols <- grep("in_out", names(merged_df), value = TRUE)
  colData <- merged_df[, c("Iavarone_ID", "Patient_ID.x", in_out_cols)]
  
  # Start loop here
  for (i in 1:length(in_out_cols)) {
    feature <- in_out_cols[i]
    
    message(paste0("Start ", feature))
    comparison <- paste0(feature)
    formula <- reformulate(c(feature, "Patient_ID.x"))
    colData_sub <- colData[, c("Iavarone_ID", "Patient_ID.x", feature)]
    colData_sub <- na.omit(colData_sub)
    countData_sub <- countData[, colnames(countData) %in% colData_sub$Iavarone_ID]
    
    # Check for DESeq issue with patient counts
    in_out_col_name <- grep("^in_out", names(colData_sub), value = TRUE)
    unique_patients_in <- unique(colData_sub[colData_sub[[in_out_col_name]] == "in", ]$Patient_ID.x)
    unique_patients_out <- unique(colData_sub[colData_sub[[in_out_col_name]] == "out", ]$Patient_ID.x)
    
    # Count the number of samples per patient for 'in' group
    samples_per_patient_in <- sapply(unique_patients_in, function(patient) {
      sum(colData_sub$Patient_ID.x == patient & colData_sub[[in_out_col_name]] == "in")
    })
    
    # Count the number of samples per patient for 'out' group
    samples_per_patient_out <- sapply(unique_patients_out, function(patient) {
      sum(colData_sub$Patient_ID.x == patient & colData_sub[[in_out_col_name]] == "out")
    })
    
    # Check conditions for both 'in' and 'out' groups
    if (length(unique_patients_in) < 2 || length(unique_patients_out) < 2) {
      message(paste0("Skipping covariate model for ", feature, " due to insufficient unique patients"))
      formula <- reformulate(feature)  # Use only the feature as the formula
    }
    
    # Check if any patient has only one sample in 'in' or 'out' group
    single_sample_in <- any(samples_per_patient_in == 1)
    single_sample_out <- any(samples_per_patient_out == 1)
    
    if (single_sample_in || single_sample_out) {
      message(paste0("Applying different formula for ", feature, " due to single sample per patient"))
      formula <- reformulate(feature)  # Use only the feature as the formula
    }
    
    # New condition: Check if there is only one patient and they have only one sample
    if ((length(unique_patients_in) == 1 && samples_per_patient_in == 1) ||
        (length(unique_patients_out) == 1 && samples_per_patient_out == 1)) {
      message(paste0("Skipping ", feature, " due to single patient with only one sample in either 'in' or 'out' group"))
      next
    }
    
    # New condition: Check if all samples are the same in the "in" or "out" group
    if (all(colData_sub[[in_out_col_name]] == "in") || all(colData_sub[[in_out_col_name]] == "out")) {
      message(paste0("Skipping ", feature, " due to all samples being in the same group"))
      next
    }
    
    dds <- DESeqDataSetFromMatrix(countData_sub, colData_sub, design = formula) # MAKE DEG LIST
    dds <- DESeq(dds)
    
    results <- nbinomWaldTest(dds, maxit = 1000)
    res <- results(dds, contrast = c(feature, "in", "out"))
    
    # Filter results for significant genes (padj < 0.05)
    sig_res <- subset(res, padj < 0.05)
    
    # all genes through pathway enrichment
    gene_names <- rownames(res)
    stat_values <- res$stat
    method_1 <- data.frame(Gene = gene_names, Statistic = stat_values)
    method_1_list <- list()
    for (gene in 1:nrow(method_1)) {
      gene_name <- method_1$Gene[gene]
      method_1_list[[gene_name]] <- method_1$Statistic[gene]
    }
    
    method_1_atomic <- unlist(method_1_list)
    method_1_atomic <- rev(sort(method_1_atomic))
    
    GSA_GO <- gseGO(geneList = method_1_atomic,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "ALL",
                    nPerm = 1000,
                    pvalueCutoff = 1,
                    verbose = FALSE)
    
    GSEA_GO_df <- GSA_GO@result
    #GSEA_GO_df_sig <- GSEA_GO_df[GSEA_GO_df$p.adjust <= 0.05,]
    GSEA_GO_df_sig <- GSEA_GO_df[GSEA_GO_df$pvalue <= 0.05,]
    output_file_name <- file.path(pathways_dir, paste0(base_file_name, "_", comparison, "_NE_GSEA.csv"))
    write.csv(GSEA_GO_df_sig, output_file_name, row.names = FALSE)
    
    # Save the significant genes with padj < 0.05 to a CSV file
    sig_genes_file_name <- file.path(method_list_dir, paste0(base_file_name, "_", comparison, "_significant_genes.csv"))
    write.csv(as.data.frame(sig_res), sig_genes_file_name, row.names = TRUE)
  }
}

# Process each filtered CSV file
for (file in filtered_files) {
  process_filtered_file(file)
}


# PLOTTING BIOLOGICAL PROCESSES ENRICHMENT BY COMMUNITY ----
# Load required libraries
library(dplyr)
library(stringr)
library(ggplot2)

# Set directories
pathways_dir <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration/Pathways_regional/NE_All"
output_dir <- "/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration/Pathways_regional"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# List all the files in the pathways directory
ffile <- list.files(pathways_dir, pattern = "\\.csv$", full.names = TRUE)

# Define a function to create a bubble plot for a given group of files
create_bubble_plot <- function(files, group_name) {
  df <- NULL
  for (file in files) {
    gsea <- read.csv(file)
    
    if (nrow(gsea) > 0) {
      #gsea <- gsea[gsea$p.adjust < 0.05,]
      gsea <- gsea[gsea$pvalue < 0.05,]
      
      # Extract the sample name
      sample_name <- sub(".*p([0-9]{2}).*", "p\\1", basename(file))
      
      tmp <- data.frame(
        Sample = sample_name,
        Signature = gsea$Description,
        NES = gsea$NES,
        FDR = -log10(gsea$p.adjust),
        Class = ".",
        stringsAsFactors = FALSE
      )
      
      tmp <- tmp[order(tmp$NES, decreasing = TRUE), ]
      tmp[1:50, "Class"] <- "select"
      tmp[nrow(tmp):nrow(tmp)-49, "Class"] <- "select"
      df <- rbind(df, tmp)
    } else {
      cat("No rows found in:", file, "\n")
    }
  }
  
  tmp <- unique(df[df$Class == "select", "Signature"])
  df <- df[df$Signature %in% tmp, ]
  
  range(df$NES) ## change in limits scale_fill_gradient2
  range(df$FDR) ## change range of scale_size_continuous
  
  table(df$Sample)
  df$Sample <- gsub("_GSEA.csv", "", df$Sample)
  df$Sample <- gsub("in_out_", "", df$Sample)
  table(df$Sample)
  unique_samples <- unique(df$Sample)
  
  df$Sample <- factor(df$Sample, levels = unique_samples)
  
  # Order the dataframe based on NES
  tmp <- df[order(df$NES, decreasing = F), "Signature"]
  tmp <- unique(tmp)
  df$Signature <- factor(df$Signature, levels = tmp)
  
  # Create Type column
  df$Type <- "Up"
  df[df$NES < 0, "Type"] <- "Down" 
  df$Type <- factor(df$Type)
  
  tmp <- names(sort(table(df$Signature)))
  df$Signature <- factor(df$Signature, levels = tmp)
  
  #df_up <- df[df$Type == "Up",]
  df_up <- df
  df_up$Sample <- factor(df_up$Sample, levels = unique(df_up$Sample[order(df_up$Signature)]))
  df_up <- df_up[df_up$Class == "select",]
  
  # Create a similarity matrix based on shared signatures
  similarity_matrix <- table(df_up$Sample, df_up$Signature)
  similarity_matrix <- similarity_matrix %*% t(similarity_matrix)
  diag(similarity_matrix) <- 0 # Remove self-similarity
  print(similarity_matrix)
  
  # Use hierarchical clustering to order samples
  hc <- hclust(dist(as.matrix(similarity_matrix)))
  
  # Determine number of clusters (you can adjust the number of clusters as needed)
  num_clusters <- length(unique(df_up$Sample))
  
  # Cut the dendrogram to form clusters
  clusters <- cutree(hc, k = num_clusters)
  df_up$Cluster <- as.factor(clusters[match(df_up$Sample, names(clusters))])
  
  # Order samples based on clustering
  ordered_samples <- rownames(similarity_matrix)[hc$order]
  df_up$Sample <- factor(df_up$Sample, levels = ordered_samples)
  
  # Output PDF file name
  output_pdf <- file.path(output_dir, paste("bubble_", "all_", group_name, ".pdf", sep = ""))
  
  # Create the bubble plot and save it as a PDF
  plot <- ggplot(df_up, aes(x = Sample, y = Signature)) + 
    geom_point(aes(size = FDR, fill = NES), alpha = 0.9, shape = 21) + 
    scale_size_continuous(limits = c(0, 10), range = c(0,10), breaks = c(0, 2, 10)) + #### FDR
    labs(x = "Community", y = "Pathway", size = "-log10(p.adjusted)", fill = "NES") + 
    theme_bw() +
    theme(
      legend.key=element_blank(), 
      #axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle = 45, vjust = 0.8, hjust = 0.9), 
      axis.text.x = element_blank(),
      #axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
      axis.text.y = element_blank(),
      legend.text = element_text(size = 10, face ="bold", colour ="black"), 
      legend.title = element_text(size = 12), 
      legend.position = "right"
    ) +  
    scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(-5, 5), midpoint = 0, breaks = c(-5,0, 5)) + #### NES
    facet_grid( ~ Cluster, scales = "free") +
    theme(strip.text.y = element_text(size = 6, face = "bold"))
  
  print(plot)
  
  #ggsave(output_pdf, plot, width = 15, height = 15)
  
  # Iterate over each cluster to create a separate plot
  for (cluster in levels(df_up$Cluster)) {
    # Filter the dataframe for the current cluster
    cluster_df <- df_up[df_up$Cluster == cluster, ]
    
    # Output PDF file name for the current cluster
    output_pdf <- file.path(output_dir, paste("bubble_", "all", "_Cluster_", cluster, group_name, ".pdf", sep = ""))
    
    # Create the bubble plot for the current cluster and save it as a PDF
    plot <- ggplot(cluster_df, aes(x = Sample, y = Signature)) + 
      geom_point(aes(size = FDR, fill = NES), alpha = 0.9, shape = 21) + 
      scale_size_continuous(limits = c(0, 10), range = c(0,10), breaks = c(0, 2, 10)) + # FDR
      labs(x = "Community", y = "", size = "-log10(p.adjusted)", fill = "NES") + 
      theme_bw() +
      theme(
        legend.key=element_blank(), 
        #axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle = 45, vjust = 0.8, hjust = 0.9), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12), 
        legend.position = "right"
      ) +  
      scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(-5, 5), midpoint = 0, breaks = c(-5, 0, 5)) + # NES
      theme(strip.text.y = element_text(size = 6, face = "bold"))
    
    #ggsave(output_pdf, plot, width = 5, height = 10)
  }
  return(df_up)
}

# Group files by BSW, UCLA, and QWS
bsw_files <- ffile[grep("BSW", ffile)]
ucla_files <- ffile[grep("UCLA", ffile)]
qws_files <- ffile[grep("QWS", ffile)]
nc_files <- ffile[grep("nc", ffile)]

#clusts_num <- 2

# Create bubble plots for each group
bsw_up <- create_bubble_plot(bsw_files, "BSW")
ucla_up <- create_bubble_plot(ucla_files, "UCLA")
qws_up <- create_bubble_plot(qws_files, "QWS")
nc_up <- create_bubble_plot(nc_files, "nc")

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# Assume df_up_NE is your input dataframe
# Step 1: Aggregate the data with count

# what calc type
df_up_all <- bsw_up

# Group by Cluster and Signature to calculate the Count and Average NES
agg_data <- df_up_all %>%
  group_by(Cluster, Signature) %>%
  summarise(
    Count = n(),
    Avg_NES = mean(NES, na.rm = TRUE),  # Calculate average NES
    .groups = 'drop'
  )

# Calculate the total count per cluster for normalization
cluster_sizes <- df_up_all %>%
  group_by(Cluster) %>%
  summarise(
    TotalCount = n(),
    .groups = 'drop'
  )

# Add Frequency column (normalized count) and include Avg_NES in the final output
agg_data <- agg_data %>%
  left_join(cluster_sizes, by = "Cluster") %>%
  mutate(Frequency = Count / TotalCount) %>%
  dplyr::select(Cluster, Signature, Count, Frequency, Avg_NES)

# Step 3: Filter top 10 signatures per cluster based on count
top_agg_data <- agg_data %>%
  group_by(Cluster) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  slice_head(n = 40) %>%
  ungroup()


# Step 2: Create a contingency table
heatmap_data <- tidyr::spread(top_agg_data, Signature, as.numeric(Avg_NES), fill = 0)

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data[, -c(1,2)])
rownames(heatmap_matrix) <- heatmap_data$Cluster

# Step 3a: Generate the heatmap using ggplot2
ggplot_figure <- ggplot(top_agg_data, aes(x = Cluster, y = Signature, fill = Avg_NES)) +
  geom_tile() +
  scale_fill_gradient(low = "purple", high = "yellow") +
  theme_classic() +
  labs(title = "Top 20 Pathways per Cluster",
       x = "Cluster",
       y = "Signature",
       fill = "Mean NES per Cluster")

# Print ggplot2 heatmap
print(ggplot_figure)

# PLOTTING IMAGING FEATURES BY COMMUNITY ----
# Load required libraries
library(dplyr)
library(tidyr)

cluster_samples <- df_up_all %>%
  group_by(Cluster) %>%
  summarise(Samples = list(unique(Sample)), .groups = 'drop')

cluster_sample_list <- setNames(cluster_samples$Samples, cluster_samples$Cluster)

clean_cluster_sample_list <- lapply(cluster_sample_list, function(samples) {
  gsub("filtered_NECE_comm_percentages_samples_by_comm_BSW_", "", samples)
})

clean_cluster_sample_list <- lapply(clean_cluster_sample_list, function(samples) {
  gsub("_NE", "", samples)
})

setwd('/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Integration')
# Read the samples_by_comm.csv dataframe
samples_by_comm <- read.csv("samples_by_comm_BSW.csv")

# Create an empty dataframe to store the filtered results
filtered_samples_by_comm <- data.frame()

# Iterate over each cluster in the cluster_sample_list
for (cluster in names(clean_cluster_sample_list)) {
  # Extract samples for the current cluster
  samples_in_cluster <- list(clean_cluster_sample_list[[cluster]])
  unlist(samples_in_cluster)
  # Filter rows in samples_by_comm where the community column matches any sample in the current cluster
  matching_rows <- samples_by_comm %>%
    filter(community %in% unlist(samples_in_cluster))
  
  # Add the Cluster column to indicate which cluster the sample belongs to
  matching_rows <- matching_rows %>%
    mutate(Cluster = cluster)
  
  # Append the filtered rows to the result dataframe
  filtered_samples_by_comm <- bind_rows(filtered_samples_by_comm, matching_rows)
}

# Print the filtered results
print(filtered_samples_by_comm)

filtered_samples_by_comm$Imgs <- gsub(" ", "", filtered_samples_by_comm$Imgs)

# Count the number of instances of each value in the imgs column for each cluster
img_counts <- filtered_samples_by_comm %>%
  # Split the imgs column into separate rows by comma
  separate_rows(Imgs, sep = ",") %>%
  # Count the occurrences of each img value for each cluster
  group_by(Cluster, Imgs) %>%
  summarise(Count = n(), .groups = 'drop')

# Compute total number of imgs entries per cluster for normalization
total_counts <- img_counts %>%
  group_by(Cluster) %>%
  summarise(Total = sum(Count), .groups = 'drop')

# Merge total counts with img_counts to compute frequencies
img_frequencies <- img_counts %>%
  left_join(total_counts, by = "Cluster") %>%
  mutate(Frequency = Count / Total) %>%
  dplyr::select(Cluster, Imgs, Frequency)

# Pivot the data to a wide format suitable for heatmap
heatmap_data <- img_frequencies %>%
  pivot_wider(names_from = Cluster, values_from = Frequency, values_fill = list(Frequency = 0))

# Print the heatmap data
print(heatmap_data)

# Plot the heatmap using ggplot2
ggplot(img_frequencies, aes(x = Cluster, y = Imgs, fill = Frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "pink", high = "blue") +
  labs(title = "Feature Frequencies",
       x = "Cluster",
       y = "Imaging Feature",
       fill = "Frequency") +
  theme_classic()


# PLOTTING IMAGING SIGNAL FOR CLUSTERS ----
library(dplyr)
library(tidyr)

samples_combined <- filtered_samples_by_comm %>%
  # Group by Cluster
  group_by(Cluster) %>%
  # Extract and remove duplicates from the Samples column
  summarise(Samples = paste(unique(samples), collapse = ","), .groups = 'drop') %>%
  # Ensure unique values in the Samples column
  mutate(Samples = sapply(Samples, function(x) {
    # Split the samples, remove duplicates, and then recombine
    unique_samples <- unique(unlist(strsplit(x, ",")))
    paste(unique_samples, collapse = ",")
  }))

#get imaging intensities
setwd("/Users/m254284/Desktop/Multiregional_/for_matt")
master <- read.delim('master.tsv', header = TRUE, sep = "\t")
library(dplyr)
library(tidyr)

# Convert Samples column into a list of vectors for easy matching
samples_list <- split(samples_combined$Samples, samples_combined$Cluster)

# Get the list of relevant imaging features from img_frequencies
relevant_imgs <- unique(img_frequencies$Imgs)

# Debugging print statements
print("Samples List:")
print(samples_list)
print("Relevant Imaging Features:")
print(relevant_imgs)

# First, select and scale the relevant columns across the entire dataset
scaled_master_all <- master %>%
  dplyr::select(all_of(c(relevant_imgs, "Iavarone_ID", "MRI.contrast.enhancing.annotation"))) %>%
  mutate(across(all_of(relevant_imgs), scale)) # Scaling across the entire dataset

scaled_master_NE <- master[master$MRI.contrast.enhancing.annotation == "NE",]

scaled_master_NE <- scaled_master_NE %>%
  dplyr::select(all_of(c(relevant_imgs, "Iavarone_ID", "MRI.contrast.enhancing.annotation"))) %>%
  mutate(across(all_of(relevant_imgs), scale)) # Scaling across the entire dataset

#select which to use
scaled_master <- scaled_master_NE

# Now, calculate the mean intensity values for each cluster
mean_intensity_per_cluster <- lapply(names(samples_list), function(cluster) {
  # Get sample IDs for the current cluster
  sample_ids <- strsplit(samples_list[[cluster]], ",")[[1]]
  sample_ids <- gsub(" ", "", sample_ids)
  
  # Filter the scaled master dataframe based on the sample IDs
  filtered_master <- scaled_master %>%
    filter(Iavarone_ID %in% sample_ids)
  
  # Remove the unnecessary columns
  filtered_master$Iavarone_ID <- NULL
  filtered_master$MRI.contrast.enhancing.annotation <- NULL
  
  # Calculate mean intensity values for each relevant imaging feature
  mean_intensities <- filtered_master %>%
    dplyr::select(all_of(relevant_imgs)) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "Imgs", values_to = "MeanIntensity") %>%
    mutate(Cluster = cluster)
  
  return(mean_intensities)
}) %>%
  bind_rows()

# Create heatmap
heatmap_plot <- ggplot(mean_intensity_per_cluster, aes(x = Cluster, y = Imgs, fill = MeanIntensity)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "magenta") +
  labs(title = "Intensity vs NE Samples",
       x = "Cluster",
       y = "Imaging Features",
       fill = "Normalized Intensity") +
  theme_minimal()

# Print the heatmap
print(heatmap_plot)

# CORRELATION OF IMAGING FEATURES WITH PATHWAYS ----
library(dplyr)
library(tidyr)

samples_combined <- filtered_samples_by_comm %>%
  # Group by Cluster
  group_by(Cluster) %>%
  # Extract and remove duplicates from the Samples column
  summarise(Samples = paste(unique(samples), collapse = ","), .groups = 'drop') %>%
  # Ensure unique values in the Samples column
  mutate(Samples = sapply(Samples, function(x) {
    # Split the samples, remove duplicates, and then recombine
    unique_samples <- unique(unlist(strsplit(x, ",")))
    paste(unique_samples, collapse = ",")
  }))

#get imaging intensities
setwd("/Users/m254284/Desktop/Multiregional_/for_matt")
master <- read.delim('master.tsv', header = TRUE, sep = "\t")
library(dplyr)
library(tidyr)

# Convert Samples column into a list of vectors for easy matching
samples_list <- split(samples_combined$Samples, samples_combined$Cluster)

# Get the list of relevant imaging features from img_frequencies
relevant_imgs <- master[, grep("Raw_Mean", colnames(master))]
relevant_imgs <- relevant_imgs[, !grepl("nii.gz", colnames(relevant_imgs))]
relevant_imgs <- relevant_imgs[, !grepl("UCLA", colnames(relevant_imgs))]
relevant_imgs <- relevant_imgs[, !grepl("QWS", colnames(relevant_imgs))]
relevant_imgs <- relevant_imgs[, !grepl("nc", colnames(relevant_imgs))]
relevant_imgs <- relevant_imgs[, !grepl("1_.DTI", colnames(relevant_imgs))]
relevant_imgs <- relevant_imgs[, !grepl("2_.DTI", colnames(relevant_imgs))]
relevant_imgs <- relevant_imgs[, !grepl("nc", colnames(relevant_imgs))]
# relevant_imgs <- relevant_imgs[, !grepl("BSW", colnames(relevant_imgs))]
relevant_imgs <- colnames(relevant_imgs)

# Calculate mean intensity values for each cluster with scaling
mean_intensity_per_cluster <- lapply(names(samples_list), function(cluster) {
  # Get sample IDs for the current cluster
  sample_ids <- strsplit(samples_list[[cluster]], ",")[[1]]
  sample_ids <- gsub(" ", "", sample_ids)
  
  selected_scaled_master <- master %>%
    dplyr::select(all_of(c(relevant_imgs, "Iavarone_ID")))
  
  # Filter master dataframe based on the sample IDs
  # filtered_master <- master %>%
  #   filter(Iavarone_ID %in% sample_ids) %>%
  #   dplyr::select(all_of(relevant_imgs))
  
  filtered_master <- selected_scaled_master %>%    
    filter(Iavarone_ID %in% sample_ids)
  
  filtered_master$Iavarone_ID <- NULL
  
  filtered_master_scale <- data.frame(scale(filtered_master))
  #log_master_filter <- data.frame(log10(filtered_master))
  
  # Calculate scaled mean intensity values for each relevant imaging feature
  mean_intensities <- filtered_master_scale %>%
    dplyr::select(all_of(relevant_imgs)) %>%  # Select only the relevant imaging features
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "Imgs", values_to = "MeanIntensity") %>%
    mutate(Cluster = cluster)
  
  return(mean_intensities)
}) %>%
  bind_rows()

# Reorient the data to have each unique Imgs value as a new column
reoriented_df <- mean_intensity_per_cluster %>%
  pivot_wider(
    names_from = Imgs,  # This will create a new column for each unique value in Imgs
    values_from = MeanIntensity  # The values for each Imgs will come from MeanIntensity
  ) %>%
  rename(Community = Cluster)  # Rename the Cluster column to Community

# View the reoriented dataframe
head(reoriented_df)

# now combine with agg_data 
length(table(agg_data$Signature))

# Loop through each unique value in agg_data$Signature
for (signature in unique(agg_data$Signature)) {
  # Add a new column to reoriented_df named after the current signature
  reoriented_df[[signature]] <- NA
  
  # Populate the new column with Avg_NES values where Cluster matches Community
  reoriented_df <- reoriented_df %>%
    rowwise() %>%
    mutate(
      !!signature := ifelse(
        Community %in% agg_data$Cluster[agg_data$Signature == signature],
        agg_data$Avg_NES[agg_data$Signature == signature & agg_data$Cluster == Community],
        NA  # Keep NA if no matching Cluster is found
      )
    )
}

# View the updated dataframe
head(reoriented_df)

# now correlate imaging features with pathways
# Load necessary library
library(dplyr)

# Identify the imaging feature columns and signature columns
imgs_cols <- grep("Raw_Mean", colnames(reoriented_df), value = TRUE)
signature_cols <- setdiff(colnames(reoriented_df), c("Community", imgs_cols))

# Initialize an empty list to store results
correlation_results <- list()

# Define a threshold for the minimum number of valid observations
min_valid_obs <- 0.1*nrow(reoriented_df) # You can adjust this based on your data

# Loop through each combination of imaging feature and signature column
for (img in imgs_cols) {
  for (sig in signature_cols) {
    # Get the vectors of the two columns
    img_values <- reoriented_df[[img]]
    sig_values <- reoriented_df[[sig]]
    
    # Count the number of finite observations
    valid_obs <- sum(is.finite(img_values) & is.finite(sig_values))
    
    # Perform the correlation test only if there are enough valid observations
    if (valid_obs >= min_valid_obs) {
      cor_test <- cor.test(img_values, sig_values, method = "spearman", use = "complete.obs")
      
      # Store the results in a data frame
      correlation_results[[paste(img, sig, sep = "_")]] <- data.frame(
        Imgs = img,
        Signature = sig,
        Correlation_Coefficient = cor_test$estimate,
        P_Value = cor_test$p.value
      )
    } else {
      # If not enough valid observations, store NA values
      correlation_results[[paste(img, sig, sep = "_")]] <- data.frame(
        Imgs = img,
        Signature = sig,
        Correlation_Coefficient = NA,
        P_Value = NA
      )
    }
  }
}

# Combine the list into a single dataframe
correlation_results_df <- bind_rows(correlation_results)
correlation_results_df <- na.omit(correlation_results_df)
#correlation_results_df <- correlation_results_df[correlation_results_df$P_Value < 0.01, ]


# View the organized results
print(correlation_results_df)

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(reshape2)

# Prepare the data for heatmap with p-values
heatmap_data <- correlation_results_df %>%
  dplyr::select(Imgs, Signature, Correlation_Coefficient, P_Value)

# Melt the data for ggplot2
heatmap_melted <- melt(heatmap_data, id.vars = c("Imgs", "Signature"))

# Separate Correlation and P-Value data
correlation_data <- heatmap_melted %>%
  filter(variable == "Correlation_Coefficient") %>%
  rename(Correlation = value)

p_value_data <- heatmap_melted %>%
  filter(variable == "P_Value") %>%
  rename(P_Value = value)

# Merge correlation data with p-value data
heatmap_final <- merge(correlation_data, p_value_data, by = c("Imgs", "Signature"))

# Create a new column to control color intensity based on p-value
heatmap_final$alpha_value <- ifelse(heatmap_final$P_Value < 0.05, 1, 0.3)

# Create the heatmap with conditional coloring
ggplot(heatmap_final, aes(x = Imgs, y = Signature, fill = Correlation, alpha = alpha_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name = "Correlation") +
  scale_alpha(range = c(0.3, 1)) + # Adjust the transparency based on p-value
  theme_classic() +
  labs(x = "Imaging Features", y = "Signature", title = "nc Heatmap with p < 0.05 Threshold") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels for readability

# sig only

# Set the p-value threshold
p_value_threshold <- 0.01

# Filter Signatures that have at least one p-value below the threshold
significant_signatures <- correlation_results_df %>%
  group_by(Signature) %>%
  filter(any(P_Value < p_value_threshold))

# Prepare the data for heatmap with only significant Signatures
significant_heatmap_data <- significant_signatures %>%
  dplyr::select(Imgs, Signature, Correlation_Coefficient, P_Value)

# Melt the data for ggplot2
significant_heatmap_melted <- melt(significant_heatmap_data, id.vars = c("Imgs", "Signature"))

# Separate Correlation and P-Value data
significant_correlation_data <- significant_heatmap_melted %>%
  filter(variable == "Correlation_Coefficient") %>%
  rename(Correlation = value)

significant_p_value_data <- significant_heatmap_melted %>%
  filter(variable == "P_Value") %>%
  rename(P_Value = value)

# Merge correlation data with p-value data
significant_heatmap_final <- merge(significant_correlation_data, significant_p_value_data, by = c("Imgs", "Signature"))

# Create a new column to control color intensity based on p-value threshold
significant_heatmap_final$alpha_value <- ifelse(significant_heatmap_final$P_Value < p_value_threshold, "p < 0.01", "ns")

# Create the heatmap with conditional coloring and filtered signatures
ggplot(significant_heatmap_final, aes(x = Imgs, y = Signature, fill = Correlation, alpha = alpha_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkgreen", high = "orange", mid = "white", midpoint = 0, name = "Spearman Correlation") +
  scale_alpha_manual(values = c("p < 0.01" = 1, "ns" = 0.3)) +  # Adjust the transparency based on p-value
  theme_classic() +
  labs(x = "Imaging Features", y = "Signature", title = "BSW Correlation Heatmap (p<0.01)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for readability
