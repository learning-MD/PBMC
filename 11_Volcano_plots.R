library(tidyverse)

# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/DGE/Prednisone_v_nonPred"

# read all csv files in the folder
files <- list.files(path, pattern = "*.csv", full.names = TRUE)

# create an empty list to store the data from each file
data_list <- list()

# loop through the files, add a cluster column, and append the data to the data_list
for (file in files) {
  data <- read_csv(file)
  data$cluster <- basename(file)
  data_list[[length(data_list) + 1]] <- data
}

# merge all files using bind_rows
combined_data <- bind_rows(data_list)
combined_data <- combined_data %>%
  mutate(cluster = str_replace_all(cluster, c("_" = " ", ".csv" = "")))

ggplot(combined_data, aes(x = logFC, y = -log10(padf), color = cluster)) +
  geom_point() +
  xlab("log2 Fold Change") +
  ylab("-log10Padj") +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed")

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/Volcano_prednisone.png", dpi = 600)

############################## High-grade CAV
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/DGE/High_grade_vs_nonCAV"

# read all csv files in the folder
files <- list.files(path, pattern = "*.csv", full.names = TRUE)

# create an empty list to store the data from each file
data_list <- list()

# loop through the files, add a cluster column, and append the data to the data_list
for (file in files) {
  data <- read_csv(file)
  data$cluster <- basename(file)
  data_list[[length(data_list) + 1]] <- data
}

# merge all files using bind_rows
combined_data <- bind_rows(data_list)
combined_data <- combined_data %>%
  mutate(cluster = str_replace_all(cluster, c("_" = " ", ".csv" = "")))

ggplot(combined_data, aes(x = logFC, y = -log10(padf), color = cluster)) +
  geom_point() +
  xlab("log2 Fold Change") +
  ylab("-log10Padj") +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed")

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/Volcano_high_grade.png", dpi = 600)

############################## mTORi - pred-adjusted
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/DGE/mTORi_v_nonmTORi"

# read all csv files in the folder
files <- list.files(path, pattern = "*.csv", full.names = TRUE)

# create an empty list to store the data from each file
data_list <- list()

# loop through the files, add a cluster column, and append the data to the data_list
for (file in files) {
  data <- read_csv(file)
  data$cluster <- basename(file)
  data_list[[length(data_list) + 1]] <- data
}

# merge all files using bind_rows
combined_data <- bind_rows(data_list)
combined_data <- combined_data %>%
  mutate(cluster = str_replace_all(cluster, c("_" = " ", ".csv" = "")))

ggplot(combined_data, aes(x = logFC, y = -log10(padf), color = cluster)) +
  geom_point() +
  xlab("log2 Fold Change") +
  ylab("-log10Padj") +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed")

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/Volcano_mTORi.png", dpi = 600)

