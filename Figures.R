library(Seurat)
library(tidyverse)
library(patchwork)
library(gtsummary)

# Table 1 Demographics
metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
metadata <- metadata %>% dplyr::select(-Sample_ID, -stanford, -acr, -amr, -high_grade_cav)
metadata <- metadata %>% dplyr::select(age, sex, dm2, everything())
metadata <- metadata %>% mutate(
  sex = case_when(
    sex == "m" ~ 1,
    TRUE ~ 0
  )
)

table_1 <- tbl_summary(metadata,
                      by = group,
                      missing = "no",
                      label = list(age ~ "Age (years)",
                                   sex ~ "Gender (male)",
                                   dm2 ~ "Diabetes mellitus",
                                   cav_grade ~ "CAV grade",
                                   hx_acr ~ "History of acute cellular rejection",
                                   hx_amr ~ "History of antibody-mediated rejection",
                                   class2_dsa ~ "Class II donor-specific antibodies",
                                   prednisone ~ "Prednisone",
                                   cni ~ "Calcineurin inhibitor",
                                   metabolite ~ "Anti-metabolite",
                                   mtori ~ "mTOR inhibitor")) %>%
  add_p() %>%
  modify_header(label = "**Demographics by CAV status**",
                stat_1 = "**Has CAV (N = 22)**",
                stat_2 = "**No CAV (N = 18)**") %>% bold_labels() 
table_1 %>%
  as_gt() %>%
  gt::gtsave(filename = "PBMC_CellRanger_6.1.2/Upload/Figures/Table_1.png")

######### Volcano plots for subclusters

###### CD4 - high-grade
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD4/High_grade"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/CD4_volcano_high_grade.png", dpi = 600)

# CD4 - CAV vs all
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD4/CAV_group"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/CD4_volcano_cav_group.png", dpi = 600)


############ CD8
### High-grade
# CD4 - high-grade
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD8/High_grade"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/CD8_volcano_high_grade.png", dpi = 600)

# CD8 - CAV vs all
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD8/CAV_group"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/CD8_volcano_cav_group.png", dpi = 600)

#################### NK cells
### High-grade
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/NK_cells/High_grade"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/NK_volcano_high_grade.png", dpi = 600)

## NK cells - CAV vs all
# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/NK_cells/CAV_group"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/NK_volcano_cav_group.png", dpi = 600)

### CAV vs non-CAV for all PBMCs

# set the path to the folder containing the csv files
path <- "PBMC_CellRanger_6.1.2/Upload/DGE/CAV_v_nonCAV"

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

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/Volcano_cav_vs_nonCAV.png", dpi = 600)
