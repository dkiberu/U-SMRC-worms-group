
setwd("/Users/th3godfather/Desktop/Fibroschot/mar_19_2025")

library(data.table)
library(tidyverse)
library(plotly)
library(readxl)
library(VariantAnnotation)
library(sequoia)



##### 50k subset #######
# read in data
pca <- read_table2("pca.eigenvec", col_names = FALSE)
eigenval <- scan("pca.eigenval")


# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
## format ids
pca$ind <- sub("^0_", "", pca$ind)

## read meta data
metadata <- read.csv("~/Desktop/U-SMRC/sibship/merged_13_03_2025.csv")
meta <- metadata[, c("vcf_name", "Treatment.timepoint", "Location")]
colnames(meta) <- c("ind", "Treatment", "Location")
meta <- meta %>%
  distinct(ind, .keep_all = TRUE)
table(meta$Treatment)
table(meta$Location)
# 1. Ensure the column is character
meta$Treatment <- as.character(meta$Treatment)

# 2. Trim leading/trailing whitespace
meta$Treatment <- trimws(meta$Treatment)

# 3. Convert textual "NA" or "N/A" to actual NA
meta$Treatment[meta$Treatment %in% c("NA", "N/A", "")] <- NA
meta$Location[meta$Location %in% c("NA", "N/A", "")] <- NA

# 4. Standardize "pre" and "post" (case-insensitive)
meta$Treatment[tolower(meta$Treatment) == "pre"]  <- "Pre"
meta$Treatment[tolower(meta$Treatment) == "post"] <- "Post"

# 5. Convert to a factor (optional, but often helpful)
meta$Treatment <- factor(meta$Treatment, levels = c("Pre","Post"))
# Now 'NA' will remain as actual missing values
table(meta$Treatment)
table(meta$Location)
#write.csv(meta, 'combined_meta', row.names = F)

## merge pca with metadata
new_pca <- merge(x = pca, y = meta, by = 'ind')

## remake dataframe
Location<- new_pca$Location
new_pca1 <- as_tibble(data.frame(new_pca, Location))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()



# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca (Location)
b <- ggplot(new_pca1, aes(PC4, PC5, col = Location, text = ind)) + 
  geom_point(size = 1.5) +  # Increase point size
  coord_equal() + 
  theme_light()

b <- b + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  labs(colour = "Village") + 
  theme(
    axis.text = element_text(size = 12),     # Font size for axis labels
    axis.title = element_text(size = 14),   # Font size for axis titles
    legend.text = element_text(size = 12),  # Font size for legend text
    legend.title = element_text(size = 14)  # Font size for legend title
  )
b <- b + theme(aspect.ratio = 1 / 1)  # Makes the plot twice as wide as tall
#b
b + ggtitle("A PCA of S.mansoni at Lake Albert showing population structure by village")
ggplotly(b)

# plot pca (Treatment status)
# Convert treatment status to factor, replacing NA with "Unknown"
#new_pca1$Treatment <- factor(new_pca1$Treatment, levels = unique(new_pca1$Treatment))
new_pca1$Treatment[is.na(new_pca1$Treatment)] <- "Unknown"
new_pca1$Treatment <- as.character(new_pca1$Treatment)
new_pca1$Treatment[is.na(new_pca1$Treatment)] <- "Unknown"
new_pca1$Treatment <- factor(new_pca1$Treatment, levels = c("Pre", "Post", "Unknown"))

b <- ggplot(new_pca1, aes(PC4, PC5, col = Treatment, text = ind)) + 
  geom_point(size = 1.5) +  # Increase point size
  scale_colour_manual(values = c("steelblue", "orange", "gray")) +
  coord_equal() + 
  theme_light()

b <- b + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  labs(colour = "Treatment status") + 
  theme(
    axis.text = element_text(size = 12),     # Font size for axis labels
    axis.title = element_text(size = 14),   # Font size for axis titles
    legend.text = element_text(size = 12),  # Font size for legend text
    legend.title = element_text(size = 14)  # Font size for legend title
  )
#b
b <- b + theme(aspect.ratio = 1 / 1)  # Makes the plot twice as wide as tall
b + ggtitle("A PCA of S.mansoni at Lake Albert showing population structure by treatment status")
ggplotly(b)



### view pcas from different pcs (Location)
generate_pca_plot <- function(data, pve, pc1, pc2) {
  # Extract numeric part from PC1 and PC2
  pc1_number <- as.numeric(gsub("\\D", "", pc1))
  pc2_number <- as.numeric(gsub("\\D", "", pc2))
  
  #ggplotly(
  ggplot(data, aes(get(pc2), get(pc1), col = Location, text = ind)) + 
    geom_point(size = 1.5) + 
    
    theme_light() + 
    xlab(paste0(pc2, " (", signif(pve$pve[pc2_number], 3), "%)")) + 
    ylab(paste0(pc1, " (", signif(pve$pve[pc1_number], 3), "%)")) + 
    ggtitle(paste("Population structure by village,", pc1, "and", pc2)) + theme(text=element_text(size=20)) +
    labs(colour="Village")
  #)
}

# Generate PCA plots for different pairs of principal components
pc_pairs <- expand.grid(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
pc_pairs <- pc_pairs[pc_pairs$Var1 != pc_pairs$Var2, ]  # Exclude pairs with the same PC
pc_pairs <- pc_pairs[!duplicated(t(apply(pc_pairs, 1, sort))), ]  # Exclude reverse combinations
for (i in 1:nrow(pc_pairs)) {
  pc1 <- as.character(pc_pairs[i, 1])
  pc2 <- as.character(pc_pairs[i, 2])
  
  print(generate_pca_plot(new_pca1, pve, pc1, pc2))
}


### view pcas from different pcs (Treatment)
generate_pca_plot <- function(data, pve, pc1, pc2) {
  # Extract numeric part from PC1 and PC2
  pc1_number <- as.numeric(gsub("\\D", "", pc1))
  pc2_number <- as.numeric(gsub("\\D", "", pc2))
  
  ggplotly(
    ggplot(data, aes(get(pc2), get(pc1), col = Treatment, text = ind)) + 
      scale_colour_manual(values = c("steelblue", "orange", "gray")) +
      geom_point(size = 3) + 
      
      theme_light() + 
      xlab(paste0(pc2, " (", signif(pve$pve[pc2_number], 3), "%)")) + 
      ylab(paste0(pc1, " (", signif(pve$pve[pc1_number], 3), "%)")) + 
      ggtitle(paste("Population structure by treatment status,", pc1, "and", pc2)) + theme(text=element_text(size=20)) +
      labs(colour="Treatment status")
  ) 
}

# Generate PCA plots for different pairs of principal components
pc_pairs <- expand.grid(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
pc_pairs <- pc_pairs[pc_pairs$Var1 != pc_pairs$Var2, ]  # Exclude pairs with the same PC
pc_pairs <- pc_pairs[!duplicated(t(apply(pc_pairs, 1, sort))), ]  # Exclude reverse combinations
for (i in 1:nrow(pc_pairs)) {
  pc1 <- as.character(pc_pairs[i, 1])
  pc2 <- as.character(pc_pairs[i, 2])
  
  print(generate_pca_plot(new_pca1, pve, pc1, pc2))
}



##### manhattan plots pi #########
## non-overlapping windows
library(dplyr)
library(ggplot2)

df <- read.table("pi_non_overlapping.windowed.pi", header = TRUE)

desired_order <- c("SM_V10_1", "SM_V10_2", "SM_V10_3", "SM_V10_4", "SM_V10_5", "SM_V10_6", "SM_V10_7", "SM_V10_WSR", "SM_V10_Z", "SM_V10_MITO")
df$CHROM <- factor(df$CHROM, levels = desired_order)

df <- df %>%
  mutate(bin_mid = (BIN_START + BIN_END) / 2)

chr_sizes <- df %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(BIN_END)) %>%
  arrange(match(CHROM, desired_order)) %>%
  mutate(cum_len = cumsum(as.numeric(chr_len)) - chr_len)

df <- df %>%
  left_join(chr_sizes, by = "CHROM") %>%
  mutate(pos_cum = bin_mid + cum_len) %>%
  arrange(match(CHROM, desired_order), bin_mid)

axisdf <- df %>%
  group_by(CHROM) %>%
  summarize(center = mean(range(pos_cum)))

a <- ggplot(df, aes(x = pos_cum, y = PI, color = CHROM)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  scale_y_continuous(labels = scales::label_number(trim = TRUE)) +
  scale_color_manual(values = rep(c("steelblue", "tomato"), length(unique(df$CHROM)))) +
  labs(x = "Chromosome", y = expression(pi), title = "Manhattan plot of Nucleotide Diversity (\u03C0), Non-overlapping windows") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
a


# with thresholds
low_threshold <- 0.001  # Known biologically meaningful threshold
#high_threshold <- 0.01

# Generate the Manhattan plot with thresholds
a <- ggplot(df, aes(x = pos_cum, y = PI, color = CHROM)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_hline(yintercept = low_threshold, color = "green", linetype = "dashed", size = 1, label = "Low Diversity") +
  #geom_hline(yintercept = high_threshold, color = "red", linetype = "dashed", size = 1, label = "High Diversity") +
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  scale_color_manual(values = rep(c("steelblue", "tomato"), length(unique(df$CHROM)))) +
  labs(x = "Chromosome", 
       y = expression(pi), 
       title = "Manhattan plot of Nucleotide Diversity (\u03C0), Non-overlapping windows") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(a)

## overlapping windows
df <- read.table("pi_overlapping.windowed.pi", header = TRUE)

desired_order <- c("SM_V10_1", "SM_V10_2", "SM_V10_3", "SM_V10_4", "SM_V10_5", "SM_V10_6", "SM_V10_7", "SM_V10_WSR", "SM_V10_Z", "SM_V10_MITO")
df$CHROM <- factor(df$CHROM, levels = desired_order)

df <- df %>%
  mutate(bin_mid = (BIN_START + BIN_END) / 2)

chr_sizes <- df %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(BIN_END)) %>%
  arrange(match(CHROM, desired_order)) %>%
  mutate(cum_len = cumsum(as.numeric(chr_len)) - chr_len)

df <- df %>%
  left_join(chr_sizes, by = "CHROM") %>%
  mutate(pos_cum = bin_mid + cum_len) %>%
  arrange(match(CHROM, desired_order), bin_mid)

axisdf <- df %>%
  group_by(CHROM) %>%
  summarize(center = mean(range(pos_cum)))

ggplot(df, aes(x = pos_cum, y = PI, color = CHROM)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  scale_color_manual(values = rep(c("steelblue", "tomato"), length(unique(df$CHROM)))) +
  labs(x = "Chromosome", y = expression(pi), title = "Manhattan plot of Nucleotide Diversity (\u03C0), overlapping windows") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

## with thresholds
#low_threshold <- 0.001  # Known biologically meaningful threshold
#high_threshold <- 0.01

# Generate the Manhattan plot with thresholds
ggplot(df, aes(x = pos_cum, y = PI, color = CHROM)) +
  geom_point(alpha = 0.8, size = 1.2) +
  #geom_hline(yintercept = low_threshold, color = "green", linetype = "dashed", size = 1, label = "Low Diversity") +
  #geom_hline(yintercept = high_threshold, color = "red", linetype = "dashed", size = 1, label = "High Diversity") +
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  scale_color_manual(values = rep(c("steelblue", "tomato"), length(unique(df$CHROM)))) +
  labs(x = "Chromosome", 
       y = expression(pi), 
       title = "Manhattan plot of Nucleotide Diversity (\u03C0), overlapping windows") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


#### fst ####
### villages
library(dplyr)
library(ggplot2)

fst <- read.table("fst_villages.weir.fst", header = TRUE)

# Explicitly set the chromosome order, moving "SM_V10_MITO" to the last
chromosome_order <- c("SM_V10_1", "SM_V10_2", "SM_V10_3", "SM_V10_4", 
                      "SM_V10_5", "SM_V10_6", "SM_V10_7", "SM_V10_WSR", 
                      "SM_V10_Z", "SM_V10_MITO")

fst <- fst %>%
  mutate(CHROM = factor(CHROM, levels = chromosome_order)) %>%
  arrange(CHROM, POS)

# Calculate cumulative positions for chromosomes
chr_sizes <- fst %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS)) %>%
  arrange(CHROM) %>%
  mutate(cum_len = cumsum(chr_len) - chr_len)

fst <- fst %>%
  left_join(chr_sizes, by = "CHROM") %>%
  mutate(pos_cum = POS + cum_len) %>%
  arrange(CHROM, POS)

# Define axis labels and positions
axisdf <- fst %>%
  group_by(CHROM) %>%
  summarize(center = mean(range(pos_cum)))

# Plot the Fst Manhattan plot
library(ggplot2)
library(plotly)

max_pos <- max(fst$WEIR_AND_COCKERHAM_FST)
# Update ggplot to include POS in the tooltip
a <- ggplot(fst, aes(x = pos_cum, y = WEIR_AND_COCKERHAM_FST, 
                     color = as.factor(CHROM), 
                     text = paste0("Chromosome: ", CHROM, 
                                   "<br>Position: ", POS, 
                                   "<br>Fst: ", round(WEIR_AND_COCKERHAM_FST, 4)))) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  labs(x = "Chromosome", y = "Fst", title = "Manhattan Plot of Fst (Villages)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(limits = c(0, max_pos))
a
# Convert to ggplotly, showing POS in the tooltip
ggplotly(a, tooltip = "text")


#################################################################################
library(ggplot2)
library(dplyr)
library(ggrepel)

# Filter points with Fst > 0.05
fst_high <- fst %>%
  filter(WEIR_AND_COCKERHAM_FST >= 0.05)

# Create the plot
a <- ggplot(fst, aes(x = pos_cum, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_point(alpha = 0.8, size = 1.5) + # Plot all points
  geom_text_repel(data = fst_high, 
                  aes(label = paste0(CHROM, ":", POS)), # Label with chromosome and position
                  size = 3, 
                  color = "black",
                  max.overlaps = 10, # Control overlap
                  box.padding = 0.3) + # Add space around labels
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  labs(x = "Chromosome", y = "Fst", title = "Manhattan Plot of Fst (Villages)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(limits = c(0, max_pos))

print(a)
#################################################################################


## treatment status
fst <- read.table("fst_treatment.weir.fst", header = TRUE)

# Explicitly set the chromosome order, moving "SM_V10_MITO" to the last
chromosome_order <- c("SM_V10_1", "SM_V10_2", "SM_V10_3", "SM_V10_4", 
                      "SM_V10_5", "SM_V10_6", "SM_V10_7", "SM_V10_WSR", 
                      "SM_V10_Z", "SM_V10_MITO")

fst <- fst %>%
  mutate(CHROM = factor(CHROM, levels = chromosome_order)) %>%
  arrange(CHROM, POS)

# Calculate cumulative positions for chromosomes
chr_sizes <- fst %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS)) %>%
  arrange(CHROM) %>%
  mutate(cum_len = cumsum(chr_len) - chr_len)

fst <- fst %>%
  left_join(chr_sizes, by = "CHROM") %>%
  mutate(pos_cum = POS + cum_len) %>%
  arrange(CHROM, POS)

# Define axis labels and positions
axisdf <- fst %>%
  group_by(CHROM) %>%
  summarize(center = mean(range(pos_cum)))

max_pos <- max(fst$WEIR_AND_COCKERHAM_FST)
# Update ggplot to include POS in the tooltip
a <- ggplot(fst, aes(x = pos_cum, y = WEIR_AND_COCKERHAM_FST, 
                     color = as.factor(CHROM), 
                     text = paste0("Chromosome: ", CHROM, 
                                   "<br>Position: ", POS, 
                                   "<br>Fst: ", round(WEIR_AND_COCKERHAM_FST, 4)))) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  labs(x = "Chromosome", y = "Fst", title = "Manhattan Plot of Fst (Treatment status)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(limits = c(0, max_pos))
a
# Convert to ggplotly, showing POS in the tooltip
ggplotly(a, tooltip = "text")


#################################################################################
library(ggplot2)
library(dplyr)
library(ggrepel)

# Filter points with Fst > 0.05
fst_high <- fst %>%
  filter(WEIR_AND_COCKERHAM_FST > 0.05)

# Create the plot
a <- ggplot(fst, aes(x = pos_cum, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_point(alpha = 0.8, size = 1.5) + # Plot all points
  geom_text_repel(data = fst_high, 
                  aes(label = paste0(CHROM, ":", POS)), # Label with chromosome and position
                  size = 3, 
                  color = "black",
                  max.overlaps = 10, # Control overlap
                  box.padding = 0.3) + # Add space around labels
  scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
  labs(x = "Chromosome", y = "Fst", title = "Manhattan Plot of Fst (Treatment)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(limits = c(0, max_pos))

print(a)
#################################################################################


#### ternary plot ####
################# Ternary plot ########################
ngsrelate <- read.table("ngsrelate_results.res", header=TRUE, sep="\t")
# Read the VCF file
vcf <- readVcf("subset_main_50k.vcf.gz")
head(ngsrelate)
head(ngsrelate[,c(1:6)])
#metadata <- read_excel("~/Desktop/U-SMRC/sibship/all_samples_info.10July.xls")
#metadata <- read.csv("~/Desktop/U-SMRC/sibship/merged_13_03_2025.csv")
#meta <- metadata[,c(5,10,11)]
names(meta) <- c("vcf_name","Treatment timepoint","Location")
meta$`Treatment timepoint` <- ifelse(grepl("pre", meta$`Treatment timepoint`, ignore.case = TRUE), "Pre",
                                     ifelse(grepl("post", meta$`Treatment timepoint`, ignore.case = TRUE), "Post",
                                            meta$`Treatment timepoint`))
meta$`Treatment timepoint`[meta$`Treatment timepoint` == "N/A"] <- NA
meta$`Treatment timepoint`[meta$`Location` == "N/A"] <- NA
meta$`Treatment timepoint`[meta$`Treatment timepoint` == ""] <- NA
# Check the updated counts
table(meta$`Treatment timepoint`)
table(meta$Location)
head(meta)

# Extract sample names from the VCF file
sample_names <- colnames(geno(vcf)$GT)


# Replace the `a` and `b` columns with sample names
ngsRelate <- ngsrelate %>%
  mutate(a = sample_names[a + 1],  # Adding 1 to match 1-based indexing in R
         b = sample_names[b + 1])  # Adding 1 to match 1-based indexing in R
#write.csv(ngsRelate, file = "50k_named.res")
head(ngsRelate)

# First, merge the meta data with ngsRelate based on the `a` column
ngsRelate <- ngsRelate %>%
  left_join(meta, by = c("a" = "vcf_name"), relationship = "many-to-many")

# Rename the new columns to avoid confusion
names(ngsRelate)[names(ngsRelate) == "Treatment timepoint"] <- "Treatment_timepoint_a"
names(ngsRelate)[names(ngsRelate) == "Location"] <- "Location_a"

#Repeat the merge for the `b` column
ngsRelate <- ngsRelate %>%
  left_join(meta, by = c("b" = "vcf_name"), relationship = "many-to-many")

# Rename the new columns
names(ngsRelate)[names(ngsRelate) == "Treatment timepoint"] <- "Treatment_timepoint_b"
names(ngsRelate)[names(ngsRelate) == "Location"] <- "Location_b"

## S1-Buhirigi, S2-Kaiso
ngsRelate <- ngsRelate %>%
  mutate(
    Location_a = case_when(
      is.na(Location_a) & substr(a, 1, 2) == "S1" ~ "Buhirigi",
      is.na(Location_a) & substr(a, 1, 2) == "S2" ~ "Kaiso",
      TRUE ~ Location_a
    ),
    Location_b = case_when(
      is.na(Location_b) & substr(b, 1, 2) == "S1" ~ "Buhirigi",
      is.na(Location_b) & substr(b, 1, 2) == "S2" ~ "Kaiso",
      TRUE ~ Location_b
    )
  )

# Create the new column based on the condition
ngsRelate <- ngsRelate %>%
  mutate(Location = ifelse(Location_a == Location_b, "Same village", "Different villages"))

# Extract base names using a more flexible pattern
ngsRelate$base_a <- sub("[-_][^-_]*$", "", ngsRelate$a)
ngsRelate$base_a <- sub("[_][^_]*$", "", ngsRelate$base_a)
ngsRelate$base_b <- sub("[-_][^-_]*$", "", ngsRelate$b)
ngsRelate$base_b <- sub("[_][^_]*$", "", ngsRelate$base_b)


# Create the Donor column based on the comparison of base_a and base_b
ngsRelate$Donor <- ifelse(ngsRelate$base_a == ngsRelate$base_b, "Same donor", "Different donors")

# Optional: Remove the temporary base_a and base_b columns if you no longer need them
#ngsRelate <- ngsRelate[, !(names(ngsRelate) %in% c("base_a", "base_b"))]

# View the updated data frame
head(ngsRelate)

# Load the necessary libraries
library(ggtern)
library(ggplot2)

#############
# Define a function to categorize points based on the given conditions
# categorize_relationship <- function(j9, j8, j7, theta) {
#   if (theta >= 0.177 & j8 >0.05 ) {
#     return("First Degree")
#  } else if (theta >= 0.0884 & theta <0.177 & j8 > 0.05)  {
#     return("Second Degree")  # e.g. Uncle/Aunt-niece/nephew, Half-siblings, Grandparent-Grandchild
#  }else if (theta >= 0.0442 & theta < 0.0884 & j8 > 0.05)  {
#     return("Third Degree")  
#  } else {
#     return("Other")  # For points that don't match any category
#  }
#}

### Original classification above 👆
# categorize_relationship <- function(j9, j8, j7, theta) {
#   if (theta >= 0.18 & j8 > 0.10) {
#     return("First Degree")  # e.g., Parent-Offspring, Full-Siblings
#   } else if (theta >= 0.09 & theta < 0.18 & j8 > 0.05) {
#     return("Second Degree")  # e.g., Half-Siblings, Grandparent-Grandchild
#   } else if (theta >= 0.045 & theta < 0.09 & j8 > 0.03) {
#     return("Third Degree")  # e.g., First Cousins
#   } else {
#     return("Other")  # Distant or unrelated individuals
#   }
# }


## staples http://dx.doi.org/10.1016/j.ajhg.20
categorize_relationship <- function(j9, j8, j7) {
  if (j9 >= 0.6 & j9 <= 0.9 & j8 <= 0.3 & j8 >= 0.1 & j7 <= 0.2) {
    return("Third Degree")  # First Cousins, Great-Grandparental, Half-Avuncular
  } else if (j9 >= 0.1 & j9 < 0.55 & j8 >= 0.2 & j8 <= 0.7 & j7 >= 0.1) {
    return("First Degree")  # Parent-Offspring, Full-Siblings
  } else if (j9 > 0.48 & j9 <= 0.7 & j8 > 0.2 & j8 < 0.7 & j7 <= 0.4) {
    return("Second Degree")  # Half-Siblings, Avuncular, Grandparent-Grandchild
  } else {
    return("Other")  # Unrelated or distant relationships
  }
}

## staples http://dx.doi.org/10.1016/j.ajhg.20
categorize_relationship <- function(j9, j8, j7) {
  if (j9 >= 0.1 & j9 < 0.6 & j8 >= 0.2 & j8 <= 0.6 & j7 >= 0.15) {
    return("First Degree")  # Parent-Offspring, Full-Siblings
  } else if (j9 > 0.4 & j9 <= 0.65 & j8 >= 0.2 & j8 < 0.5 & j7 <= 0.15) {
    return("Second Degree")  # Half-Siblings, Avuncular, Grandparent-Grandchild
  } else if (j9 >= 0.6 & j8 <= 0.3 & j8 >= 0.1 & j7 <= 0.2) {
      return("Third Degree")  # First Cousins, Great-Grandparental, Half-Avuncular
  } else {
    return("Other")  # Unrelated or distant relationships
  }
}
##
# Apply the function to create a new column in the data
#ngsRelate$category <- mapply(categorize_relationship, ngsRelate$J9, ngsRelate$J8, ngsRelate$J7)
ngsRelate$category <- mapply(categorize_relationship, ngsRelate$J9, ngsRelate$J8, ngsRelate$J7)

# Ensure the category column is a factor
ngsRelate$category <- factor(ngsRelate$category, levels = c("First Degree", "Second Degree", "Third Degree", "Other"))
table(ngsRelate$category)
# Create the ternary plot
ternary_plot <- ggtern(data = ngsRelate, aes(x = J9, y = J8, z = J7, shape = Donor, color = category)) +
  geom_point() +
  theme_bw() +
  labs(
    #title = "Lake Albert S.mansoni genetic relatedness analysis",
    x = "k0",
    y = "k1",
    z = "k2",
    color = "Relationship"
  ) +
  scale_color_manual(values = c("First Degree" = "#974181", "Second Degree" = "#d57946" , "Third Degree" = "#7daf5c", "Other" = "#4e9ab7")) +
  theme_showarrows() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white"),  # White background for countries
    plot.background = element_rect(fill = "white", color = NA),  # Skyblue for water
    plot.title = element_text(size = 24, face = "bold"),  # Title
    axis.title.x = element_text(size = 20),  # X-axis title
    axis.title.y = element_text(size = 20),  # Y-axis title
    
    axis.text = element_text(size = 18),  # Axis labels
    legend.title = element_text(size = 20),  # Legend title
    legend.text = element_text(size = 18)  # Legend text
  )

# Display the plot
print(ternary_plot)

table(ngsRelate[ngsRelate$base_a==ngsRelate$base_b & ngsRelate$Treatment_timepoint_a!=ngsRelate$Treatment_timepoint_b,'category'])


############################ kaiso ###############################
ngsRelate_k <- ngsRelate[ngsRelate$Location_a=='Kaiso' & ngsRelate$Location_b=='Kaiso',]
# Apply the function to create a new column in the data
#ngsRelate_k$category <- mapply(categorize_relationship, ngsRelate_k$J9, ngsRelate_k$J8, ngsRelate_k$J7)
ngsRelate_k$category <- mapply(categorize_relationship, ngsRelate_k$J9, ngsRelate_k$J8, ngsRelate_k$J7)

# Ensure the category column is a factor
ngsRelate_k$category <- factor(ngsRelate_k$category, levels = c("First Degree", "Second Degree", "Third Degree", "Other"))

# Create the ternary plot
ternary_plot <- ggtern(data = ngsRelate_k, aes(x = J9, y = J8, z = J7, shape = Donor, color = category)) +
  geom_point() +
  theme_bw() +
  labs(
    #title = "Lake Albert S.mansoni genetic relatedness analysis",
    x = "k0",
    y = "k1",
    z = "k2",
    color = "Relationship"
  ) +
  scale_color_manual(values = c("First Degree" = "#974181", "Second Degree" = "#d57946" , "Third Degree" = "#7daf5c", "Other" = "#4e9ab7")) +
  theme_showarrows() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white"),  # White background for countries
    plot.background = element_rect(fill = "white", color = NA),  # Skyblue for water
    plot.title = element_text(size = 24, face = "bold"),  # Title
    axis.title.x = element_text(size = 20),  # X-axis title
    axis.title.y = element_text(size = 20),  # Y-axis title
    
    axis.text = element_text(size = 18),  # Axis labels
    legend.title = element_text(size = 20),  # Legend title
    legend.text = element_text(size = 18)  # Legend text
  )

# Display the plot
print(ternary_plot)


############################ Buhirigi ###############################
ngsRelate_b <- ngsRelate[ngsRelate$Location_a=='Buhirigi' & ngsRelate$Location_b=='Buhirigi',]
# Apply the function to create a new column in the data
#ngsRelate$category <- mapply(categorize_relationship, ngsRelate_b$J9, ngsRelate_b$J8, ngsRelate_b$J7)
ngsRelate_b$category <- mapply(categorize_relationship, ngsRelate_b$J9, ngsRelate_b$J8, ngsRelate_b$J7)

# Ensure the category column is a factor
ngsRelate_b$category <- factor(ngsRelate_b$category, levels = c("First Degree", "Second Degree", "Third Degree", "Other"))

# Create the ternary plot
ternary_plot <- ggtern(data = ngsRelate_b, aes(x = J9, y = J8, z = J7, shape = Donor, color = category)) +
  geom_point() +
  theme_bw() +
  labs(
    #title = "Lake Albert S.mansoni genetic relatedness analysis",
    x = "k0",
    y = "k1",
    z = "k2",
    color = "Relationship"
  ) +
  scale_color_manual(values = c("First Degree" = "#974181", "Second Degree" = "#d57946" , "Third Degree" = "#7daf5c", "Other" = "#4e9ab7")) +
  theme_showarrows() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white"),  # White background for countries
    plot.background = element_rect(fill = "white", color = NA),  # Skyblue for water
    plot.title = element_text(size = 24, face = "bold"),  # Title
    axis.title.x = element_text(size = 20),  # X-axis title
    axis.title.y = element_text(size = 20),  # Y-axis title
    
    axis.text = element_text(size = 18),  # Axis labels
    legend.title = element_text(size = 20),  # Legend title
    legend.text = element_text(size = 18)  # Legend text
  )

# Display the plot
print(ternary_plot)


############## circos plot ###################
# Load required libraries
library(dplyr)
library(circlize)
library(ComplexHeatmap)

### 1. Preprocess Data and Create Sample Metadata ###

# Subset the columns and set treatment status
circos_df <- ngsRelate[, c("a","b","Treatment_timepoint_a","Treatment_timepoint_b","Location_a","Location_b","Donor","category","base_a","base_b")]
circos_df <- circos_df %>%
  mutate(Treatment = if_else(Treatment_timepoint_a == Treatment_timepoint_b,
                             "Same", 
                             "Different"))

# Remove rows with missing data
dim(circos_df)
#circos_df <- circos_df[complete.cases(circos_df), ]
dim(circos_df)
# filter out "Other" category
circos_df <- circos_df %>% filter(category != "Other")

# Create base identifiers by stripping the trailing part from sample IDs
# circos_df$base_a <- ifelse(grepl("_", circos_df$a),
#                            sub("_.*", "", circos_df$a),
#                            sub("-[A-Za-z]+$", "", circos_df$a))
# circos_df$base_b <- ifelse(grepl("_", circos_df$b),
#                            sub("_.*", "", circos_df$b),
#                            sub("-[A-Za-z]+$", "", circos_df$b))

circos_df <- circos_df %>%
  mutate(Treatment = if_else(Treatment_timepoint_a == Treatment_timepoint_b, "Same", "Different")) %>%
  filter(!is.na(a), !is.na(b), category != "Other")

# Define color mappings
treatment_colors   <- c("Pre" = "red", "Post" = "green", "Unknown" = "gray")
relatedness_colors <- c("First Degree" = "#974181", 
                        "Second Degree" = "#d57946", 
                        "Third Degree" = "steelblue")
village_colors     <- c("Buhirigi" = "purple", "Kaiso" = "blue")

# Create a mapping from sample to its base name (individual)
sample_to_base <- circos_df %>%
  dplyr::distinct(a, base_a) %>%
  dplyr::rename(sample = a, base_name = base_a) %>%
  dplyr::bind_rows(
    circos_df %>%
      dplyr::distinct(b, base_b) %>%
      dplyr::rename(sample = b, base_name = base_b)
  ) %>%
  distinct()

# Create a mapping from sample to its village location (from both columns)
location_info <- bind_rows(
  circos_df %>% dplyr::select(sample = a, location = Location_a),
  circos_df %>% dplyr::select(sample = b, location = Location_b)
) %>% distinct()

# Combine the above into sample metadata
sample_meta <- sample_to_base %>%
  left_join(location_info, by = "sample")

# Create treatment info from both columns
treatment_info <- bind_rows(
  circos_df %>% dplyr::select(sample = a, treatment = Treatment_timepoint_a),
  circos_df %>% dplyr::select(sample = b, treatment = Treatment_timepoint_b)
) %>% dplyr::distinct() %>%
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# Merge treatment info into sample metadata.
sample_meta <- sample_meta %>%
  left_join(treatment_info, by = "sample") %>%
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment),
         treatment = factor(treatment, levels = c("Pre", "Post", "Unknown")))

# Order individuals by village first, then by sample
sample_meta <- sample_meta %>%
  arrange(location, sample) %>%
  mutate(individual = factor(sample, levels = unique(sample)))

# --- Compute Overall Treatment Status per Individual ---
# For each individual (base_name), if all samples have the same treatment then use that;
individual_status <- sample_meta %>%
  mutate(treatment = as.character(treatment)) %>%
  group_by(base_name) %>%
  summarise(overall_treatment = if_else(n_distinct(treatment) == 1,
                                        treatment[1],
                                        "Unknown"),
            .groups = "drop") %>%
  mutate(overall_treatment = factor(overall_treatment,
                                    levels = c("Post", "Pre", "Unknown")))

# Merge the overall treatment status into sample metadata
sample_meta <- sample_meta %>%
  left_join(individual_status, by = "base_name")

# Order individuals by village, then by overall treatment status, then by individual (base_name) and sample
sample_meta <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  mutate(individual = factor(base_name, levels = unique(base_name)))

# Extract the ordered sample names (sectors for the circos plot)
ordered_samples <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  pull(sample)

### 2. Initialize Circos and Add Inner Tracks ###

# Clear any existing circos plot and set parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 0.5, 
           track.margin = c(0.01, 0.01),
           cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = ordered_samples, xlim = c(0, 1))

# Outer ring: Village track (colored by village)
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    s <- get.cell.meta.data("sector.index")
    v <- sample_meta$location[sample_meta$sample == s]
    circos.rect(0, 0, 1, 1, col = village_colors[v], border = NA)
  },
  bg.border    = NA,
  track.height = 0.05
)

# Second ring: Treatment track (colored by treatment)
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    s <- get.cell.meta.data("sector.index")
    t <- as.character(sample_meta$treatment[sample_meta$sample == s])
    circos.rect(0, 0, 1, 1, col = treatment_colors[t], border = NA)
  },
  bg.border    = NA,
  track.height = 0.05
)

### 3. Draw the Relatedness Links ###

for (i in seq_len(nrow(circos_df))) {
  link_type <- ifelse(circos_df$base_a[i] == circos_df$base_b[i], 1, 3)
  circos.link(
    circos_df$a[i], 0.5,
    circos_df$b[i], 0.5,
    col = relatedness_colors[circos_df$category[i]],
    lty = link_type,
    border = NA
  )
}

### 4. Add Outer Layer Spanning Samples from the Same Individual with Inward-Facing Ticks ###

# Compute group boundaries: for each individual (base_name), find the first and last sample
# group_boundaries <- sample_to_base %>%
#   arrange(base_name, sample) %>%
#   group_by(base_name) %>%
#   summarise(
#     first_sample = first(sample),
#     last_sample  = last(sample)
#   ) %>%
#   ungroup()

# group_boundaries <- sample_to_base %>%
#   arrange(base_name, sample) %>%
#   group_by(base_name) %>%
#   summarise(
#     first_sample = sample[1],
#     last_sample  = sample[length(sample)],
#     .groups = "drop"
#   )

group_boundaries <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  group_by(base_name) %>%
  summarise(
    first_sample = dplyr::first(sample),
    last_sample  = dplyr::last(sample),
    .groups = "drop"
  )
# Add an extra (empty) track for drawing the outer group arcs
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {})

# Get the current track index for the outer annotation
current_track <- get.current.track.index()

# For each group (individual), draw an arc spanning from the first to the last sample,
# and add ticks at both boundaries facing inward.
for(i in seq_len(nrow(group_boundaries))) {
  first_sector <- group_boundaries$first_sample[i]
  last_sector  <- group_boundaries$last_sample[i]
  
  # Get angular boundaries (in degrees) for the first and last sectors
  start_angle <- get.cell.meta.data("cell.start.degree", sector.index = first_sector, track.index = current_track)
  end_angle   <- get.cell.meta.data("cell.end.degree", sector.index = last_sector, track.index = current_track)
  
  # Set the radius for the outer arc
  r <- 1.05
  
  # Create a sequence of angles (in radians) for the arc
  theta <- seq(start_angle, end_angle, length.out = 100) * pi/180
  x_arc <- r * cos(theta)
  y_arc <- r * sin(theta)
  
  # Draw the arc spanning the samples of this individual
  lines(x_arc, y_arc, col = "black", lwd = 2)
  
  # Define tick length (adjust this value as needed)
  tick_length <- 0.05
  
  # Convert boundary angles to radians
  theta_start <- start_angle * pi / 180
  theta_end   <- end_angle * pi / 180
  
  # Draw tick at the beginning of the span (tick now drawn inward)
  lines(
    x = c(r * cos(theta_start), (r - tick_length) * cos(theta_start)),
    y = c(r * sin(theta_start), (r - tick_length) * sin(theta_start)),
    col = "black", lwd = 2
  )
  
  # Draw tick at the end of the span (tick now drawn inward)
  lines(
    x = c(r * cos(theta_end), (r - tick_length) * cos(theta_end)),
    y = c(r * sin(theta_end), (r - tick_length) * sin(theta_end)),
    col = "black", lwd = 2
  )
  
  # (Optional) To add labels for the individual, uncomment and adjust the following:
  # mid_angle <- (start_angle + end_angle) / 2 * pi/180
  # x_text <- (r - tick_length - 0.02) * cos(mid_angle)
  # y_text <- (r - tick_length - 0.02) * sin(mid_angle)
  # text(x_text, y_text, labels = group_boundaries$base_name[i],
  #      srt = (start_angle + end_angle) / 2, cex = 0.8, adj = c(0.5, 0.5))
}

### 5. Finalize the Plot and Add Legends ###

# (Optional) If you wish to clear the circos plot at the end, uncomment the next line.
# circos.clear()

lgd_village <- Legend(
  title    = "Village",
  at       = names(village_colors),
  legend_gp = gpar(fill = village_colors)
)
lgd_treatment <- Legend(
  title    = "Treatment status",
  at       = names(treatment_colors),
  legend_gp = gpar(fill = treatment_colors)
)
lgd_relatedness <- Legend(
  title    = "Relatedness (Links)",
  at       = names(relatedness_colors),
  legend_gp = gpar(fill = relatedness_colors)
)
lgd_links <- Legend(
  title    = "Donor Relationship (Line Style)",
  at       = c("Same Donor", "Different Donors"),
  type     = "lines",
  legend_gp = gpar(col = "black", lty = c(1, 3), lwd = 2)
)

all_legends <- packLegend(lgd_village, lgd_treatment, lgd_relatedness, lgd_links,
                          direction = "vertical")
draw(all_legends, x = unit(1, "npc") - unit(1, "cm"), just = "right")

## network
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)

# Create edges from circos_df
edges <- circos_df %>%
  select(from = a, to = b, category, Donor, Treatment)

# Create nodes from unique samples
nodes <- data.frame(name = unique(c(edges$from, edges$to)))

# Create node-level attributes: treatment and village
node_info <- bind_rows(
  circos_df %>% select(sample = a, treatment = Treatment_timepoint_a, village = Location_a),
  circos_df %>% select(sample = b, treatment = Treatment_timepoint_b, village = Location_b)
) %>% distinct() %>%
  # Change NA treatment values to "Unknown"
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# Build the graph and join node attributes
g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  left_join(node_info, by = c("name" = "sample"))

# Define color palettes for edges, treatment, and village
edge_colors <- c("First Degree" = "#974181", "Second Degree" = "#d57946", "Third Degree" = "steelblue")
node_treatment_colors <- c("Pre" = "red", "Post" = "green", "Unknown" = "gray")
node_village_colors   <- c("Buhirigi" = "purple", "Kaiso" = "blue")

set.seed(42)
ggraph(g, layout = "fr") +
  # Edges: map category to color and Donor to line type
  geom_edge_link(aes(color = category, linetype = Donor), width = 1, show.legend = TRUE) +
  scale_edge_color_manual(values = edge_colors) +
  scale_edge_linetype_manual(values = c("Same donor" = "solid", "Different donors" = "dotted")) +
  # Nodes: fill by treatment and border color by village
  geom_node_point(aes(fill = treatment, color = village), size = 3, shape = 21,stroke=2) +
  scale_fill_manual(values = node_treatment_colors) +
  scale_color_manual(values = node_village_colors) +
  theme_graph(base_size = 14) +
  guides(
    edge_color = guide_legend(title = "Relatedness"),
    edge_linetype = guide_legend(title = "Donor"),
    fill = guide_legend(title = "Treatment status"),
    color = guide_legend(title = "Village")
  ) +
  ggtitle("Network Graph of Samples with Village Information")

table(circos_df[,c('Donor','Treatment','category')])

### circos - Kaiso ###
# Load required libraries
library(dplyr)
library(circlize)
library(ComplexHeatmap)

### 1. Data Preparation ###
# Subset and prepare the data (adjust column indices/names as needed)
circos_df <- ngsRelate_k[, c("a","b","Treatment_timepoint_a","Treatment_timepoint_b","Location_a","Location_b","Donor","category","base_a","base_b")]
circos_df <- circos_df %>%
  mutate(Treatment = if_else(Treatment_timepoint_a == Treatment_timepoint_b,
                             "Same", 
                             "Different"))

# Keep only rows with non-missing location (either Location_a or Location_b)
circos_df <- circos_df[ ( !is.na(circos_df$Location_a) | !is.na(circos_df$Location_b) ), ]
# Remove rows with category "Other"
circos_df <- circos_df %>% filter(category != "Other")

# Compute the base names for each sample (from columns 'a' and 'b')
circos_df$base_a <- sub("[-_][^-_]*$", "", circos_df$a)
circos_df$base_a <- sub("[_][^_]*$", "", circos_df$base_a)
circos_df$base_b <- sub("[-_][^-_]*$", "", circos_df$b)
circos_df$base_b <- sub("[_][^_]*$", "", circos_df$base_b)
# (Filter again for category != "Other" if needed)
circos_df <- circos_df %>% filter(category != "Other")

# Define color mappings
treatment_colors   <- c("Pre" = "red", "Post" = "green", "Unknown" = "gray")
relatedness_colors <- c("First Degree"  = "#974181", 
                        "Second Degree" = "#d57946",
                        "Third Degree"  = "steelblue")

# Create a mapping from sample to its base name (here the base name represents the individual)
sample_to_base <- circos_df %>%
  distinct(a, base_a) %>%
  rename(sample = a, base_name = base_a) %>%
  bind_rows(
    circos_df %>%
      distinct(b, base_b) %>%
      rename(sample = b, base_name = base_b)
  ) %>%
  distinct()

# Create a mapping from sample to its village location (from both columns)
location_info <- bind_rows(
  circos_df %>% select(sample = a, location = Location_a),
  circos_df %>% select(sample = b, location = Location_b)
) %>% distinct()

# Merge sample_to_base with location info to create sample metadata
sample_meta <- sample_to_base %>%
  left_join(location_info, by = "sample")

# Create treatment information from both sample columns
treatment_info <- bind_rows(
  circos_df %>% select(sample = a, treatment = Treatment_timepoint_a),
  circos_df %>% select(sample = b, treatment = Treatment_timepoint_b)
) %>% distinct() %>%
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# Merge treatment info into sample metadata and set sample-level treatment as a factor
sample_meta <- sample_meta %>%
  left_join(treatment_info, by = "sample") %>%
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# --- Compute Overall Treatment Status per Individual ---
# For each individual (base_name), if all samples have the same treatment then use that;
# The overall factor levels are set as:
# Post, Pre, Unknown
individual_status <- sample_meta %>%
  mutate(treatment = as.character(treatment)) %>%
  group_by(base_name) %>%
  summarise(overall_treatment = if_else(n_distinct(treatment) == 1,
                                        treatment[1],
                                        "Unknown"),
            .groups = "drop") %>%
  mutate(overall_treatment = factor(overall_treatment,
                                    levels = c("Post", "Pre", "Unknown")))

# Merge the overall treatment status into sample metadata
sample_meta <- sample_meta %>%
  left_join(individual_status, by = "base_name")

# Order individuals by village, then by overall treatment status, then by individual (base_name) and sample
sample_meta <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  mutate(individual = factor(base_name, levels = unique(base_name)))

# Extract the ordered sample names (sectors for the circos plot)
ordered_samples <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  pull(sample)

### 2. Initialize Circos and Add Inner Tracks ###
circos.clear()
circos.par(start.degree = 90, gap.degree = 0.5, 
           track.margin = c(0.01, 0.01),
           cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = ordered_samples, xlim = c(0, 1))

# Add a track for Treatment Status (inner circle)
circos.track(ylim = c(0, 1),
             panel.fun = function(x, y) {
               s <- get.cell.meta.data("sector.index")
               # Look up treatment for the sample
               t <- treatment_info$treatment[treatment_info$sample == s]
               circos.rect(0, 0, 1, 1, col = treatment_colors[t], border = NA)
             },
             bg.border = NA, track.height = 0.05)

### 3. Draw the Relatedness Links ###
for (i in seq_len(nrow(circos_df))) {
  link_type <- ifelse(circos_df$Donor[i] == "Same donor", 1, 3)
  circos.link(
    circos_df$a[i], 0.5,
    circos_df$b[i], 0.5,
    col = relatedness_colors[circos_df$category[i]],
    lty = link_type,
    border = NA
  )
}

### 4. Add Outer Layer Spanning Samples from the Same Individual ###
# # Compute group boundaries: for each individual (base_name), find the first and last sample
# group_boundaries <- sample_to_base %>%
#   arrange(base_name, sample) %>%
#   group_by(base_name) %>%
#   summarise(
#     first_sample = sample[1],
#     last_sample  = sample[length(sample)],
#     .groups = "drop"
#   )

group_boundaries <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  group_by(base_name) %>%   # or group_by(base_name, overall_treatment) if you want to separate different treatment statuses for the same donor
  summarise(
    first_sample = dplyr::first(sample),
    last_sample  = dplyr::last(sample),
    .groups = "drop"
  )

# Add an extra (empty) track for drawing the outer group arcs
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {})

# Get the current track index for the outer annotation
current_track <- get.current.track.index()

# For each group (individual), draw an arc spanning from the first to the last sample,
# and add ticks at both boundaries facing inward.
for(i in seq_len(nrow(group_boundaries))) {
  first_sector <- group_boundaries$first_sample[i]
  last_sector  <- group_boundaries$last_sample[i]
  
  # Get angular boundaries (in degrees) for the first and last sectors
  start_angle <- get.cell.meta.data("cell.start.degree", sector.index = first_sector, track.index = current_track)
  end_angle   <- get.cell.meta.data("cell.end.degree", sector.index = last_sector, track.index = current_track)
  
  # Set the radius for the outer arc
  r <- 1.05
  
  # Create a sequence of angles (in radians) for the arc
  theta <- seq(start_angle, end_angle, length.out = 100) * pi/180
  x_arc <- r * cos(theta)
  y_arc <- r * sin(theta)
  
  # Draw the arc spanning the samples of this individual
  lines(x_arc, y_arc, col = "black", lwd = 2)
  
  # Define tick length (adjust as needed)
  tick_length <- 0.05
  
  # Convert boundary angles to radians
  theta_start <- start_angle * pi / 180
  theta_end   <- end_angle * pi / 180
  
  # Draw tick at the beginning of the span (tick drawn inward)
  lines(
    x = c(r * cos(theta_start), (r - tick_length) * cos(theta_start)),
    y = c(r * sin(theta_start), (r - tick_length) * sin(theta_start)),
    col = "black", lwd = 2
  )
  
  # Draw tick at the end of the span (tick drawn inward)
  lines(
    x = c(r * cos(theta_end), (r - tick_length) * cos(theta_end)),
    y = c(r * sin(theta_end), (r - tick_length) * sin(theta_end)),
    col = "black", lwd = 2
  )
  
  # (Optional) To add labels for the individual, uncomment and adjust:
  # mid_angle <- (start_angle + end_angle) / 2 * pi/180
  # x_text <- (r - tick_length - 0.02) * cos(mid_angle)
  # y_text <- (r - tick_length - 0.02) * sin(mid_angle)
  # text(x_text, y_text, labels = group_boundaries$base_name[i],
  #      srt = (start_angle + end_angle) / 2, cex = 0.8, adj = c(0.5, 0.5))
}

### 5. Finalize the Plot and Add Legends ###
# Clear the circos settings (optional, does not remove drawn graphics)
circos.clear()

# Create legends
lgd_treatment <- Legend(
  title = "Treatment Status",
  at = names(treatment_colors),
  legend_gp = gpar(fill = treatment_colors)
)
lgd_relatedness <- Legend(
  title = "Relatedness (Links)",
  at = names(relatedness_colors),
  legend_gp = gpar(fill = relatedness_colors)
)
lgd_links <- Legend(
  title = "Donor Relationship (Line Style)",
  at = c("Same Donor", "Different Donors"),
  type = "lines",
  legend_gp = gpar(col = "black", lty = c(1, 3), lwd = 2)
)

all_legends <- packLegend(lgd_treatment, lgd_relatedness, lgd_links,
                          direction = "vertical")
draw(all_legends, x = unit(1, "npc") - unit(1, "cm"), just = "right")

## network
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)

# Create edges from circos_df
edges <- circos_df %>%
  select(from = a, to = b, category, Donor, Treatment)

# Create nodes from unique samples
nodes <- data.frame(name = unique(c(edges$from, edges$to)))

# Create node-level attributes: treatment and village
node_info <- bind_rows(
  circos_df %>% select(sample = a, treatment = Treatment_timepoint_a, village = Location_a),
  circos_df %>% select(sample = b, treatment = Treatment_timepoint_b, village = Location_b)
) %>% distinct() %>%
  # Change NA treatment values to "Unknown"
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# Build the graph and join node attributes
g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  left_join(node_info, by = c("name" = "sample"))

# Define color palettes for edges, treatment, and village
edge_colors <- c("First Degree" = "#974181", "Second Degree" = "#d57946", "Third Degree" = "steelblue")
node_treatment_colors <- c("Pre" = "red", "Post" = "green", "Unknown" = "gray")
node_village_colors   <- c("Buhirigi" = "purple", "Kaiso" = "blue")

set.seed(42)
ggraph(g, layout = "fr") +
  # Edges: map category to color and Donor to line type
  geom_edge_link(aes(color = category, linetype = Donor), width = 1, show.legend = TRUE) +
  scale_edge_color_manual(values = edge_colors) +
  scale_edge_linetype_manual(values = c("Same donor" = "solid", "Different donors" = "dotted")) +
  # Nodes: fill by treatment and border color by village
  geom_node_point(aes(fill = treatment, color = village), size = 3, shape = 21,stroke=2) +
  scale_fill_manual(values = node_treatment_colors) +
  scale_color_manual(values = node_village_colors) +
  theme_graph(base_size = 14) +
  guides(
    edge_color = guide_legend(title = "Relatedness"),
    edge_linetype = guide_legend(title = "Donor"),
    fill = guide_legend(title = "Treatment"),
    color = guide_legend(title = "Village")
  ) +
  ggtitle("Network Graph of Kaiso Samples")

table(circos_df[,c('Donor','Treatment','category')])

### circos - Buhirigi ###
# Load required libraries
library(dplyr)
library(circlize)
library(ComplexHeatmap)

### 1. Data Preparation ###
# Subset and prepare the data (adjust column indices/names as needed)
circos_df <- ngsRelate_b[, c("a","b","Treatment_timepoint_a","Treatment_timepoint_b","Location_a","Location_b","Donor","category","base_a","base_b")]
circos_df <- circos_df %>%
  mutate(Treatment = if_else(Treatment_timepoint_a == Treatment_timepoint_b,
                             "Same", 
                             "Different"))

# Keep only rows with non-missing location (either Location_a or Location_b)
circos_df <- circos_df[ ( !is.na(circos_df$Location_a) | !is.na(circos_df$Location_b) ), ]
# Remove rows with category "Other"
circos_df <- circos_df %>% filter(category != "Other")

# Compute the base names for each sample (from columns 'a' and 'b')
circos_df$base_a <- sub("[-_][^-_]*$", "", circos_df$a)
circos_df$base_a <- sub("[_][^_]*$", "", circos_df$base_a)
circos_df$base_b <- sub("[-_][^-_]*$", "", circos_df$b)
circos_df$base_b <- sub("[_][^_]*$", "", circos_df$base_b)
# (Filter again for category != "Other" if needed)
circos_df <- circos_df %>% filter(category != "Other")

# Define color mappings
treatment_colors   <- c("Pre" = "red", "Post" = "green", "Unknown" = "gray")
relatedness_colors <- c("First Degree"  = "#974181", 
                        "Second Degree" = "#d57946",
                        "Third Degree"  = "steelblue")

# Create a mapping from sample to its base name (here the base name represents the individual)
sample_to_base <- circos_df %>%
  distinct(a, base_a) %>%
  rename(sample = a, base_name = base_a) %>%
  bind_rows(
    circos_df %>%
      distinct(b, base_b) %>%
      rename(sample = b, base_name = base_b)
  ) %>%
  distinct()

# Create a mapping from sample to its village location (from both columns)
location_info <- bind_rows(
  circos_df %>% select(sample = a, location = Location_a),
  circos_df %>% select(sample = b, location = Location_b)
) %>% distinct()

# Merge sample_to_base with location info to create sample metadata
sample_meta <- sample_to_base %>%
  left_join(location_info, by = "sample")

# Create treatment information from both sample columns
treatment_info <- bind_rows(
  circos_df %>% select(sample = a, treatment = Treatment_timepoint_a),
  circos_df %>% select(sample = b, treatment = Treatment_timepoint_b)
) %>% distinct() %>%
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# Merge treatment info into sample metadata and set sample-level treatment as a factor
sample_meta <- sample_meta %>%
  left_join(treatment_info, by = "sample") %>%
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# --- Compute Overall Treatment Status per Individual ---
# For each individual (base_name), if all samples have the same treatment then use that;
# otherwise, assign "Mixed". The overall factor levels are set as:
# Post, Pre, Unknown, then Mixed.
individual_status <- sample_meta %>%
  mutate(treatment = as.character(treatment)) %>%
  group_by(base_name) %>%
  summarise(overall_treatment = if_else(n_distinct(treatment) == 1,
                                        treatment[1],
                                        "Unknown"),
            .groups = "drop") %>%
  mutate(overall_treatment = factor(overall_treatment,
                                    levels = c("Post", "Pre", "Unknown")))

# Merge the overall treatment status into sample metadata
sample_meta <- sample_meta %>%
  left_join(individual_status, by = "base_name")

# Order individuals by village, then by overall treatment status, then by individual (base_name) and sample
sample_meta <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  mutate(individual = factor(base_name, levels = unique(base_name)))

# Extract the ordered sample names (sectors for the circos plot)
ordered_samples <- sample_meta %>%
  arrange(location, overall_treatment, base_name, sample) %>%
  pull(sample)

### 2. Initialize Circos and Add Inner Tracks ###
circos.clear()
circos.par(start.degree = 90, gap.degree = 0.5, 
           track.margin = c(0.01, 0.01),
           cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = ordered_samples, xlim = c(0, 1))

# Add a track for Treatment Status (inner circle)
circos.track(ylim = c(0, 1),
             panel.fun = function(x, y) {
               s <- get.cell.meta.data("sector.index")
               # Look up treatment for the sample
               t <- treatment_info$treatment[treatment_info$sample == s]
               circos.rect(0, 0, 1, 1, col = treatment_colors[t], border = NA)
             },
             bg.border = NA, track.height = 0.05)

### 3. Draw the Relatedness Links ###
for (i in seq_len(nrow(circos_df))) {
  link_type <- ifelse(circos_df$Donor[i] == "Same donor", 1, 3)
  circos.link(
    circos_df$a[i], 0.5,
    circos_df$b[i], 0.5,
    col = relatedness_colors[circos_df$category[i]],
    lty = link_type,
    border = NA
  )
}

### 4. Add Outer Layer Spanning Samples from the Same Individual ###
# Compute group boundaries: for each individual (base_name), find the first and last sample
group_boundaries <- sample_to_base %>%
  arrange(base_name, sample) %>%
  group_by(base_name) %>%
  summarise(
    first_sample = sample[1],
    last_sample  = sample[length(sample)],
    .groups = "drop"
  )

# Add an extra (empty) track for drawing the outer group arcs
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {})

# Get the current track index for the outer annotation
current_track <- get.current.track.index()

# For each group (individual), draw an arc spanning from the first to the last sample,
# and add ticks at both boundaries facing inward.
for(i in seq_len(nrow(group_boundaries))) {
  first_sector <- group_boundaries$first_sample[i]
  last_sector  <- group_boundaries$last_sample[i]
  
  # Get angular boundaries (in degrees) for the first and last sectors
  start_angle <- get.cell.meta.data("cell.start.degree", sector.index = first_sector, track.index = current_track)
  end_angle   <- get.cell.meta.data("cell.end.degree", sector.index = last_sector, track.index = current_track)
  
  # Set the radius for the outer arc
  r <- 1.05
  
  # Create a sequence of angles (in radians) for the arc
  theta <- seq(start_angle, end_angle, length.out = 100) * pi/180
  x_arc <- r * cos(theta)
  y_arc <- r * sin(theta)
  
  # Draw the arc spanning the samples of this individual
  lines(x_arc, y_arc, col = "black", lwd = 2)
  
  # Define tick length (adjust as needed)
  tick_length <- 0.05
  
  # Convert boundary angles to radians
  theta_start <- start_angle * pi / 180
  theta_end   <- end_angle * pi / 180
  
  # Draw tick at the beginning of the span (tick drawn inward)
  lines(
    x = c(r * cos(theta_start), (r - tick_length) * cos(theta_start)),
    y = c(r * sin(theta_start), (r - tick_length) * sin(theta_start)),
    col = "black", lwd = 2
  )
  
  # Draw tick at the end of the span (tick drawn inward)
  lines(
    x = c(r * cos(theta_end), (r - tick_length) * cos(theta_end)),
    y = c(r * sin(theta_end), (r - tick_length) * sin(theta_end)),
    col = "black", lwd = 2
  )
  
  # (Optional) To add labels for the individual, uncomment and adjust:
  # mid_angle <- (start_angle + end_angle) / 2 * pi/180
  # x_text <- (r - tick_length - 0.02) * cos(mid_angle)
  # y_text <- (r - tick_length - 0.02) * sin(mid_angle)
  # text(x_text, y_text, labels = group_boundaries$base_name[i],
  #      srt = (start_angle + end_angle) / 2, cex = 0.8, adj = c(0.5, 0.5))
}

### 5. Finalize the Plot and Add Legends ###
# Clear the circos settings (optional, does not remove drawn graphics)
circos.clear()

# Create legends
lgd_treatment <- Legend(
  title = "Treatment Status",
  at = names(treatment_colors),
  legend_gp = gpar(fill = treatment_colors)
)
lgd_relatedness <- Legend(
  title = "Relatedness (Links)",
  at = names(relatedness_colors),
  legend_gp = gpar(fill = relatedness_colors)
)
lgd_links <- Legend(
  title = "Donor Relationship (Line Style)",
  at = c("Same Donor", "Different Donors"),
  type = "lines",
  legend_gp = gpar(col = "black", lty = c(1, 3), lwd = 2)
)

all_legends <- packLegend(lgd_treatment, lgd_relatedness, lgd_links,
                          direction = "vertical")
draw(all_legends, x = unit(1, "npc") - unit(1, "cm"), just = "right")

## network
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)

# Create edges from circos_df
edges <- circos_df %>%
  select(from = a, to = b, category, Donor, Treatment)

# Create nodes from unique samples
nodes <- data.frame(name = unique(c(edges$from, edges$to)))

# Create node-level attributes: treatment and village
node_info <- bind_rows(
  circos_df %>% select(sample = a, treatment = Treatment_timepoint_a, village = Location_a),
  circos_df %>% select(sample = b, treatment = Treatment_timepoint_b, village = Location_b)
) %>% distinct() %>%
  # Change NA treatment values to "Unknown"
  mutate(treatment = if_else(is.na(treatment), "Unknown", treatment))

# Build the graph and join node attributes
g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  left_join(node_info, by = c("name" = "sample"))

# Define color palettes for edges, treatment, and village
edge_colors <- c("First Degree" = "#974181", "Second Degree" = "#d57946", "Third Degree" = "steelblue")
node_treatment_colors <- c("Pre" = "red", "Post" = "green", "Unknown" = "gray")
node_village_colors   <- c("Buhirigi" = "purple", "Kaiso" = "blue")

set.seed(42)
ggraph(g, layout = "fr") +
  # Edges: map category to color and Donor to line type
  geom_edge_link(aes(color = category, linetype = Donor), width = 1, show.legend = TRUE) +
  scale_edge_color_manual(values = edge_colors) +
  scale_edge_linetype_manual(values = c("Same donor" = "solid", "Different donors" = "dotted")) +
  # Nodes: fill by treatment and border color by village
  geom_node_point(aes(fill = treatment, color = village), size = 3, shape = 21,stroke=2) +
  scale_fill_manual(values = node_treatment_colors) +
  scale_color_manual(values = node_village_colors) +
  theme_graph(base_size = 14) +
  guides(
    edge_color = guide_legend(title = "Relatedness"),
    edge_linetype = guide_legend(title = "Donor"),
    fill = guide_legend(title = "Treatment"),
    color = guide_legend(title = "Village")
  ) +
  ggtitle("Network Graph of Buhirigi Samples")

table(circos_df[,c('Donor','Treatment','category')])

#### tree ####
# Load required libraries
library(ggtree)
library(ggplot2)
library(dplyr)
library(ape)

# Read the tree file
tree <- read.tree("phylo.treefile")

# Filter metadata for samples with non-missing Location
meta_location <- meta %>%
  filter(!is.na(Location))

# Prune the tree to include only tips with non-missing Location
tree_location <- drop.tip(tree, setdiff(tree$tip.label, meta_location$vcf_name))

# Match pruned tree tips with metadata for Location
meta_tree_location <- tibble(label = tree_location$tip.label) %>%
  left_join(meta_location, by = c("label" = "vcf_name"))

# Plot the circular tree colored by Location
p_location <- ggtree(tree_location, layout = "circular") %<+% meta_tree_location +
  geom_tippoint(aes(color = Location), size = 2, na.rm = TRUE) +
  ggtitle("Phylogenetic Tree Colored by Location") +
  theme(legend.position = "right")

# Filter metadata for samples with non-missing Treatment timepoint
meta_treatment <- meta %>%
  filter(!is.na(`Treatment timepoint`))

# Prune the tree to include only tips with non-missing Treatment timepoint
tree_treatment <- drop.tip(tree, setdiff(tree$tip.label, meta_treatment$vcf_name))

# Match pruned tree tips with metadata for Treatment timepoint
meta_tree_treatment <- tibble(label = tree_treatment$tip.label) %>%
  left_join(meta_treatment, by = c("label" = "vcf_name"))

# Plot the circular tree colored by Treatment timepoint with NA handling
p_treatment <- ggtree(tree_treatment, layout = "circular") %<+% meta_tree_treatment +
  geom_tippoint(aes(color = `Treatment timepoint`), size = 2, na.rm = TRUE) +
  scale_color_manual(
    values = c("Pre" = "orange", "Post" = "steelblue", "NA" = "gray"),
    na.translate = TRUE
  ) +
  ggtitle("Phylogenetic Tree Colored by Treatment Timepoint") +
  theme(legend.position = "right")

# Display the trees
print(p_location)
print(p_treatment)



#### Admixture
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)  # For arranging plots side by side

# Function for generating ADMIXTURE plots
generate_admixture_plot <- function(k_file, fam_file, meta = NULL, group_var, k_value) {
  # Read ADMIXTURE .Q file
  q_data <- read.table(k_file, header = FALSE)
  
  # Read .fam file
  fam_data <- read.table(fam_file, header = FALSE, 
                         col.names = c("FID", "IID", "PID", "MID", "Sex", "Phenotype"))
  
  # Add sample IDs
  q_data$Sample <- fam_data$IID
  
  if (group_var == "Location") {
    # Assign location based on sample ID prefix, then order Kaiso before Buhirigi
    q_data <- q_data %>%
      mutate(Location = case_when(
        grepl("^S1", Sample) ~ "Buhirigi",
        grepl("^S2", Sample) ~ "Kaiso",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(Location)) %>%
      mutate(Location = factor(Location, levels = c("Kaiso", "Buhirigi")))
    
    group_col <- "Location"
    
  } else if (group_var == "Treatment") {
    # Merge treatment information from metadata
    q_data <- q_data %>%
      left_join(meta %>% dplyr::select(vcf_name, `Treatment timepoint`), 
                by = c("Sample" = "vcf_name")) %>%
      filter(!is.na(`Treatment timepoint`)) %>%
      # Normalize 'Treatment timepoint' to ensure consistent matching
      mutate(`Treatment timepoint` = tolower(trimws(`Treatment timepoint`))) %>%
      mutate(`Treatment timepoint` = case_when(
        .data$`Treatment timepoint` == "pre" ~ "pre",
        .data$`Treatment timepoint` == "post" ~ "post",
        TRUE ~ .data$`Treatment timepoint`  # fallback for any unexpected value
      )) %>%
      # Force the order of 'pre' to come before 'post'
      mutate(`Treatment timepoint` = factor(`Treatment timepoint`, 
                                            levels = c("pre", "post")))
    
    group_col <- "Treatment timepoint"
  }
  
  # --- Sorting step ---
  # Sort samples within each group based on the proportion in the first cluster (V1).
  # You can change "V1" to another column if you wish to sort by a different cluster.
  q_data <- q_data %>% 
    group_by(!!sym(group_col)) %>% 
    arrange(desc(V1)) %>% 
    mutate(Sample = factor(Sample, levels = Sample)) %>% 
    ungroup()
  
  # Convert data to long format
  q_long <- melt(q_data, 
                 id.vars = c("Sample", group_col), 
                 variable.name = "Cluster", 
                 value.name = "Proportion")
  
  # Rename clusters dynamically
  num_clusters <- length(unique(q_long$Cluster))
  
  # Shades of blue (you can change the palette as you like)
  shades <- colorRampPalette(c("lightblue", "blue"))(num_clusters)
  
  # Optionally reassign palette if in Treatment mode
  if (group_var == "Treatment") {
    shades <- colorRampPalette(c("lightblue", "blue"))(num_clusters)
  }
  
  # Factor cluster labels
  q_long$Cluster <- factor(
    q_long$Cluster,
    levels = paste0("V", seq_len(num_clusters)),
    labels = paste0("Cluster ", seq_len(num_clusters))
  )
  
  # Generate plot grouped by Location or Treatment
  p <- ggplot(q_long, aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_wrap(
      as.formula(paste0("~`", group_col, "`")), 
      scales = "free_x", 
      strip.position = "bottom"
    ) +
    scale_fill_manual(values = shades) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      plot.title.position = "plot",
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    ) +
    ggtitle(paste("K =", k_value)) +
    theme(plot.title = element_text(hjust = 1))
  
  return(p)
}


# Define file paths
fam_file <- "king_data.fam"  # .fam file
q_files <- list("king_data.2.Q", "king_data.3.Q", "king_data.4.Q")  # ADMIXTURE .Q files
k_values <- c(2, 3, 4)  # Values of K

# Generate plots for Location
plots_location <- lapply(seq_along(q_files), function(i) {
  generate_admixture_plot(
    k_file   = q_files[[i]],
    fam_file = fam_file,
    group_var = "Location",
    k_value = k_values[i]
  )
})

# Generate plots for Treatment
plots_treatment <- lapply(seq_along(q_files), function(i) {
  generate_admixture_plot(
    k_file   = q_files[[i]],
    fam_file = fam_file,
    meta     = meta,
    group_var = "Treatment",
    k_value = k_values[i]
  )
})

# Combine Location plots into a single vertical stack
combined_plot_location <- wrap_plots(plots_location, ncol = 1) &
  labs(x = "Village", y = "Ancestry Proportion") &
  theme(
    plot.margin = margin(10, 30, 30, 30),  # Adjust margin for alignment
    plot.tag.position = "right"
  )

# Combine Treatment plots into a single vertical stack
combined_plot_treatment <- wrap_plots(plots_treatment, ncol = 1) &
  labs(x = "Treatment status", y = "Ancestry Proportion") &
  theme(
    plot.margin = margin(10, 30, 30, 30),  # Adjust margin for alignment
    plot.tag.position = "right"
  )

# Display the Location-based grid
print(combined_plot_location)

# Display the Treatment-based grid
print(combined_plot_treatment)


## Shannon Diversity Index
# Shannon Diversity Index function
calculate_shannon_index <- function(data, group_var) {
  data %>%
    group_by(!!sym(group_var), Cluster) %>%
    summarize(mean_proportion = mean(Proportion, na.rm = TRUE)) %>%
    group_by(!!sym(group_var)) %>%
    summarize(
      shannon_index = -sum(mean_proportion * log(mean_proportion), na.rm = TRUE)
    )
}

# Extract q_long for Location
q_long_location <- do.call(rbind, lapply(1:length(q_files), function(i) {
  q_data <- read.table(q_files[[i]], header = FALSE)
  fam_data <- read.table(fam_file, header = FALSE, col.names = c("FID", "IID", "PID", "MID", "Sex", "Phenotype"))
  q_data$Sample <- fam_data$IID
  q_data <- q_data %>%
    mutate(Location = case_when(
      grepl("^S1", Sample) ~ "Buhirigi",
      grepl("^S2", Sample) ~ "Kaiso",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Location))
  melt(q_data, id.vars = c("Sample", "Location"), variable.name = "Cluster", value.name = "Proportion")
}))



# Shannon's Diversity and Equitability Index for Location
# Function to calculate Shannon metrics for a single q_data
calculate_shannon_all_k <- function(k_file, fam_file, group_var, meta = NULL) {
  q_data <- read.table(k_file, header = FALSE)
  fam_data <- read.table(fam_file, header = FALSE, col.names = c("FID", "IID", "PID", "MID", "Sex", "Phenotype"))
  q_data$Sample <- fam_data$IID
  
  if (group_var == "Location") {
    q_data <- q_data %>%
      mutate(Location = case_when(
        grepl("^S1", Sample) ~ "Buhirigi",
        grepl("^S2", Sample) ~ "Kaiso",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(Location))
    group_col <- "Location"
  } else if (group_var == "Treatment") {
    q_data <- q_data %>%
      left_join(meta %>% dplyr::select(vcf_name, `Treatment timepoint`), 
                by = c("Sample" = "vcf_name")) %>%
      filter(!is.na(`Treatment timepoint`))
    group_col <- "Treatment timepoint"
  }
  
  q_long <- melt(q_data, id.vars = c("Sample", group_col), 
                 variable.name = "Cluster", value.name = "Proportion")
  
  num_clusters <- length(unique(q_long$Cluster))
  
  q_long %>%
    group_by(!!sym(group_col), Cluster) %>%
    summarize(mean_proportion = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
    group_by(!!sym(group_col)) %>%
    summarize(
      shannon_index = -sum(mean_proportion * log(mean_proportion), na.rm = TRUE),
      shannon_equitability = (-sum(mean_proportion * log(mean_proportion), na.rm = TRUE)) / log(num_clusters)
    )
}

# Iterate over all K values for Location
shannon_location_all_k <- lapply(1:length(q_files), function(i) {
  calculate_shannon_all_k(q_files[[i]], fam_file, group_var = "Location")
}) %>%
  bind_rows(.id = "K") %>%
  mutate(K = as.numeric(K))

print(shannon_location_all_k)

# Iterate over all K values for Treatment
shannon_treatment_all_k <- lapply(1:length(q_files), function(i) {
  calculate_shannon_all_k(q_files[[i]], fam_file, group_var = "Treatment", meta = meta)
}) %>%
  bind_rows(.id = "K") %>%
  mutate(K = as.numeric(K))

print(shannon_treatment_all_k)


## King
# R example: scatter-plot kinship vs IBS0
library(data.table) 
library(ggplot2)
kin <- fread("king_kinship.kin0")               # tab-delimited by default
ggplot(kin, aes(Kinship, IBS0)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(0.177, 0.088, 0.044), linetype = "dashed") +
  labs(x = "Kinship coefficient", y = "Proportion IBS0",
       title = "KING pairwise relatedness: Kinship vs IBS0") +
  theme_bw()
