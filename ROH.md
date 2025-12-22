#### Calculate the ROH using Plink

###### Convert VCF to PLINK binary format (BED/BIM/FAM)
```bash
plink --vcf /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/ROH/Guam_rail_biallelic_snps.vcf.gz \
      --make-bed \
 --allow-extra-chr \
      --out /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/ROH/Guam_rail
```

###### Calculate the ROH using Plink
```bash
plink --bfile /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/ROH/Guam_rail \
      --homozyg \
      --homozyg-window-snp 50 \
      --homozyg-snp 50 \
      --homozyg-kb 500 \
      --homozyg-gap 1000 \
      --homozyg-density 50 \
      --homozyg-window-missing 5 \
      --homozyg-window-het 3 \
--allow-extra-chr \
      --out /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/ROH/Guam_rail/Guam_rail_ROH
```


#### Visualization of ROH 

##### ---------------------------------------
##### ROH Analysis Script for Dama Gazelle
##### Working Directory: F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/Plink_Final
##### ---------------------------------------
```bash
#######################################################
# ROH Analysis Script for Guam Rail
#######################################################

# -----------------------------
# 1. Load required packages
# -----------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges) 
library(scales) 

# -----------------------------
# 2. Set Parameters
# -----------------------------
# Update this path to your actual folder
setwd("F:/Collaborative_Projects/Guam_Rail/ROH/Plink_Final")

genome_size <- 1.16e9  # 1.16 Gb
samples <- c("FMNH390989", "N23-0063", "N23-0568")
output_dir <- getwd()

# Shared aesthetics
individual_colors <- c(
  "FMNH390989" = "forestgreen",
  "N23-0063"   = "tomato",
  "N23-0568"   = "maroon"
)

individual_shapes <- c(
  "FMNH390989" = 15,
  "N23-0063"   = 16,
  "N23-0568"   = 17
)

# Create a lookup table for samples
sample_info <- tibble(
  IID = samples,
  Species = "Guam Rail"
)

# -----------------------------
# 3. Load Data
# -----------------------------
roh_segments <- fread(file.path(output_dir, "Guam_rail_ROH.hom"))
roh_indiv    <- fread(file.path(output_dir, "Guam_rail_ROH.hom.indiv"))

# --------------------------------------------------------------------------------------
# PLOT 1: Percent Genome in ROH per Individual (Proper Boxplot via Scaffolds)
# --------------------------------------------------------------------------------------
# We use segments to calculate ROH % per scaffold so we have enough data for a "box"

scaffold_lengths <- roh_segments %>%
  group_by(CHR) %>%
  summarise(Scaff_End = max(POS2), .groups = "drop")

roh_dist_data <- roh_segments %>%
  filter(IID %in% samples) %>%
  group_by(IID, CHR) %>%
  summarise(Sum_KB = sum(KB), .groups = "drop") %>%
  left_join(scaffold_lengths, by = "CHR") %>%
  mutate(Percent_ROH_Scaff = (Sum_KB * 1000 / Scaff_End) * 100)

p1 <- ggplot(roh_dist_data, aes(x = IID, y = Percent_ROH_Scaff, fill = IID)) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.4, color = "black") +
  scale_fill_manual(values = individual_colors) +
  labs(
    x = "Individual",
    y = "% ROH per Scaffold"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    #Axis titles (labels)
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    
    #Axis tick values
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    
    legend.position = "none"
  )

print(p1)


# --------------------------------------------------------------------------------------
# PLOT 2: Percent Genome in ROH by Size Category (Barplot)
# --------------------------------------------------------------------------------------
roh_by_cat <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000    ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000   ~ "1-5Mb",
    KB >= 5000 & KB < 10000  ~ "5-10Mb",
    KB >= 10000              ~ ">10Mb"
  )) %>%
  filter(!is.na(ROH_category)) %>%
  group_by(IID, ROH_category) %>%
  summarise(
    Total_Mb = sum(KB / 1000),
    .groups = "drop"
  ) %>%
  mutate(
    Percent_Genome = (Total_Mb * 1e6 / genome_size) * 100
  )

p2 <- ggplot(roh_by_cat, aes(x = IID, y = Percent_Genome, fill = ROH_category)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c(
    "0.1-1Mb"  = "#e6ab02",
    "1-5Mb"    = "#d95f02",
    "5-10Mb"   = "#1b9e77",
    ">10Mb"    = "#7570b3"
  )) +
  labs(
    x = "Individual",
    y = "% Genome in ROH",
    fill = "Size Category"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    # Axis titles
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    # Axis tick values
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    
    # Legend inside TOP-LEFT
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold")
  )
p2

# --------------------------------------------------------------------------------------
# PLOT 3: Number of ROH vs Total ROH Length (Scatterplot)
# --------------------------------------------------------------------------------------
# 1. Prepare detailed data
roh_detailed_comp <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000   ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000  ~ "1-5Mb",
    KB >= 5000              ~ ">5Mb"
  )) %>%
  filter(!is.na(ROH_category)) %>%
  group_by(IID, ROH_category) %>%
  summarise(
    Num_ROH = n(),
    Sum_ROH_Mb = sum(KB / 1000),
    .groups = "drop"
  ) %>%
  mutate(
    ROH_category = factor(
      ROH_category,
      levels = c("0.1-1Mb", "1-5Mb", ">5Mb")
    )
  )

# 2. Create the Plot
p3_comparative <- ggplot(
  roh_detailed_comp,
  aes(x = Sum_ROH_Mb, y = Num_ROH, color = IID, group = IID)
) +
  theme_bw() +
  
  # Connect categories
  geom_path(size = 1.2, alpha = 0.5, linetype = "dashed") +
  
  # Points
  geom_point(aes(shape = IID), size = 7, stroke = 1.5) +
  
  # Category labels
  geom_text(
    aes(label = ROH_category),
    hjust = -0.2, vjust = -1,
    size = 4.5, fontface = "bold", color = "black"
  ) +
  
  scale_color_manual(values = individual_colors) +
  scale_shape_manual(values = individual_shapes) +
  
  #Zoom BOTH axes to 0–110
  coord_cartesian(xlim = c(0, 110), ylim = c(0, 110)) +
  
  labs(
    x = "Total ROH Length (Mb)",
    y = "Number of ROH Segments"
  ) +
  
  theme(
    axis.text.x = element_text(size = 15, color = "black", face = "bold"),
    axis.text.y = element_text(size = 15, color = "black", face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
    
    # Legend inside TOP-LEFT
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(1.2, "cm"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

print(p3_comparative)

# Save the plot
ggsave(
  "Plot3_ROH_Comparative_Profile.jpeg",
  p3_comparative,
  width = 12,
  height = 8,
  dpi = 300
)

# --------------------------------------------------------------------------------------
# PLOT 4: Cumulative % Genome in ROH vs ROH Size (CDF)
# --------------------------------------------------------------------------------------
merge_overlaps_chr <- function(df) {
  df %>%
    group_by(CHR) %>%
    group_modify(~{
      ir <- IRanges(start = .x$POS1, end = .x$POS2)
      merged <- reduce(ir)
      tibble(ROH_Mb = width(merged) / 1e6)
    }) %>%
    ungroup()
}

roh_merged <- roh_segments %>%
  filter(IID %in% samples) %>%
  group_by(IID) %>%
  group_modify(~ merge_overlaps_chr(.x)) %>%
  ungroup()

roh_cdf <- roh_merged %>%
  arrange(IID, ROH_Mb) %>%
  group_by(IID) %>%
  mutate(
    Cumulative_bp = cumsum(ROH_Mb * 1e6),
    Percent_Cumulative = (Cumulative_bp / genome_size) * 100
  )

p4 <- ggplot(roh_cdf, aes(x = ROH_Mb, y = Percent_Cumulative, color = IID)) +
  geom_line(linewidth = 1.5) +
  scale_x_log10(labels = label_comma()) +
  scale_color_manual(values = individual_colors) +
  labs(
    x = "ROH Size (Mb) (log10)",
    y = "Cumulative % Genome in ROH"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    # Axis titles
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    
    #Axis tick values
    axis.text.x  = element_text(size = 15, color = "black"),
    axis.text.y  = element_text(size = 15, color = "black"),
    
    #Legend text
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    
    # Keep legend inside (as you had)
    legend.position = c(0.15, 0.85)
  )

print(p4)

# -----------------------------
# 5. Save all Plots
# -----------------------------
ggsave("Plot1_ROH_Boxplot.jpeg", p1, width=6, height=5, dpi=300)
ggsave("Plot2_ROH_Barplot.jpeg", p2, width=8, height=5, dpi=300)
ggsave("Plot3_ROH_Scatter.jpeg", p3, width=7, height=5, dpi=300)
ggsave("Plot4_ROH_CDF.jpeg", p4, width=8, height=6, dpi=300)
```

