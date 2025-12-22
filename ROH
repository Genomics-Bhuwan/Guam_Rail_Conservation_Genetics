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
# Loading the packages
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges) # for merging overlapping ROHs
library(scales) # Needed for better log-scale breaks

# -----------------------------
# Set Working Directory & Parameters
# -----------------------------
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/Plink_Final")

# Genome size set to 3 billion base pairs (3e9 bp)
genome_size <- 3e9
samples <- c("SRR17134085","SRR17134086","SRR17129394", # Addra
             "SRR17134087","SRR17134088") # Mhorr
output_dir <- getwd()

# Define shared aesthetics
subspecies_colors <- c("Addra"="skyblue", "Mhorr"="orange")
individual_shapes <- c(
  "SRR17134085"=15, # Solid Square
  "SRR17134086"=16, # Solid Circle
  "SRR17129394"=17, # Solid Triangle
  "SRR17134087"=18, # Solid Diamond
  "SRR17134088"=8  # Star
)

# --- Load Segment Data (Used for Plots 3, 4, 5) ---
# ASSUMPTION: This file is the corrected PLINK output with ROH >= 10 KB
# NOTE: If your 10kb data is in a different file, update the name here!
roh_segments <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom"))

# --------------------------------------------------------------------------------------
# 1. Percent Genome in ROH per Subspecies (Boxplot)
# --------------------------------------------------------------------------------------
# Uses .hom.indiv (Assumed corrected to include 10kb ROHs)
roh_indiv <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom.indiv"))

roh_sub <- roh_indiv %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"),
                          "Addra", "Mhorr"),
         Percent_ROH = (KB * 1000 / genome_size) * 100)

p1 <- ggplot(roh_sub, aes(x=Species, y=Percent_ROH, fill=Species)) +
  geom_boxplot() +
  geom_jitter(width=0.1, size=3) +
  labs(x="Sub-species",
       y="% Genome in ROH",
       fill="Sub-species") +
  scale_fill_manual(values=subspecies_colors) +
  scale_y_continuous(breaks = c(0, 4, 8, 12, 16), limits = c(0, 16)) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill="white", colour="black", size=0.6),
    axis.text.y = element_text(color="black", size=13),
    axis.text.x = element_text(color="black", size=13),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )

print(p1)
ggsave(file.path(output_dir, "Percent_ROH_by_subspecies.jpeg"), p1, width=6, height=4, dpi=300)
ggsave(file.path(output_dir, "Percent_ROH_by_subspecies.pdf"), p1, width=6, height=4)

# --------------------------------------------------------------------------------------
# 2.  Percent Genome in ROH by Size Category (Barplot) - FILTER FIXED
# --------------------------------------------------------------------------------------
# Use the loaded roh_segments data
roh_segments_3 <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"),
                          "Addra", "Mhorr"),
         ROH_Mb = KB / 1000)

roh_segments_3$CHR <- factor(roh_segments_3$CHR, levels=unique(roh_segments_3$CHR))

roh_segments_filtered <- roh_segments_3 %>%
  # --- CRITICAL FIX: Removed the line: filter(KB * 1000 > 100000) ---
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000 ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000 ~ "1-5Mb",
    KB >= 5000 & KB < 10000 ~ "5-10Mb"
  )) %>%
  filter(!is.na(ROH_category)) # Filter out ROHs < 100 KB or > 10 Mb


roh_by_cat <- roh_segments_filtered %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(Total_ROH_Mb = sum(KB / 1000), .groups="drop") %>%
  mutate(Percent_Genome_ROH = (Total_ROH_Mb * 1e6 / genome_size) * 100)

p3 <- ggplot(roh_by_cat, aes(x=IID, y=Percent_Genome_ROH, fill=ROH_category)) +
  geom_bar(stat="identity", color="black") +
  facet_wrap(~Species, scales="free_x") +
  labs(x="Individual",
       y="% Genome in ROH",
       fill="ROH size category") +
  scale_fill_manual(values=c("0.1-1Mb"="#e6ab02", "1-5Mb"="#d95f02", "5-10Mb"="#1b9e77")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.text.y = element_text(color="black", size=12),
    axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill="white", colour="black", size=0.6),
    strip.text = element_text(color="black", size=14, face="bold")
  )

print(p3)
ggsave(file.path(output_dir, "Percent_ROH_by_size_category.jpeg"), p3, width=8, height=5, dpi=300)
ggsave(file.path(output_dir, "Percent_ROH_by_size_category.pdf"), p3, width=8, height=5)


# --------------------------------------------------------------------------------------
# 3. Number of ROH vs Total ROH (Mb) per Size Category (Scatterplot) - FILTER FIXED
# --------------------------------------------------------------------------------------
# Use the loaded roh_segments data

roh_threecat <- roh_segments %>%
  # --- CRITICAL FIX: Removed the line: filter(KB * 1000 > 100000) ---
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000 ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000 ~ "1-5Mb",
    KB >= 5000 & KB < 10000 ~ "5-10Mb",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(ROH_category))


roh_count_sum <- roh_threecat %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(Num_ROH = n(),
            Sum_ROH_Mb = sum(ROH_Mb),
            .groups="drop")

p5 <- ggplot(roh_count_sum, aes(x=Sum_ROH_Mb, y=Num_ROH)) +
  geom_point(aes(shape=IID, color=Species), size=4, alpha=0.9) +
  facet_wrap(~ROH_category) +
  labs(title="Number of ROH vs Total ROH Length (Mb) per Size Category",
       x="Total ROH Length (Mb)",
       y="Number of ROH",
       shape="Individual Sample",
       color="Sub-species") +
  scale_color_manual(values=subspecies_colors) +
  scale_shape_manual(values=individual_shapes) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill="white", colour="black", size=0.6)
  )

print(p5)
ggsave(file.path(output_dir, "ROH_Count_vs_Sum_shapes.jpeg"), p5, width=9, height=5, dpi=300)
ggsave(file.path(output_dir, "ROH_Count_vs_Sum_shapes.pdf"), p5, width=9, height=5)



  # ----------------------------------------------------
# 4.  Cumulative % Genome in ROH vs ROH Size (PLAIN LINE CDF) - CORRECTED
# ----------------------------------------------------
# ----------------------------------------------------
# ROH Cumulative % Genome vs ROH Size (PLAIN LINE CDF)
# ----------------------------------------------------
library(IRanges)
library(dplyr)
library(ggplot2)
library(scales)

# -----------------------------
# 1. Merge ROHs per chromosome
# -----------------------------
merge_overlaps_chr <- function(df) {
  df %>%
    group_by(CHR) %>%
    group_modify(~{
      ir <- IRanges(start = .x$POS1, end = .x$POS2)
      merged <- reduce(ir)
      tibble(ROH_Mb = width(merged) / 1e6)
    }) %>% ungroup()
}

# -----------------------------
# 2. Sample and species info
# -----------------------------
sample_info <- tibble(
  IID = c("SRR17129394","SRR17134085","SRR17134086","SRR17134087","SRR17134088"),
  Species = c("Addra","Addra","Addra","Mhorr","Mhorr")
)

# Optional: assign colors per individual (blue shades = Addra, orange shades = Mhorr)
addra_colors <- c("#1f78b4", "#6baed6", "#a6cee3")  # 3 Addra
mhorr_colors <- c("#ff7f00", "#fdbf6f")            # 2 Mhorr

individual_colors <- c(
  setNames(addra_colors, sample_info$IID[sample_info$Species=="Addra"]),
  setNames(mhorr_colors, sample_info$IID[sample_info$Species=="Mhorr"])
)

# -----------------------------
# 3. Merge ROHs genome-wide and attach species
# -----------------------------
roh_merged <- roh_segments %>%
  filter(IID %in% sample_info$IID) %>%
  group_by(IID) %>%
  group_modify(~ merge_overlaps_chr(.x)) %>%
  ungroup() %>%
  left_join(sample_info, by="IID")

# -----------------------------
# 4. Compute CDF
# -----------------------------
roh_cdf <- roh_merged %>%
  arrange(IID, ROH_Mb) %>%
  group_by(IID, Species) %>%
  mutate(
    Total_ROH_bp = sum(ROH_Mb * 1e6),
    Cumulative_bp = cumsum(ROH_Mb * 1e6),
    Percent_Cumulative = (Cumulative_bp / genome_size) * 100,
    Label = paste(IID, "(", Species, ")", sep="")   # For legend
  ) %>%
  ungroup()

# -----------------------------
# 5. Plot CDF with legend slightly left
# -----------------------------
p4 <- ggplot(roh_cdf, aes(x = ROH_Mb, y = Percent_Cumulative,
                          group = IID, color = IID)) +
  geom_line(linewidth = 1.2) +          # lines only
  scale_x_log10(labels = label_comma()) +
  scale_color_manual(values = individual_colors,
                     labels = setNames(roh_cdf$Label[match(individual_colors %>% names(), roh_cdf$IID)],
                                       names(individual_colors)),
                     name = "Individual (Sub-species)") +
  labs(
    x = "ROH Size (Mb) (log10 scale)",
    y = "Cumulative % Genome in ROH"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = c(0.22, 0.88),      # moved slightly up and left
    legend.background = element_rect(fill = "white", color = "black", size = 0.6),
    legend.box = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

# -----------------------------
# 6. Print and save
# -----------------------------
print(p4)

ggsave("ROH_CDF_plot_species_colors_legend_adjusted_left.jpeg", plot = p4, width = 12, height = 8, dpi = 300)
ggsave("ROH_CDF_plot_species_colors_legend_adjusted_left.pdf", plot = p4, width = 12, height = 8)
```

