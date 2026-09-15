## Author: Ibrahim Fagbohun
## Date: 04/04/26
## Aim: Analyse data from the survival experiment study.  

## Get and set working directory 
getwd()
setwd("/Users/ibrah/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/My Research/Data")

## Load relevant libraries 
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(vegan)
library(MASS)
library(broom)
library(purrr)
library(emmeans)
library(ggalluvial)
library(glmmTMB)
library(patchwork)

## Load relevant data 
Survival_Experiment_Data <- read.csv("https://raw.githubusercontent.com/Ibrahim-Fagbohun/Survival_Experiment/main/Survival_Experiment_Data.csv")
Functional_Trait_Data <- read_excel("~/OneDrive - The Pennsylvania State University/My Research/Data/CBT_Survival_Data/Cleaned_CBT_Survival_Data/Functional_Trait_Data.xlsx")

###############################################################################################################################
## Objective One (Effect of site type on community composition)
## Separate full data into metadata and taxa data 
Taxa_Data <- Full_Survival_Experiment_Data [, 6:74]
Metadata <- Full_Survival_Experiment_Data [, 1:5]

## For the NMDS and PERMANOVA analysis, separate full data into Fall and Spring Data to avoid strong seasonal variations from overwhelming the analysis. 
Fall_Samples <- Full_Survival_Experiment_Data %>%
  filter(Sampling_Time == "Fall")

Spring_Samples <- Full_Survival_Experiment_Data %>%
  filter(Sampling_Time == "Spring")

## Separate fall dataset into Community matrix and metadata 
Fall_Metadata <- Fall_Samples [, 1:5]
Fall_Community_Matrix <- Fall_Samples [, 6:74]

## Clean fall data with the following steps
## Remove all taxa that sum up to zero because that they were only present in the spring
Fall_Community_Matrix_Cleaned <- Fall_Community_Matrix[, colSums(Fall_Community_Matrix) != 0] ## 27 taxa removed. 

## Remove rare taxa, i.e taxa that are not present in at least 3 samples. 
Fall_Community_Matrix_Cleaned <- Fall_Community_Matrix_Cleaned[, colSums(Fall_Community_Matrix_Cleaned > 0) >= 3] ## 17 taxa removed. 

## Do a square root transformation to reduce the impact of highly dominant taxa 
Fall_Community_Matrix_Transformed <- sqrt(Fall_Community_Matrix_Cleaned)

## Run NMDS 
set.seed(123)
Full_Fall_NMDS <- metaMDS(
  Fall_Community_Matrix_Transformed,
  distance = "bray",
  k = 3,
  trymax = 1000
) ## Converged after 20 iterations. 

## Extract Stress Value 
Full_Fall_NMDS$stress ## 0.149

## Extract other scores 
Full_Fall_NMDS_Scores <- as.data.frame(scores(Full_Fall_NMDS, display = "sites"))

## Merge NMDS Scores with meta data 
Full_Fall_NMDS_Data <- cbind(Full_Fall_NMDS_Scores, Fall_Metadata)

## Compute centroids for site type 
Fall_Centroids <- aggregate(cbind(NMDS1, NMDS2, NMDS3) ~ Site_Type,
                            data = Full_Fall_NMDS_Data, FUN = mean)

## Run a PERMANOVA to formally test for differences between pre and post transplant community during the spring and fall seasons. 
## Build bray curtis distance matrix for PERMANOVA 
Fall_Bray_Distance <- vegdist(Fall_Community_Matrix_Transformed, method = "bray")

## Run PERMANOVA 
Fall_Permanova <- adonis2(
  Fall_Bray_Distance ~ Site_Type,
  data = Fall_Metadata,
  permutations = 999
)

Fall_Permanova ## F = 7.6626, p = 0.001, R2 = 0.1185

## Test for homogeneity dispersion to be sure that the PERMANOVA results are not inflated by unequal variations 
Fall_Perm_Disp <- betadisper(Fall_Bray_Distance, Fall_Metadata$Site_Type)
anova(Fall_Perm_Disp) ## F= 3.5959, p = 0.063

# Recode factor level before plotting NMDS 
Full_Fall_NMDS_Data$Site_Type <- fct_recode(
  Full_Fall_NMDS_Data$Site_Type,
  "Pre-transplant" = "Donor",
  "Post-transplant" = "Recipient"
)

## Recode factor levels for Centroids too 
Fall_Centroids$Site_Type <- recode(
  Fall_Centroids$Site_Type,
  "Donor" = "Pre-transplant",
  "Recipient" = "Post-transplant"
)

## Build the NMDS plot 
## Primary axis
Fall_NMDS_Main_Axis <- ggplot(
  Full_Fall_NMDS_Data,
  aes(x = NMDS1, y = NMDS2, color = Site_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Site_Type (your main inference groups)
  stat_ellipse(aes(group = Site_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Site_Type
  geom_point(
    data = Fall_Centroids,
    aes(x = NMDS1, y = NMDS2, color = Site_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Donor vs Recipient
  scale_color_manual(values = c("Pre-transplant" = "darkred",
                                "Post-transplant" = "black")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS2") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) + 
  theme(legend.position = "none")

## Secondary axis 
Fall_NMDS_Secondary_Axis <- ggplot(
  Full_Fall_NMDS_Data,
  aes(x = NMDS1, y = NMDS3, color = Site_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Site_Type (your main inference groups)
  stat_ellipse(aes(group = Site_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Site_Type
  geom_point(
    data = Fall_Centroids,
    aes(x = NMDS1, y = NMDS3, color = Site_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Donor vs Recipient
  scale_color_manual(values = c("Pre-transplant" = "darkred",
                                "Post-transplant" = "black")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS3") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  )

print(Fall_NMDS_Main_Axis)
print(Fall_NMDS_Secondary_Axis)

## Combine figures with patchwork 
library(patchwork)

Fall_Combined_Plot <- Fall_NMDS_Main_Axis + Fall_NMDS_Secondary_Axis + 
  plot_layout(ncol = 2)
print(Fall_Combined_Plot)

####### Spring 
## Separate spring dataset into Community matrix and metadata 
Spring_Metadata <- Spring_Samples [, 1:5]
Spring_Community_Matrix <- Spring_Samples [, 6:74]

## Clean spring data with the following steps
## Remove all taxa that sum up to zero because that they were only present in the fall 
Spring_Community_Matrix_Cleaned <- Spring_Community_Matrix[, colSums(Spring_Community_Matrix) != 0] ## 16 taxa removed. 

## Remove rare taxa, i.e taxa that are not present in at least 3 samples. 
Spring_Community_Matrix_Cleaned <- Spring_Community_Matrix_Cleaned[, colSums(Spring_Community_Matrix_Cleaned > 0) >= 3] ## 25 taxa removed. 

## Do a square root transformation to reduce the impact of highly dominant taxa 
Spring_Community_Matrix_Transformed <- sqrt(Spring_Community_Matrix_Cleaned)

## Run NMDS 
set.seed(123)
Full_Spring_NMDS <- metaMDS(
  Spring_Community_Matrix_Transformed,
  distance = "bray",
  k = 3,
  trymax = 1000
) ## Converged after 20 iterations 

## Extract Stress Value 
Full_Spring_NMDS$stress ## 0.163

## Extract other scores 
Full_Spring_NMDS_Scores <- as.data.frame(scores(Full_Spring_NMDS, display = "sites"))

## Merge NMDS Scores with meta data 
Full_Spring_NMDS_Data <- cbind(Full_Spring_NMDS_Scores, Spring_Metadata)

## Compute centroids for site type 
Spring_Centroids <- aggregate(cbind(NMDS1, NMDS2, NMDS3) ~ Site_Type,
                              data = Full_Spring_NMDS_Data, FUN = mean)

## Run a PERMANOVA to formally test for differences between pre and post transplant community during the spring season 
## Build matrix for PERMANOVA 
Spring_Bray_Distance <- vegdist(Spring_Community_Matrix_Transformed, method = "bray")

## Run PERMANOVA 
Spring_Permanova <- adonis2(
  Spring_Bray_Distance ~ Site_Type,
  data = Spring_Metadata,
  permutations = 999
)

Spring_Permanova ## F = 14.977, p = 0.001, R2 = 0.20523

## Test for homogeneity dispersion to be sure that the PERMANOVA results are not inflated by unequal variations 
Spring_Perm_Disp <- betadisper(Spring_Bray_Distance, Spring_Metadata$Site_Type)
anova(Spring_Perm_Disp) ## F= 1.4626, p = 0.2314

# Recode factor level before plotting NMDS 
Full_Spring_NMDS_Data$Site_Type <- fct_recode(
  Full_Spring_NMDS_Data$Site_Type,
  "Pre-transplant" = "Donor",
  "Post-transplant" = "Recipient"
)

## Recode factor levels for Centroids too 
Spring_Centroids$Site_Type <- recode(
  Spring_Centroids$Site_Type,
  "Donor" = "Pre-transplant",
  "Recipient" = "Post-transplant"
)

## Build the NMDS plot 
## Primary axis
Spring_NMDS_Main_Axis <- ggplot(
  Full_Spring_NMDS_Data,
  aes(x = NMDS1, y = NMDS2, color = Site_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Site_Type (your main inference groups)
  stat_ellipse(aes(group = Site_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Site_Type
  geom_point(
    data = Spring_Centroids,
    aes(x = NMDS1, y = NMDS2, color = Site_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Donor vs Recipient
  scale_color_manual(values = c("Pre-transplant" = "darkred",
                                "Post-transplant" = "black")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS2") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) + 
  theme(legend.position = "none")

## Secondary axis 
Spring_NMDS_Secondary_Axis <- ggplot(
  Full_Spring_NMDS_Data,
  aes(x = NMDS1, y = NMDS3, color = Site_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Site_Type (your main inference groups)
  stat_ellipse(aes(group = Site_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Site_Type
  geom_point(
    data = Spring_Centroids,
    aes(x = NMDS1, y = NMDS3, color = Site_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Donor vs Recipient
  scale_color_manual(values = c("Pre-transplant" = "darkred",
                                "Post-transplant" = "black")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS3") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  )

print(Spring_NMDS_Main_Axis)
print(Spring_NMDS_Secondary_Axis)

## Combine figures with patchwork 
library(patchwork)

Spring_Combined_Plot <- Spring_NMDS_Main_Axis + Spring_NMDS_Secondary_Axis + 
  plot_layout(ncol = 2)
print(Spring_Combined_Plot)

##########################################################################################################################
## Objective 2a (Compare benthic macroinvertebrate abundance between Donor (Pre-transplant) and recipient sites (Post-transplant)
## Remove all taxa without established tolerance category 
Taxa_Data_Cleaned_II <- Taxa_Data %>%
  dplyr::select(-Chironomidae, -Entomobryidae, -Hydrachnidia, -Crambidae, -Neoleptophlebia, -Teloganopsis)

## Join taxa matrix back with Metadata to prepare for long pivot 
Taxa_Data_Cleaned_II_Full <- dplyr::bind_cols(Metadata, Taxa_Data_Cleaned_II)

## Pivot long the community matrix data so we can add tolerance information 
Taxa_Data_Cleaned_II_Long <- Taxa_Data_Cleaned_II_Full %>%
  pivot_longer(
    cols = -c(1:5),   # pivot all taxa columns
    names_to = "Lowest_Taxonomic_Unit",
    values_to = "Abundance"
  )

## Add Tolerance Data 
Taxa_Data_Cleaned_II_Long <- Taxa_Data_Cleaned_II_Long %>%
  left_join(Functional_Trait_Data %>% dplyr::select(Lowest_Taxonomic_Unit, Tolerance_Category), by = "Lowest_Taxonomic_Unit")

## Group by Cages, Site, Substrate, and Season 
Taxa_Data_Cleaned_II_Long <- Taxa_Data_Cleaned_II_Long %>%
  group_by(
    Cage_ID,
    Site_ID,
    Site_Type,
    Tolerance_Category, 
    Sampling_Time,
    Substrate_Type
  ) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

## Pivot The data wide so that we have the tolerance categories as columns 
Tolerance_Data <- Taxa_Data_Cleaned_II_Long %>%
  pivot_wider(
    names_from = Tolerance_Category,
    values_from = Abundance,
    values_fill = 0
  )

## Add watershed data 
Tolerance_Data <- Tolerance_Data %>%
  mutate(
    Watershed = case_when(
      Site_ID %in% c("Donor_02", "Recipient_02") ~ "Potomac_River",
      Site_ID %in% c("Donor_11", "Recipient_10") ~ "Little_Patuxent",
      Site_ID %in% c("Donor_13", "Recipient_11") ~ "Loch_Raven",
      TRUE ~ NA_character_
    )
  )

## Pivot the data long 
Tolerance_Data_Long <- Tolerance_Data %>%
  pivot_longer(
    cols = c(Sensitive, Moderately_Sensitive, Tolerant),
    names_to = "Sensitivity",
    values_to = "Abundance"
  )

## Fit a mixed GLM with site and sensitivity as the fixed factor
Full_Survival_Model <- glmmTMB(
  Abundance ~ Site_Type * Sensitivity,
  family = nbinom2(),
  data = Tolerance_Data_Long
)

summary(Full_Survival_Model)

## Run a pairwise comparison 
emmeans(Full_Survival_Model, pairwise ~ Site_Type | Sensitivity)

## Check the ANOVA Summary 
library(car)
Anova(Full_Survival_Model)

## Extract model coefficients and bind into a dataframe 
Full_Survival_Model_Coefficients <- emmeans(
  Full_Survival_Model,
  ~ Site_Type | Sensitivity,
  type = "response"   # back-transform from log scale
)
Full_Survival_Model_Coefficients_df <- as.data.frame(Full_Survival_Model_Coefficients)

## Order sensitivity categories 
Full_Survival_Model_Coefficients_df$Sensitivity <- factor(
  Full_Survival_Model_Coefficients_df$Sensitivity,
  levels = c("Sensitive", "Moderately_Sensitive", "Tolerant")
)

## Recode Donor and recipient as Pre-transplant and post-transplant before creating the plot 
Tolerance_Data_Long$Site_Type <- dplyr::recode(
  Tolerance_Data_Long$Site_Type,
  "Donor" = "Pre-transplant",
  "Recipient" = "Post-transplant"
)

Full_Survival_Model_Coefficients_df$Site_Type <- dplyr::recode(
  Full_Survival_Model_Coefficients_df$Site_Type,
  "Donor" = "Pre-transplant",
  "Recipient" = "Post-transplant"
)

## Build mean plot 
ggplot(Full_Survival_Model_Coefficients_df,
       aes(x = Sensitivity,
           y = response,
           color = Site_Type)) +
  
  ## Raw jittered points
  geom_jitter(
    data = Tolerance_Data_Long,
    aes(x = Sensitivity,
        y = Abundance,
        color = Site_Type),
    alpha = 0.35,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.01,
      dodge.width = 0.6
    )
  ) +
  
  ## Predicted means
  geom_point(
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  
  ## Error bars
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15,
    linewidth = 1,
    position = position_dodge(width = 0.6)
  ) +
  
  scale_color_manual(values = c(
    "Pre-transplant" = "darkred",
    "Post-transplant" = "black"
  )) +
  
  labs(
    x = "Sensitivity Category",
    y = "Abundance"
  ) +
  
  coord_cartesian(ylim = c(0, 110)) +
  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_blank()
  )

#######################################################################################################################################################
## Objective 2b (Model how macroinvertebrate abundance pre and post-transplant vary across selected watersheds)
## Subset the original data into sensitive, moderately sensitive, and tolerant. 
Sensitive_Data <- Tolerance_Data_Long %>%
  filter(Sensitivity == "Sensitive")

Moderately_Sensitive_Data <- Tolerance_Data_Long %>% 
  filter(Sensitivity == "Moderately_Sensitive")

Tolerant_Data <- Tolerance_Data_Long %>% 
  filter(Sensitivity == "Tolerant")

## Fit the models. 
Sensitive_Watershed_Model <- glmmTMB(
  Abundance ~ Site_Type * Watershed,
  family = nbinom2(),
  data = Sensitive_Data
)

Moderately_Sensitive_Watershed_Model <- glmmTMB(
  Abundance ~ Site_Type * Watershed,
  family = nbinom2(),
  data = Moderately_Sensitive_Data
)

Tolerant_Watershed_Model <- glmmTMB(
  Abundance ~ Site_Type * Watershed,
  family = nbinom2(),
  data = Tolerant_Data
)

## Get GLM Summaries 
summary(Sensitive_Watershed_Model)
summary(Moderately_Sensitive_Watershed_Model)
summary(Tolerant_Watershed_Model)

## Get the ANOVA summary for the models 
Anova(Sensitive_Watershed_Model)
Anova(Moderately_Sensitive_Watershed_Model)
Anova(Tolerant_Watershed_Model)


## Get Pairwise Contrasts 
emmeans(Sensitive_Watershed_Model, pairwise ~ Site_Type | Watershed)
emmeans(Moderately_Sensitive_Watershed_Model, pairwise ~ Site_Type | Watershed)
emmeans(Tolerant_Watershed_Model, pairwise ~ Site_Type | Watershed)

## Extract model coefficients and bind into a dataframe 
Sensitive_Watershed_Model_Coefficients <- emmeans(
  Sensitive_Watershed_Model,
  ~ Site_Type | Watershed,
  type = "response"   # back-transform from log scale
)
Sensitive_Watershed_Model_Coefficients_df <- as.data.frame(Sensitive_Watershed_Model_Coefficients)
pairs(Sensitive_Watershed_Model_Coefficients)


####
Moderately_Sensitive_Watershed_Model_Coefficients <- emmeans(
  Moderately_Sensitive_Watershed_Model,
  ~ Site_Type | Watershed,
  type = "response"   # back-transform from log scale
)
Moderately_Sensitive_Watershed_Model_Coefficients_df <- as.data.frame(Moderately_Sensitive_Watershed_Model_Coefficients)
pairs(Moderately_Sensitive_Watershed_Model_Coefficients)

###
Tolerant_Watershed_Model_Coefficients <- emmeans(
  Tolerant_Watershed_Model,
  ~ Site_Type | Watershed,
  type = "response"   # back-transform from log scale
)
Tolerant_Watershed_Model_Coefficients_df <- as.data.frame(Tolerant_Watershed_Model_Coefficients)
pairs(Tolerant_Watershed_Model_Coefficients)

## Build mean plots 
pd <- position_dodge(width = 0.4)

## Reorder factor levels before plotting 
Moderately_Sensitive_Data$Site_Type <- factor(
  Moderately_Sensitive_Data$Site_Type,
  levels = c("Pre-transplant", "Post-transplant")
)

Moderately_Sensitive_Watershed_Model_Coefficients_df$Site_Type <- factor(
  Moderately_Sensitive_Watershed_Model_Coefficients_df$Site_Type,
  levels = c("Pre-transplant", "Post-transplant")
)


Sensitive_Data$Site_Type <- factor(
  Sensitive_Data$Site_Type,
  levels = c("Pre-transplant", "Post-transplant")
)

Sensitive_Watershed_Model_Coefficients_df$Site_Type <- factor(
  Sensitive_Watershed_Model_Coefficients_df$Site_Type,
  levels = c("Pre-transplant", "Post-transplant")
)

Tolerant_Data$Site_Type <- factor(
  Tolerant_Data$Site_Type,
  levels = c("Pre-transplant", "Post-transplant")
)

Tolerant_Watershed_Model_Coefficients_df$Site_Type <- factor(
  Tolerant_Watershed_Model_Coefficients_df$Site_Type,
  levels = c("Pre-transplant", "Post-transplant")
)

## Moderately Sensitive 
Moderately_Sensitive_Plot <- ggplot(
  Moderately_Sensitive_Watershed_Model_Coefficients_df,
  aes(x = Watershed, y = response, color = Site_Type)
) +
  
  geom_jitter(
    data = Moderately_Sensitive_Data,
    aes(x = Watershed, y = Abundance, color = Site_Type),
    alpha = 0.35,
    size = 2,
    position = pd
  ) +
  
  geom_point(size = 4, position = pd) +
  
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15,
    linewidth = 1,
    position = pd
  ) +
  
  scale_color_manual(values = c(
    "Pre-transplant" = "darkred",
    "Post-transplant" = "black"
  )) +
  
  labs(
    x = "Watersheds",
    y = "Mod_Sensitive_Taxa_Abundance"
  ) +
  
  coord_cartesian(ylim = c(0, 170)) +
  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_blank()
  )
print(Moderately_Sensitive_Plot)

## Sensitive 
Sensitive_Plot <- ggplot(
  Sensitive_Watershed_Model_Coefficients_df,
  aes(x = Watershed, y = response, color = Site_Type)
) +
  
  geom_jitter(
    data = Sensitive_Data,
    aes(x = Watershed, y = Abundance, color = Site_Type),
    alpha = 0.35,
    size = 2,
    position = pd
  ) +
  
  geom_point(size = 4, position = pd) +
  
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15,
    linewidth = 1,
    position = pd
  ) +
  
  scale_color_manual(values = c(
    "Pre-transplant" = "darkred",
    "Post-transplant" = "black"
  )) +
  
  labs(
    x = "Watersheds",
    y = "Sensitive_Taxa_Abundance"
  ) +
  
  coord_cartesian(ylim = c(0, 75)) +
  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_blank()
  )
print(Sensitive_Plot)

## Tolerant 
Tolerant_Plot <- ggplot(
  Tolerant_Watershed_Model_Coefficients_df,
  aes(x = Watershed, y = response, color = Site_Type)
) +
  
  geom_jitter(
    data = Tolerant_Data,
    aes(x = Watershed, y = Abundance, color = Site_Type),
    alpha = 0.35,
    size = 2,
    position = pd
  ) +
  
  geom_point(size = 4, position = pd) +
  
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15,
    linewidth = 1,
    position = pd
  ) +
  
  scale_color_manual(values = c(
    "Pre-transplant" = "darkred",
    "Post-transplant" = "black"
  )) +
  
  labs(
    x = "Watersheds",
    y = "Tolerant_Taxa_Abundance"
  ) +
  
  coord_cartesian(ylim = c(0, 80)) +
  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_blank()
  )
print(Tolerant_Plot)

## Use Cowplot to combine the plots 
library(cowplot)

legend <- get_legend(
  Sensitive_Plot +
    theme(legend.position = "right")
)

Sensitive_no_legend <- Sensitive_Plot + theme(legend.position = "none")
Moderately_no_legend <- Moderately_Sensitive_Plot + theme(legend.position = "none")
Tolerant_no_legend <- Tolerant_Plot + theme(legend.position = "none")

top_row <- plot_grid(
  Sensitive_no_legend,
  Tolerant_no_legend,
  ncol = 2
)

bottom_row <- plot_grid(
  NULL,
  Moderately_no_legend,
  NULL,
  ncol = 3,
  rel_widths = c(1, 2, 1)
)

plots_only <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 0.9)
)

Full_Watershed_Plot <- plot_grid(
  plots_only,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.15)
)

print(Full_Watershed_Plot)

##################################################################################################################################
## Objective 2c: Model individual taxa performance 
## Run a GLM for individual macroinvertebrate taxa using site_type as a fixed predictor. 
## Split data into spring and fall because most taxa were only present in one season.  
Fall_Taxon_Data   <- Taxa_Data_Cleaned_II_Full %>% filter(Sampling_Time == "Fall")
Spring_Taxon_Data <- Taxa_Data_Cleaned_II_Full %>% filter(Sampling_Time == "Spring")

## Remove all taxa from the fall and spring dataset that did not appear in at least five tubes. 
Fall_Taxon_Data <- Fall_Taxon_Data %>%
  dplyr::select(1:5, where(~ sum(.x > 0, na.rm = TRUE) >= 5))

Spring_Taxon_Data <- Spring_Taxon_Data %>% 
  dplyr::select(1:5, where(~sum(.x > 0, na.rm = TRUE) >= 5))

## Pivot your community matrix data long 
Fall_Taxon_Long <- Fall_Taxon_Data %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = "Lowest_Taxonomic_Unit",
    values_to = "Abundance"
  )

Spring_Taxon_Long <- Spring_Taxon_Data %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = "Lowest_Taxonomic_Unit",
    values_to = "Abundance"
  )

## Fit fall models 
Fall_GLMs <- Fall_Taxon_Long %>%
  group_by(Lowest_Taxonomic_Unit) %>%
  group_map(~ {
    m <- glm.nb(Abundance ~ Site_Type, data = .x)
    broom::tidy(m) %>%
      mutate(Lowest_Taxonomic_Unit = .y$Lowest_Taxonomic_Unit)
  }) %>%
  bind_rows() %>%
  filter(term == "Site_TypeRecipient")

## Fit spring models 
Spring_GLMs <- Spring_Taxon_Long %>%
  group_by(Lowest_Taxonomic_Unit) %>%
  group_map(~ {
    m <- glm.nb(Abundance ~ Site_Type, data = .x)
    broom::tidy(m) %>%
      mutate(Lowest_Taxonomic_Unit = .y$Lowest_Taxonomic_Unit)
  }) %>%
  bind_rows() %>%
  filter(term == "Site_TypeRecipient")

## Extract fall coefficients 
Fall_Coefficients <- Fall_GLMs %>%
  filter(term == "Site_TypeRecipient") %>%
  mutate(
    Ratio = exp(estimate),
    
    ## Lower bounds for ratio 
    Ratio_Lower = exp(estimate - 1.96 * std.error),
    Ratio_Upper = exp(estimate + 1.96 * std.error),
    
    ## Bounds for estimate (real CI)
    Lower = estimate - 1.96 * std.error,
    Upper = estimate + 1.96 * std.error,
    
    ## Percentage change and bounds 
    Percent_Change = (Ratio - 1) * 100,
    Percent_Change_Lower = (Ratio_Lower - 1) * 100,
    Percent_Change_Upper = (Ratio_Upper - 1) * 100
  )

## Remove Plectrocnemia (Has a ridiculously high CI)
Fall_Coefficients <- Fall_Coefficients[-15, ]

Fall_Coefficients <- Fall_Coefficients %>%
  arrange(desc(estimate)) %>% 
  mutate(Lowest_Taxonomic_Unit = factor(Lowest_Taxonomic_Unit,
                                        levels = rev(Lowest_Taxonomic_Unit))) ## Rearrange taxa order 

## Add Tolerance_category data to the dataframe 
Fall_Coefficients <- Fall_Coefficients %>%
  left_join(
    Functional_Trait_Data %>% 
      dplyr::select(Lowest_Taxonomic_Unit, Tolerance_Category),
    by = "Lowest_Taxonomic_Unit"
  )

## Reorder Lowest taxonomic unit 
Fall_Coefficients <- Fall_Coefficients %>%
  mutate(
    Lowest_Taxonomic_Unit = factor(
      Lowest_Taxonomic_Unit,
      levels = Fall_Coefficients$Lowest_Taxonomic_Unit[order(Fall_Coefficients$estimate, decreasing = TRUE)]
    )
  )

Fall_Coefficients <- Fall_Coefficients %>%
  dplyr::mutate(
    Lowest_Taxonomic_Unit = forcats::fct_reorder(
      Lowest_Taxonomic_Unit,
      estimate,
      .desc = FALSE   # increasing order
    )
  )

## Build Fall Plot (Log_ Scale plot)
ggplot(Fall_Coefficients,
       aes(x = estimate, y = Lowest_Taxonomic_Unit)) +
  
  # Points stay fixed color
  geom_point(size = 4, color = "black") +
  
  # CI lines get mapped color
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper, color = Tolerance_Category),
    height = 0.2,
    linewidth = 2
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Sensitive" = "#1F968BFF",
    "Moderately_Sensitive" = "#F98E3B",
    "Tolerant" = "#453781FF"
  )) +
  
  labs(
    x = "Effect of Transplant (log scale)",
    y = "Taxon"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 19)
  )

## Extract Spring Coefficients 
Spring_Coefficients <- Spring_GLMs %>%
  filter(term == "Site_TypeRecipient") %>%
  mutate(
    Ratio = exp(estimate),
    
    ## Lower bounds for ratio 
    Ratio_Lower = exp(estimate - 1.96 * std.error),
    Ratio_Upper = exp(estimate + 1.96 * std.error),
    
    ## Bounds for estimate (real CI)
    Lower = estimate - 1.96 * std.error,
    Upper = estimate + 1.96 * std.error,
    
    ## Percentage change and bounds 
    Percent_Change = (Ratio - 1) * 100,
    Percent_Change_Lower = (Ratio_Lower - 1) * 100,
    Percent_Change_Upper = (Ratio_Upper - 1) * 100
  )

## Add Tolerance_category data to the dataframe 
Spring_Coefficients <- Spring_Coefficients %>%
  left_join(
    Functional_Trait_Data %>% 
      dplyr::select(Lowest_Taxonomic_Unit, Tolerance_Category),
    by = "Lowest_Taxonomic_Unit"
  )

## Remove Hemerodromia and Probezzia because they have ridiculously high standard errors and CI. 
Spring_Coefficients <- Spring_Coefficients[-9, ] ## Hemerodromia
Spring_Coefficients <- Spring_Coefficients[-18, ] ## Probezzia

## Reorder Lowest taxonomic unit 
Spring_Coefficients <- Spring_Coefficients %>%
  mutate(
    Lowest_Taxonomic_Unit = factor(
      Lowest_Taxonomic_Unit,
      levels = Spring_Coefficients$Lowest_Taxonomic_Unit[order(Spring_Coefficients$estimate, decreasing = TRUE)]
    )
  )

Spring_Coefficients <- Spring_Coefficients %>%
  dplyr::mutate(
    Lowest_Taxonomic_Unit = forcats::fct_reorder(
      Lowest_Taxonomic_Unit,
      estimate,
      .desc = FALSE   # increasing order
    )
  )

## Build Spring Coefficient plot 
ggplot(Spring_Coefficients,
       aes(x = estimate, y = Lowest_Taxonomic_Unit)) +
  
  # Points stay fixed color
  geom_point(size = 4, color = "black") +
  
  # CI lines get mapped color
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper, color = Tolerance_Category),
    height = 0.2,
    linewidth = 2
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Sensitive" = "#1F968BFF",
    "Moderately_Sensitive" = "#F98E3B",
    "Tolerant" = "#453781FF"
  )) +
  
  labs(
    x = "Effect of Transplant (log scale)",
    y = "Taxon"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 19)
  )

##################################################################################################################################
## Objective 3 : Compare community composition, diversity and abundance between leaf and rock substrates. 
## Clean up data for this analysis 
## We only need donor data for this, so lets extract donor site data from the full dataset 

Substrate_Analysis_Data <- Full_Survival_Experiment_Data %>% 
  filter(Site_Type == "Donor")

### Separate Donor site data into fall and spring data 
Fall_Substrate_Data <- Substrate_Analysis_Data %>%
  filter(Sampling_Time == "Fall")

Spring_Substrate_Data <- Substrate_Analysis_Data %>% 
  filter(Sampling_Time == "Spring")

## Separate into matrix and metadata 
Fall_Substrate_Matrix <- Fall_Substrate_Data[, 6:74]
Fall_Substrate_Metadata <- Fall_Substrate_Data[, 1:5]
Spring_Substrate_Matrix <- Spring_Substrate_Data[, 6:74]
Spring_Substrate_Metadata <- Spring_Substrate_Data[, 1:5]

## Clean Fall and Spring data with the following steps 
## Remove columns that sums up to zero (They represent taxa not present for that season)
Fall_Substrate_Matrix_Cleaned <- Fall_Substrate_Matrix[, colSums(Fall_Substrate_Matrix) != 0] ## 28 taxa removed 

## Remove taxa that are not present in at least 3 samples 
Fall_Substrate_Matrix_Cleaned <- Fall_Substrate_Matrix_Cleaned[, colSums(Fall_Substrate_Matrix_Cleaned > 0) >= 3] ## 24 taxa removed

## Transform data to reduce the effect of highly dominant taxa 
Fall_Substrate_Matrix_Transformed <- sqrt(Fall_Substrate_Matrix_Cleaned)

## Compute bray-curtis distance and run permdsip
Fall_Substrate_Bray_Curtis <- vegdist(Fall_Substrate_Matrix_Transformed, method = "bray") ## Calculates bray curtis distance 
Fall_Substrate_Matrix_Dispersion <- betadisper(Fall_Substrate_Bray_Curtis, Fall_Substrate_Metadata$Substrate_Type)
anova(Fall_Substrate_Matrix_Dispersion) ## F = 1.2052, p= 0.2816

## Run PERMANOVA
Fall_Substrate_PERMANOVA <- adonis2(Fall_Substrate_Bray_Curtis ~ Substrate_Type, data = Fall_Substrate_Metadata)
Fall_Substrate_PERMANOVA ## F = 1.9371, p = 0.107

## Run and plot NMDS 
Fall_Substrate_NMDS <- metaMDS(
  Fall_Substrate_Matrix_Transformed,
  distance = "bray",
  k = 3,
  trymax = 1000
) ## Converged after 35 iterations 

Fall_Substrate_NMDS$stress ## 0.11

## Extract NMDS Scores and calculate centroids 
Fall_Substrate_NMDS_Scores <- as.data.frame(scores(Fall_Substrate_NMDS, display = "sites"))
Fall_Substrate_NMDS_Scores$Substrate_Type <- Fall_Substrate_Metadata$Substrate_Type

## Compute Centroids 
Fall_Substrate_Centroids <- Fall_Substrate_NMDS_Scores %>%
  group_by(Substrate_Type) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    NMDS3 = mean(NMDS3)
  )

## Run substrate analysis NMDS for fall 
## Primary axis
Fall_Substrate_NMDS_Main_Axis <- ggplot(
  Fall_Substrate_NMDS_Scores,
  aes(x = NMDS1, y = NMDS2, color = Substrate_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Substrate_Type
  stat_ellipse(aes(group = Substrate_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Substrate_Type
  geom_point(
    data = Fall_Substrate_Centroids,
    aes(x = NMDS1, y = NMDS2, color = Substrate_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Leaf vs Rock
  scale_color_manual(values = c("Leaf" = "#556B2F",
                                "Rock" = "goldenrod")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS2") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) + 
  theme(legend.position = "none")

## Secondary axis 
Fall_Substrate_NMDS_Secondary_Axis <- ggplot(
  Fall_Substrate_NMDS_Scores,
  aes(x = NMDS1, y = NMDS3, color = Substrate_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Site_Type (your main inference groups)
  stat_ellipse(aes(group = Substrate_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Site_Type
  geom_point(
    data = Fall_Substrate_Centroids,
    aes(x = NMDS1, y = NMDS3, color = Substrate_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Donor vs Recipient
  scale_color_manual(values = c("Leaf" = "#556B2F",
                                "Rock" = "goldenrod")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS3") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  )

print(Fall_Substrate_NMDS_Main_Axis)
print(Fall_Substrate_NMDS_Secondary_Axis)

## Combine figures with patchwork 
library(patchwork)

Fall_Substrate_Combined_Plot <- Fall_Substrate_NMDS_Main_Axis + Fall_Substrate_NMDS_Secondary_Axis + 
  plot_layout(ncol = 2)
print(Fall_Substrate_Combined_Plot)

################################################# Spring 
## Clean Spring data with the following steps 
## Remove columns that sums up to zero (They represent taxa not present for that season)
Spring_Substrate_Matrix_Cleaned <- Spring_Substrate_Matrix[, colSums(Spring_Substrate_Matrix) != 0] ## 28 taxa removed 

## Remove taxa that are not present in at least 3 samples 
Spring_Substrate_Matrix_Cleaned <- Spring_Substrate_Matrix_Cleaned[, colSums(Spring_Substrate_Matrix_Cleaned > 0) >= 3] ## 22 taxa removed

## Transform data to reduce the effect of highly dominaant taxa 
Spring_Substrate_Matrix_Transformed <- sqrt(Spring_Substrate_Matrix_Cleaned)

## Compute bray-curtis distance and run permdsip
Spring_Substrate_Bray_Curtis <- vegdist(Spring_Substrate_Matrix_Transformed, method = "bray") ## Calculates bray curtis distance 
Spring_Substrate_Matrix_Dispersion <- betadisper(Spring_Substrate_Bray_Curtis, Spring_Substrate_Metadata$Substrate_Type)
anova(Spring_Substrate_Matrix_Dispersion) ## F = 1.9242, p= 0.1763

## Run PERMANOVA
Spring_Substrate_PERMANOVA <- adonis2(Spring_Substrate_Bray_Curtis ~ Substrate_Type, data = Spring_Substrate_Metadata)
Spring_Substrate_PERMANOVA ## F = 1.5663, p = 0.144

## Run and plot NMDS 
Spring_Substrate_NMDS <- metaMDS(
  Spring_Substrate_Matrix_Transformed,
  distance = "bray",
  k = 3,
  trymax = 1000
)

Spring_Substrate_NMDS$stress ## 0.12

## Extract NMDS Scores and calculate centroids 
Spring_Substrate_NMDS_Scores <- as.data.frame(scores(Spring_Substrate_NMDS, display = "sites"))
Spring_Substrate_NMDS_Scores$Substrate_Type <- Spring_Substrate_Metadata$Substrate_Type

## Compute Centroids 
Spring_Substrate_Centroids <- Spring_Substrate_NMDS_Scores %>%
  group_by(Substrate_Type) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    NMDS3 = mean(NMDS3)
  )

## Run substrate analysis NMDS for fall 
## Primary Axis
Spring_Substrate_NMDS_Main_Axis <- ggplot(
  Spring_Substrate_NMDS_Scores,
  aes(x = NMDS1, y = NMDS2, color = Substrate_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Substrate_Type
  stat_ellipse(aes(group = Substrate_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Substrate_Type
  geom_point(
    data = Spring_Substrate_Centroids,
    aes(x = NMDS1, y = NMDS2, color = Substrate_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Leaf vs Rock
  scale_color_manual(values = c("Leaf" = "#556B2F",
                                "Rock" = "goldenrod")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS2") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) + 
  theme(legend.position = "none")

## Secondary axis
Spring_Substrate_NMDS_Secondary_Axis <- ggplot(
  Spring_Substrate_NMDS_Scores,
  aes(x = NMDS1, y = NMDS3, color = Substrate_Type)
) +
  
  # Points
  geom_point(size = 4, stroke = 1, alpha = 0.9) +
  
  # 95% ellipses by Site_Type (your main inference groups)
  stat_ellipse(aes(group = Substrate_Type),
               level = 0.95,
               linewidth = 1) +
  
  # Centroids for Site_Type
  geom_point(
    data = Spring_Substrate_Centroids,
    aes(x = NMDS1, y = NMDS3, color = Substrate_Type),
    size = 7,
    shape = 23,
    stroke = 1.5
  ) +
  
  # Color palette for Donor vs Recipient
  scale_color_manual(values = c("Leaf" = "#556B2F",
                                "Rock" = "goldenrod")) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  labs(x = "NMDS1",
       y = "NMDS3") +
  
  theme(
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  ) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  )

print(Spring_Substrate_NMDS_Main_Axis)
print(Spring_Substrate_NMDS_Secondary_Axis)

## Combine figures with patchwork 
Spring_Substrate_Combined_Plot <- Spring_Substrate_NMDS_Main_Axis + Spring_Substrate_NMDS_Secondary_Axis + 
  plot_layout(ncol = 2)
print(Spring_Substrate_Combined_Plot)

###########################################################################################################################
## Objective 3b: Compare taxa richness, %EPT, Shannon diversity, and abundance 
Substrate_Analysis_Matrix <- Substrate_Analysis_Data [, 6:74]
Substrate_Analysis_Metadata <- Substrate_Analysis_Data[, 1:5]

## Compute taxa richness per cage 
Taxa_Richness <- specnumber(Substrate_Analysis_Matrix)

## Combine with metadata 
Taxa_Richness_Data <- data.frame(
  Richness = Taxa_Richness,
  Substrate_Type = Substrate_Analysis_Metadata$Substrate_Type,
  Sampling_Time = Substrate_Analysis_Metadata$Sampling_Time,
  Site_Type = Substrate_Analysis_Metadata$Site_Type
)

## Fit a taxa richness model 
Taxa_Richness_Model <- lm(Richness ~ Substrate_Type * Sampling_Time, data = Taxa_Richness_Data)
summary(Taxa_Richness_Model) ## Site - season interaction not significant p = 0.7204, refit a simpler model without the interaction and season 

Taxa_Richness_Model_Simpler <- lm(Richness ~ Substrate_Type, data = Taxa_Richness_Data)
summary(Taxa_Richness_Model_Simpler)
Anova(Taxa_Richness_Model_Simpler, type = "II")

##Check normality of residuals 
res <- residuals(Taxa_Richness_Model_Simpler)
fit <- fitted(EPT_Model)
qqnorm(res)
qqline(res, col = "red", lwd = 2)
hist(res, breaks = 10, col = "grey", main = "Residuals Histogram")
shapiro.test(res) ## Residuals are normally distributed 

## Build box_plot_stats 
Taxa_Richness_box_stats <- Taxa_Richness_Data %>%
  group_by(Substrate_Type) %>%
  summarise(
    ymin  = min(Richness),
    lower = quantile(Richness, 0.25),
    middle = mean(Richness),      ## mean instead of median
    upper = quantile(Richness, 0.75),
    ymax  = max(Richness)
  )
## Build box plot
Taxa_Richness_Plot <- ggplot() +
  geom_boxplot(
    data = Taxa_Richness_box_stats,
    aes(x = Substrate_Type,
        ymin = ymin,
        lower = lower,
        middle = middle,   ## mean
        upper = upper,
        ymax = ymax,
        fill = Substrate_Type),
    stat = "identity",
    alpha = 0.7
  ) +
  
  geom_jitter(
    data = Taxa_Richness_Data,
    aes(x = Substrate_Type, y = Richness, color = Substrate_Type),
    width = 0.15,
    alpha = 0.5
  ) +
  
  scale_fill_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  scale_color_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  
  theme_bw(base_size = 16)

## Compute Shannon Diversity per cage 
Shannon_Diversity <- diversity(Substrate_Analysis_Matrix, index = "shannon") ## Shannon diversity computed for all 60 cages 
## Combine with Metadata 
Shannon_Data <- data.frame(
  Shannon = Shannon_Diversity,
  Substrate_Type = Substrate_Analysis_Metadata$Substrate_Type,
  Sampling_Time = Substrate_Analysis_Metadata$Sampling_Time,
  Site_Type = Substrate_Analysis_Metadata$Site_Type
)

## Fit a linear model to compare Shannon diversity between leaf and rock (Add site-type as a fixed factor to the model to see if the effect is consistent before and after transplant)
Shannon_Model <- lm(Shannon ~ Substrate_Type * Sampling_Time, data = Shannon_Data)
summary(Shannon_Model)
anova(Shannon_Model)
Shannon_Model_Simpler <- lm(Shannon ~ Substrate_Type, data = Shannon_Data)
summary(Shannon_Model_Simpler)

## Check the normality of the residuals 
par(mfrow = c(2,2))
plot(Shannon_Model) ## Residuals are normally distributed 
par(mfrow = c(1,1)) ## reset the canvas 

## Build boxplot stats to override median with mean. 
Shannon_box_stats <- Shannon_Data %>%
  group_by(Substrate_Type) %>%
  summarise(
    ymin  = min(Shannon),
    lower = quantile(Shannon, 0.25),
    middle = mean(Shannon),      # mean instead of median
    upper = quantile(Shannon, 0.75),
    ymax  = max(Shannon)
  )

## Build box plot
Shannon_Plot <- ggplot() +
  geom_boxplot(
    data = Shannon_box_stats,
    aes(x = Substrate_Type,
        ymin = ymin,
        lower = lower,
        middle = middle,   # <-- mean
        upper = upper,
        ymax = ymax,
        fill = Substrate_Type),
    stat = "identity",
    alpha = 0.7
  ) +
  
  geom_jitter(
    data = Shannon_Data,
    aes(x = Substrate_Type, y = Shannon_Diversity, color = Substrate_Type),
    width = 0.15,
    alpha = 0.5
  ) +
  
  scale_fill_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  scale_color_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  
  theme_bw(base_size = 16)
    
## Fit a linear model to test substrate effect on macroinvertebrate abundance  
Abundance_Data <- rowSums(Substrate_Analysis_Matrix)

Abundance_Data <- data.frame(
  Substrate_Type   = Substrate_Analysis_Metadata$Substrate_Type,
  Sampling_Time = Substrate_Analysis_Metadata$Sampling_Time,
  Abundance  = Abundance_Data
)

## Do a quick visual normality check 
hist(Abundance_Data$Abundance) ## Data appears normally distributed, residuals should be too 

## Fit a linear model 
Abundance_Model <- lm(Abundance ~ Substrate_Type * Sampling_Time, data = Abundance_Data)
summary(Abundance_Model)
Abundance_Model_Simpler <- lm(Abundance ~ Substrate_Type, data = Abundance_Data)
summary(Abundance_Model_Simpler)
Anova (Abundance_Model_Simpler, type = "II")

## Check residuals 
par(mfrow = c(1,2))
hist(residuals(Abundance_Model_Simpler))
qqnorm(residuals(Abundance_Model_Simpler))
plot(Abundance_Model_Simpler, which = 1) ## Homogeneity looks good 

## Build boxplot stats to override median with mean. 
box_stats <- Abundance_Data %>%
  group_by(Substrate_Type) %>%
  summarise(
    ymin  = min(Abundance),
    lower = quantile(Abundance, 0.25),
    middle = mean(Abundance),      ## mean instead of median
    upper = quantile(Abundance, 0.75),
    ymax  = max(Abundance)
  )

## Build box plot
Abundance_Plot <- ggplot() +
  geom_boxplot(
    data = box_stats,
    aes(x = Substrate_Type,
        ymin = ymin,
        lower = lower,
        middle = middle,   ## mean
        upper = upper,
        ymax = ymax,
        fill = Substrate_Type),
    stat = "identity",
    alpha = 0.7
  ) +
  
  geom_jitter(
    data = Abundance_Data,
    aes(x = Substrate_Type, y = Abundance, color = Substrate_Type),
    width = 0.15,
    alpha = 0.5
  ) +
  
  scale_fill_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  scale_color_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  
  theme_bw(base_size = 16)

## Fit a linear model to test substrate effect on EPT_Richness
## Create EPT taxa list 
EPT_Taxa <- c("Allocapnia", "Maccaffertium", "Dolophilodes", "Cheumatopsyche", "Isonychia", "Plectrocnemia", "Chimarra", "Ameletus", "Diplectrona", "Eccoptura", "Hydropsyche",
              "Baetisca", "Isoperla", "Psilotreta", "Ptilostomis", "Molanna", "Leptophlebia", "Nyctiophylax", "Rhyacophila", "Neureclipsis", "Eurylophella", "Glossosoma", "Neophylax",
              "Hydatophylax", "Lepidostoma", "Baetis", "Ironoquia", "Timpanoga", "Polycentropus", "Amphinemura", "Penelomax", "Ephemerella")

EPT_Richness <- specnumber(Substrate_Analysis_Matrix[, EPT_Taxa])

EPT_Taxa_Data <- data.frame(
  Substrate_Type   = Substrate_Analysis_Metadata$Substrate_Type,
  Sampling_Time = Substrate_Analysis_Metadata$Sampling_Time,
  EPT_Richness  = EPT_Richness
)

## Fit Substrate - Percent_EPT Model 
EPT_Model <- glm(EPT_Richness ~ Substrate_Type,
                          data = EPT_Taxa_Data,
                          family = poisson)

summary(EPT_Model)
anova(EPT_Model)
library(AER)
dispersiontest(EPT_Model) ## Test to see if poisson is the right fit 

## Check normality of residuals 
res <- residuals(EPT_Model)
fit <- fitted(EPT_Model)
qqnorm(res)
qqline(res, col = "red", lwd = 2)
hist(res, breaks = 10, col = "grey", main = "Residuals Histogram")
shapiro.test(res) ## Residuals are not normally distributed. Use Poisson GLM.
library(AER)
dispersiontest(EPT_Model) ## Test to see if poisson is the right fit

## Build boxplot stats to override median with mean. 
EPT_box_stats <- EPT_Taxa_Data %>%
  group_by(Substrate_Type) %>%
  summarise(
    ymin  = min(EPT_Richness),
    lower = quantile(EPT_Richness, 0.25),
    middle = mean(EPT_Richness),      ## mean instead of median
    upper = quantile(EPT_Richness, 0.75),
    ymax  = max(EPT_Richness)
  )

## Build box plot
EPT_Plot <- ggplot() +
  geom_boxplot(
    data = EPT_box_stats,
    aes(x = Substrate_Type,
        ymin = ymin,
        lower = lower,
        middle = middle,   ## mean
        upper = upper,
        ymax = ymax,
        fill = Substrate_Type),
    stat = "identity",
    alpha = 0.7
  ) +
  
  geom_jitter(
    data = EPT_Taxa_Data,
    aes(x = Substrate_Type, y = EPT_Richness, color = Substrate_Type),
    width = 0.15,
    alpha = 0.5
  ) +
  
  scale_fill_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  scale_color_manual(values = c("Leaf" = "#556B2F", "Rock" = "goldenrod")) +
  
  theme_bw(base_size = 16)

## Merge figures with patchwork 
Substrate_Models <- (Taxa_Richness_Plot | Shannon_Plot) /
  (EPT_Plot | Abundance_Plot) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

Substrate_Models
