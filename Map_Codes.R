## Author: Ibrahim Fagbohun 
## Date: 04/04/26
# Aim: Provide a map showing study sites (Donor and recipients and the land use around the areas)

## Load relevant libraries 
library(sf)
library(ggplot2)
library(dplyr)
library(leaflet)
library(readxl)
library(ggthemes)
library(usmap)
library(viridisLite)
library(openxlsx)
library(ggnewscale)
library(lwgeom)
library(ggspatial)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)

## Load relevant data 
Survival_Experiment_Sites <- read_excel("~/OneDrive - The Pennsylvania State University/My Research/Data/CBT_Survival_Data/Cleaned_CBT_Survival_Data/Final_Survival_Experiment_Sites.xlsx")
Maryland_Watersheds <- st_read("Maryland_Watersheds/Maryland_Watersheds.shp") ## Maryland Watershed boundaries downloaded from Maryland DNR Website. 
Maryland_Flowlines <- st_read("Maryland_Flowlines/Rivers_and_Streams_-_Detailed.shp") ## Maryland flowlines downloaded from Maryland DNR Website. 
Maryland_Land_Use <- st_read("Maryland_Land_Use/Maryland_Land_Use_Land_Cover_-_Land_Use_Land_Cover_2010.shp") ## Land use data downloaded from Maryland DNR Website. 
Maryland_Counties <- st_read("Maryland_Counties/Maryland_Physical_Boundaries_-_County_Boundaries_(Detailed).shp") ## County boundaries downloaded from Maryland DNR website. 

## Convert needed dataset to sf simple feature (sf) object 
Survival_Experiment_Sites_sf <- st_as_sf(Survival_Experiment_Sites, coords = c("Longitude", "Latitude"), crs = 4326)

## Confirm the format that all dataset are in 
st_crs(Maryland_Watersheds) ##EPSG:3857
st_crs(Maryland_Flowlines) ##EPSG:3857
st_crs(Maryland_Land_Use) ##EPSG:3857
st_crs(Survival_Experiment_Sites_sf) ##EPSG:4326 

## Transform all layer into UTM 18N 
Maryland_Watersheds_utm <- st_transform(Maryland_Watersheds, 26918)
Maryland_Flowlines_utm  <- st_transform(Maryland_Flowlines, 26918)
Maryland_Land_Use_utm   <- st_transform(Maryland_Land_Use, 26918)
Survival_Experiment_Sites_utm <- st_transform(Survival_Experiment_Sites_sf, 26918)
Maryland_Counties_utm <- st_transform(Maryland_Counties, 26918)

## Regroup Land use Categories 
Maryland_Land_Use_utm <- Maryland_Land_Use_utm %>%
  mutate(
    Relevant_Categories = case_when(
      Descriptio %in% c(
        "Low Density Residential", "Very Low Density Residential",
        "Medium Density Residential", "High Density Residential",
        "Commercial", "Industrial", "Institutional",
        "Transportation", "Other Developed Lands"
      ) ~ "Urban",
      
      Descriptio == "Agriculture" ~ "Agriculture",
      Descriptio == "Forest" ~ "Forest",
      
      TRUE ~ NA_character_   # drop Water, Wetlands, Barren, etc.
    )
  )

## Filter Maryland counties to Howard, Baltimore, and Montgomery 
Study_Counties <- Maryland_Counties %>%
  filter(COUNTY %in% c("Baltimore", "Howard", "Montgomery"))

# Combine counties 
Study_Region <- Study_Counties %>%
  st_union() %>%
  st_as_sf()

## Extract Target watersheds from the entire watershed data 
Target_Watersheds <- Maryland_Watersheds_utm %>%
  filter(mde8name %in% c(
    "Seneca Creek",
    "Loch Raven Reservoir",
    "Jones Falls",
    "Anacostia River",
    "Potomac River MO Cnty",
    "Rock Creek",
    "Little Patuxent River"
  ))

## Extract the Chesapeake Bay Polygon 
Chesapeake_Bay <- Maryland_Watersheds_utm %>%
  filter(grepl("Chesapeake", mde8name, ignore.case = TRUE))

Chesapeake_Bay_utm <- st_transform(Chesapeake_Bay, 26918) ## Transform to the version as others 

## Combine the watersheds 
Target_Watershed_Region <- Target_Watersheds %>%
  st_union() %>%
  st_as_sf()

## Clip flow lines to the target watershed 
Flowlines_Target <- st_intersection(Maryland_Flowlines_utm, Target_Watershed_Region)

## Clip sites to watershed region 
Sites_Target <- st_intersection(Survival_Experiment_Sites_utm, Target_Watershed_Region)

## Clip land use to target watershed
LandUse_Target <- Maryland_Land_Use_utm |>
  dplyr::filter(COUNTY %in% c("BACO", "HOWA", "MONT"),
                !is.na(Relevant_Categories))

## Assign colors for land use categories 
Land_Use_Colors <- c(
  "Forest"      = "#A8D5A2",
  "Agriculture" = "#F9E79F",
  "Urban"       = "#F6B8B0"
)

## Plot everything in the county map 
Main_Map <- ggplot() +
  geom_sf(data = Chesapeake_Bay_utm, fill = "#5DADE2", color = NA) +
  geom_sf(data = LandUse_Target, aes(fill = Relevant_Categories), color = NA) +
  scale_fill_manual(values = Land_Use_Colors, name = "Land Use") +
  geom_sf(data = Study_Counties, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = Flowlines_Target, color = "#2C73D2", linewidth = 0.4) +
  geom_sf(data = Sites_Target, aes(color = Site_Type), size = 5) +
  scale_color_manual(values = c("Donor" = "red", "Recipient" = "black"),
                     name = "Site Type") +
  coord_sf(
    xlim = c(-77.58, -76.0),
    ylim = c(38.70, 39.70),
    crs = 4326,
    expand = FALSE
  ) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"), width = unit(1.2, "cm")
  ) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   line_width = 0.8, text_cex = 1.0) +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold")
  )

print(Main_Map)

## Build the Maryland inset map 
Maryland_Inset <- ggplot() +
  geom_sf(data = Maryland_Counties, fill = "gray", color = "white", linewidth = 0.3) +
  geom_sf(data = Study_Counties, fill = "#F98E3B", color = "black", linewidth = 0.4) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )
print((Maryland_Inset))

## Build a United States Inset Map
install.packages("USAboundaries")
library(USAboundaries)

# Get US states
us_states <- us_states()
print(us_states)

# Remove Alaska (and Hawaii if you want)
us_mainland <- us_states[!us_states$state_name %in% c("Alaska", "Puerto Rico", "Hawaii"), ]

# Extract Maryland
maryland <- us_mainland[us_mainland$state_name == "Maryland", ]

US_Inset <- ggplot() +
  geom_sf(data = us_mainland, fill = "grey90", color = "white") +
  geom_sf(data = maryland, fill = "#F98E3B", color = "black", size = 0.6) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
print(US_Inset)

## Combine Maps 
library(cowplot)
library(ggspatial)

Final_Map <- ggdraw() +
  draw_plot(Main_Map) +
  draw_plot(Maryland_Inset, 
            x = 0.07, y = 0.67,   # adjust inset position
            width = 0.30, height = 0.40) + 
  draw_plot(US_Inset, 
            x = 0.05, y = -0.05,   # adjust inset position
            width = 0.30, height = 0.40)

print(Final_Map)
