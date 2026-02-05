#script to pull supplemental phytoplankton colony lengths for n = 17 filamentous autotrophs
#looking for published phyto size data for filamentous (and 1 tiny) phytos

pacman::p_load(dplyr, tidyverse)

#read in list of colonial autotrophs that need supplemental lengths
phyto_list <- read.csv("data/colonial_autotrophs.csv")

#Phytoplankton morpho-functional trait dataset from French water-bodies
size1 <- read.csv("data/supp_lengths/FRENCH_PHYTOPLANKTON_TRAITS.csv", sep=";") |>
  dplyr::select(Taxa_Name, Genus, Min_Length, max_Lenght) |>
  filter(Genus %in% phyto_list$target_taxon) |>
  group_by(Genus) |>
  summarise(max_length = mean(max_Lenght), #length of the cell or colony
            min_length = mean(Min_Length)) |>
  rename(target_taxon = Genus) # 14/17
#metadata: https://www.nature.com/articles/s41597-021-00814-0/tables/3

# from Rimet and Durat 2018 of >1200 taxa commonly observed in temperate lakes
size2 <- read.csv("data/supp_lengths/Rimet_and_Druart_2018.csv", stringsAsFactors = FALSE) |>
  filter(!Genus...species.name %in% "Cyst of péridinien") |>
  mutate(Genus = trimws(Genus),
         Genus = if_else(Genus == "", word(Genus...species.name, 1), Genus)) |>
  filter(Genus %in% phyto_list$target_taxon) |>
  group_by(Genus) |>
  summarise(mean_length2 = mean(coalesce(Colony.length..takes.into.account.mucilage..µm, 
                                         Cell.length.µm), na.rm = TRUE),.groups = "drop") |>
  rename(target_taxon = Genus) # 16/17
#https://www.limnology-journal.org/articles/limn/full_html/2018/01/limn170074/limn170074.html
#the paper says that the lengths were arbitrarily assigned 100 um to use for biovolume, but there are plenty of lengths that are >100 um

#LTER phytoplankton epilimnetic data 1984-2015 --> Phycotech
size3 <- read.csv("data/supp_lengths/cascade_phytoplankton_v0.1.csv_upload.csv") |>
  select(lakename, genus, description, gal_dimension) |>
  filter(genus %in% phyto_list$target_taxon) |>
  group_by(genus) |>
  summarise(mean_length3 = mean(gal_dimension)) |>
  rename(target_taxon = genus) # 9/17
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.353.5

# USGS Finger Lake phytos 2019-2020
size4 <-read.csv("data/supp_lengths/FingerLakesHABs_Phytoplankton_Data.csv", skip=2) |>
  select(GENUS, AVERAGE_GALD.UM) |>
  filter(!GENUS %in% "--", #first remove rows that do not have corresponding genus
         GENUS %in% phyto_list$target_taxon) |>
  group_by(GENUS) |>
  summarise(mean_length4 = mean(AVERAGE_GALD.UM)) |> 
  rename(target_taxon = GENUS) # 11/17
#https://www.sciencebase.gov/catalog/item/60ae4e58d34e4043c8538b8a
# note from readme --> "AVERAGE_GALD.UM, Greatest Axial Linear Dimension; data are calculated raw values and are not rounded to USGS significant figures, in micrometers."

#combine dfs --> all 17 genera have values from at least 2 dfs to average
phyto_size <- phyto_list |>
  left_join(size1, by = "target_taxon") |>
  left_join(size2, by = "target_taxon") |>
  left_join(size3, by = "target_taxon") |>
  left_join(size4, by = "target_taxon") |>
  mutate(midpoint = (min_length + max_length) / 2,
         avg_length_um = rowMeans(across(c(
           mean_length2, mean_length3, 
           mean_length4, midpoint)), na.rm = TRUE)) |>
  select(target_taxon, avg_length_um) 
#write.csv(phyto_size, "data/supp_lengths/collated_mean_phyto_colony_lengths.csv", row.names=F)
