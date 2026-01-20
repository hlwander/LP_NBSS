#script to pull supplemental phytoplankton colony lengths for n = 17 filamentous autotrophs
#looking for published phyto size data for filamentous (and 1 tiny) phytos

pacman::p_load(dplyr, tidyverse)

bea_lengths <- read.csv("data/phyto.trait.length.merged.csv", row.names=1) |>
  mutate(gald_um = as.numeric(Bea.s.value.in.um), #select Bea's value for colonial species; greatest axial linear dimension
         target_taxon = sub(" .*", "", genus.species)) |>  
  dplyr::select(target_taxon, gald_um) |>
  group_by(target_taxon) |>
  summarise(gald_avg_um = mean(gald_um, na.rm=T), .groups = "drop") |>
  bind_rows(tibble(target_taxon = "Lyngbya", 
                   gald_avg_um = 200)) 

phyto_list <- bea_lengths |>
  filter(!is.nan(gald_avg_um)) |>
  dplyr::select(target_taxon) #listing the genera that Bea flagged bc they are either colonial or v tiny

#read in the density/biomass file; units: biomass = mg/m3, density = #/L, GALD = um
phyto_dens_biom <- read.csv("data/LP_allphyto.full_01Oct2025.csv", row.names=1) |>
  mutate(target_taxon = sub("\\s.*", "", Basionym)) |>
  rename(biomass_mg_m3 = biomass,
         density_nopL = density,
         gald_um = gald) |>
  mutate(target_taxon = ifelse(target_taxon == "A.", "Anabaena", 
                               ifelse(target_taxon == "M.", "Microcystis", target_taxon))) |>
  filter(!(biomass_mg_m3 == 0 & density_nopL == 0)) |>
  arrange(lake_id) |>
  dplyr::select(target_taxon, lake_id, biomass_mg_m3, density_nopL, gald_um) |>
  distinct() |>
  mutate(target_taxon = case_when(
    target_taxon == "Cyrptomonas" ~ "Cryptomonas",
    target_taxon == "Crucegenia" ~ "Crucigenia",
    target_taxon == "Scenesdesmus" ~ "Scenedesmus",
    target_taxon == "Mallomomas" ~ "Mallomonas",
    target_taxon == "Chamydomonas" ~ "Chlamydomonas",
    target_taxon == "Pedridium" ~ "Peridinium",
    target_taxon == "Perdinium" ~ "Peridinium",
    target_taxon == "Perididium" ~ "Peridinium",
    target_taxon == "Cosmariun" ~ "Cosmarium",
    target_taxon == "Trachemolomas" ~ "Trachelomonas",
    target_taxon == "Trachelomomas" ~ "Trachelomonas",
    target_taxon ==  "Anabaena" ~ "Dolichospermum",
    target_taxon == "Micratinium" ~ "Micractinium",
    target_taxon == "Planktolyngya" ~ "Planktolyngbya",
    target_taxon == "Planktonema" ~ "Planctonema",
    target_taxon == "Schroderia" ~ "Schroederia",
    target_taxon == "Haptophyte" ~ "Chrysochromulina", #id'ed as "Haptophyte (Erkenia/Chrysochromulina)", but Erkenia does not seem to be a fully independent genus according to google
    target_taxon == "Limnothrix/Jaaginema" ~ "Limnothrix", #lumping this taxon in with Limnothrix because Jaaginema is not present in dataset
    target_taxon == "Unidentified" ~ "Other",   
    TRUE ~ target_taxon  # keep everything else as is
  )) #keeping Cyclotella/Stephanodiscus as a distinct taxon bc both genera are already present in dataset Limnothrix/Jaaginema

phyto_dens_biom_genera <- phyto_dens_biom |>
  group_by(target_taxon, lake_id) |> 
  summarise(biomass_mg_m3 = sum(biomass_mg_m3, na.rm = TRUE), 
            density_nopL = sum(density_nopL, na.rm = TRUE), 
            gald_um = mean(gald_um, na.rm = TRUE),
            .groups = "drop") 

#caclulate % contribution per lake
phyto_pct_per_lake <- phyto_dens_biom_genera |>
  group_by(lake_id) |>
  mutate(total_biomass_lake = sum(biomass_mg_m3, na.rm = TRUE),
         pct_biomass = if_else(total_biomass_lake > 0,
                               100 * biomass_mg_m3 / total_biomass_lake, 0)) |>
  ungroup() |>
  group_by(target_taxon) |>
  summarise(n_lakes_present = sum(biomass_mg_m3 > 1e-6, na.rm = TRUE),
            max_pct_across_lakes = max(pct_biomass, na.rm = TRUE),
            mean_pct_across_lakes = mean(pct_biomass, na.rm = TRUE),
            total_biomass_across_dataset = sum(biomass_mg_m3, na.rm = TRUE),
            .groups = "drop")

#identify taxa with max_pct < 1% AND occur in fewer than 5 lakes
rare_taxa <- phyto_pct_per_lake |>
  filter(max_pct_across_lakes < 1 |
           n_lakes_present < 5) |>
  arrange(n_lakes_present, max_pct_across_lakes)

#remove rare taxa
final_phyto_dens_biom <- phyto_dens_biom_genera |>
  filter(!target_taxon %in% rare_taxa$target_taxon)

#now combine dfs
phytos_all <- final_phyto_dens_biom |>
  left_join(phyto_list, by = "target_taxon") |>
  mutate(lake_id = str_replace(lake_id, "^id", ""),   # remove "id" at start
         lake_id = str_replace(lake_id, "\\.", "-"),        # replace "." with "-"
         lake_id = str_extract(lake_id, "^[0-9]{2}-[0-9]{3}"))   # replace "." with "-"

#only select the phytos that will be in the analysis
phyto_list <- phyto_list |> #n = 17
  filter(target_taxon %in% phytos_all$target_taxon)

#------------------------------------------------------------------------------
#supplemental phyto lengths

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
#write.csv(phyto_size, "data/size/published_phytoplankton_colony_lengths.csv", row.names=F)


#delete later
#quick check against Bea's assigned values
phyto_size <- phyto_size |>
  left_join(bea_lengths, by = "target_taxon") |>
  rename(bea_length = gald_avg_um)

