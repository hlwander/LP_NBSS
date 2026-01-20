# calculating size spectra slopes using functional feeding groups

pacman::p_load(tidyverse, dplyr, ggplot2, ggmap,
               rnaturalearth, rnaturalearthdata, 
               ARTool, ggpubr, mgcv, ggpp)

#tidyverse, broom, purrr, ggplot2, stringr, mice,
#car, ARTool, corrplot, magrittr, vegan, mgcv, dplyr,
#rnaturalearth, rnaturalearthdata, ggmap, multcompView, 
#multcomp, emmeans, ggpubr, scales, ggridges, forcats,
#viridis, gratia, ggpp

# calculate trophic state using chla (file modified from CP where NAs were replaced using miss forest or ecozonal means)
pred <- read.csv("data/lake_chla_ecozone.csv") |>
  rename(lake_id = lakepulse_id) |>
  mutate(trophic_state = case_when( is.na(chla_day) ~ NA_character_, 
                                    chla_day < 2 ~ "Oligotrophic",
                                    chla_day < 7 ~ "Mesotrophic",
                                    chla_day < 30 ~ "Eutrophic",
                                    TRUE ~ "Hypereutrophic")) |>
  filter(!is.na(trophic_state)) |> #removing lakes with NA chla values
  mutate(trophic_state = if_else(trophic_state == "Hypereutrophic", 
                                 "Eutrophic", trophic_state),
         trophic_state = factor(trophic_state, levels = c(
           "Oligotrophic", "Mesotrophic", "Eutrophic"))) 

# read in supplemental phyto length data
phyto_supp_lengths <- read.csv("data/supp_lengths/published_phytoplankton_colony_lengths.csv") |>
  rename(supp_length_um = avg_length_um)

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
  )) |>  #keeping Cyclotella/Stephanodiscus as a distinct taxon bc both genera are already present in dataset Limnothrix/Jaaginema
  dplyr::select(target_taxon, lake_id, biomass_mg_m3, density_nopL, gald_um)

#calculate % contribution per lake
phyto_pct_per_lake <- phyto_dens_biom |>
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
final_phyto_dens_biom <- phyto_dens_biom |>
  filter(!target_taxon %in% rare_taxa$target_taxon)

#bring in zoop size data 
zoops_raw <- read.csv("data/raw_zoops.csv") |> #was called all_rawdata.csv
  rename(lake_id = ID_lakepulse,
         density_NopL = X....L,
         biomass_ugL = species.biomass..µg.d.w..L.) |>
  filter(biomass_ugL > 0 & !is.na(biomass_ugL),
         X.individuals.counted >= 100, #following CP methods
         !comments %in% "Unidentified sample either 09-420 or 11-343.") |>
  mutate(target_taxon = str_to_sentence(genus),
         target_taxon = ifelse(species %in% c("sicilis", "sicilis copepodid"),
                               "Leptodiaptomus sicilis", 
                               target_taxon)) |>
  pivot_longer(cols = c(D1:D10), names_to = "rep",
               values_to = "length_mm") |>
  mutate(length_um = length_mm * 1000,
         lake_id = trimws(lake_id)) |> 
  dplyr::select(lake_id, target_taxon, density_NopL, biomass_ugL, length_um) |>
  distinct()
#not summarizing across genera bc we want the length variability for size spectra

#calculate % contribution per lake
zoop_pct_per_lake <- zoops_raw |>
  group_by(lake_id) |>
  mutate(total_biomass_lake = sum(biomass_ugL, na.rm = TRUE),
         pct_biomass = if_else(total_biomass_lake > 0,
                               100 * biomass_ugL / total_biomass_lake, 0)) |>
  ungroup() |>
  group_by(target_taxon) |>
  summarise(n_lakes_present = sum(biomass_ugL > 1e-6, na.rm = TRUE),
            max_pct_across_lakes = max(pct_biomass, na.rm = TRUE),
            mean_pct_across_lakes = mean(pct_biomass, na.rm = TRUE),
            total_biomass_across_dataset = sum(biomass_ugL, na.rm = TRUE),
            .groups = "drop")

#identify taxa with max_pct < 1% AND occur in fewer than 5 lakes
zoop_rare_taxa <- zoop_pct_per_lake |>
  filter(max_pct_across_lakes < 1 |
           n_lakes_present < 5) |>
  arrange(n_lakes_present, max_pct_across_lakes)
#none to drop

#remove rare taxa and rename cols
final_zoop_df <- zoops_raw |>
  filter(!target_taxon %in% zoop_rare_taxa$target_taxon) |>
  rename(biomass = biomass_ugL, #dropping units for merge w/ phytos
         density_nopL = density_NopL) 

#now combine dfs
phytos_all <- final_phyto_dens_biom |>
  left_join(phyto_supp_lengths, by = "target_taxon") |>
  mutate(lake_id = str_replace(lake_id, "^id", ""),   # remove "id" at start
         lake_id = str_replace(lake_id, "\\.", "-"),        # replace "." with "-"
         lake_id = str_extract(lake_id, "^[0-9]{2}-[0-9]{3}")) |>  # replace "." with "-"
  mutate(supp_length_um = ifelse(is.na(supp_length_um), gald_um, supp_length_um)) |>
  dplyr::select(-gald_um) |>
  rename(length_um = supp_length_um,
         biomass = biomass_mg_m3) #dropping units for merge w/ zoops

#merge with the zoop df
all_plankton <- bind_rows(phytos_all, final_zoop_df)

#log the gald col
all_plankton <- all_plankton |> 
  mutate(log_length = log10(as.numeric(length_um)))

#define bin breaks
breaks <- seq(floor(min(all_plankton$log_length, na.rm = TRUE) / 0.1) * 0.1,
              ceiling(max(all_plankton$log_length, na.rm = TRUE) / 0.1) * 0.1, by = 0.1)

# pull in phyto nutritional strategies
strategy <- read.csv("data/nanoplanktonnutritionstrategiesdb_v2_fbae035.csv") |>
  mutate(target_taxon = Genus) |> 
  bind_rows(tibble(target_taxon = c("Cyanobacteria", "Diatoms", "Haptophyte", "Eutetramorus"), #bc all cyaos and diatoms are autotrophs and haptophytes are mixotrophs
                   Final_Nutrition_Strategy = c("Autotroph", "Autotroph", "Mixotroph", "Autotroph"))) |>
  rename(trophic_group = Final_Nutrition_Strategy)
#adding Eutetramorus because Leong and Chang 2023 say that E. planctonicus is autotrophic (and cite	D’Alessandro et al. (2018))

# pull in zoop functional groups
zoop_function <- read.csv("data/Functionaltraitsmatrix.csv", sep = ";") |>
  mutate(target_taxon = str_to_sentence(word(Species.name, 1, sep = "\\."))) |>
  rename(trophic_group = trophic.group) |>
  dplyr::select(target_taxon, trophic_group) |>
  filter(target_taxon %in% final_zoop_df$target_taxon) |>
  mutate(trophic_group = if_else(trophic_group %in% c(
    "Herbivore", "Omnivore/Herbivore"), "Herbivore", "Non-herbivore")) |>
  mutate(target_taxon = ifelse(trophic_group == "Non-herbivore" & 
                                 target_taxon == "Leptodiaptomus",
                               "Leptodiaptomus sicilis", 
                               target_taxon)) |>
  distinct()
#herbivores vs. everyone else
#Lepto is both omnivore (siicilis) and herbivore

#join plankton functional groups
fun_groups <- bind_rows(strategy, zoop_function) |>
  dplyr::select(target_taxon, trophic_group)

# join functional groups w/ plankton data
# some taxa have multiple sizes per lake
plankton_with_strategy <- all_plankton |>
  left_join(fun_groups,  by = "target_taxon") |>
  mutate(trophic_group = ifelse(target_taxon == "Chydorus ", "Non-herbivore", 
                                ifelse(target_taxon == "Cyclotella/Stephanodiscus",
                                       "Mixotroph", trophic_group)))
#write.csv(plankton_with_strategy, "output/plankton_final.csv", row.names = FALSE)
#some species of Stephanodiscus are mixotrophic so grouping as such

#bin taxa and sum biomass per bin per lake
spectra_lake_bins <- plankton_with_strategy |>
  mutate(size_bin = cut(log_length, breaks = breaks, include.lowest = TRUE),
         bin_id = as.integer(size_bin),
         log_lower = breaks[bin_id],
         log_upper = breaks[bin_id + 1],
         size_lower = 10^log_lower, #linear units
         size_upper = 10^log_upper, #linear units
         bin_width = size_upper - size_lower,
         bin_mid = (log_lower + log_upper) / 2) |>
  filter(!is.na(biomass)) |>
  group_by(lake_id, trophic_group, bin_mid, size_lower, size_upper, bin_width) |>
  summarise(total_biomass = sum(biomass, na.rm = TRUE), .groups = "drop") |>
  mutate(nbss = total_biomass / bin_width,
         log_nbss = ifelse(nbss > 0, log10(nbss), NA_real_),
         log_size = bin_mid) |>
  ungroup() |>
  mutate(trophic_group = factor(trophic_group, levels = c(
    "Autotroph", "Mixotroph", "Herbivore", "Non-herbivore")),
    plankton = ifelse(trophic_group %in% c("Mixotroph","Autotroph"),
                      "Phytoplankton", "Zooplankton"))

#calculate NBSS as slope of log NBSS vs log size bin
fits <- spectra_lake_bins |>
  filter(!is.na(log_nbss), !is.na(log_size))  |>
  group_by(lake_id, trophic_group) |>  
  summarise(n_bins = n(), model = list(
    if (n() >= 5) {
      tryCatch(lm(log_nbss ~ log_size, data = cur_data()), error = function(e) NULL)
    } else NULL), .groups = "drop") |>
  mutate(model_class = map(model, ~ {if (is.null(.x)) {"NULL"
  } else {paste(class(.x), collapse = "|")}
  }) |> unlist(),
  ok = map_lgl(model, ~ !is.null(.x) && inherits(.x, "lm")))

# Keep only successful fits
good <- fits |> filter(ok) |> #keeping at least 5 size class bins
  mutate(tidy = map(model, broom::tidy), # 599/659 lakes 
         glance = map(model, broom::glance))

# Extract tidy coefficients and glance stats, then pivot coefficients to columns
coeffs_wide <- good |>
  mutate(tidy = map(model, broom::tidy),
         glance = map(model, broom::glance)) |>
  dplyr::select(lake_id, trophic_group, n_bins, tidy, glance) |>
  unnest(tidy) |>
  pivot_wider(id_cols = c(lake_id, trophic_group, n_bins),
              names_from = term,
              values_from = c(estimate, std.error, p.value),
              names_glue = "{.value}__{term}")

#rename cols
names(coeffs_wide) <- names(coeffs_wide) %>%
  gsub("^estimate__\\(Intercept\\)$", "intercept", .) %>%
  gsub("^std.error__\\(Intercept\\)$", "intercept_se", .) %>%
  gsub("^p.value__\\(Intercept\\)$", "intercept_p", .) %>%
  gsub("^estimate__log_size$", "slope", .) %>%
  gsub("^std.error__log_size$", "slope_se", .) %>%
  gsub("^p.value__log_size$", "slope_p", .)

# Unpack glance stats into columns and join with coeffs
glance_wide <- good |>
  dplyr::select(lake_id, trophic_group, glance) |>
  unnest_wider(glance)

#this drops 11 lakes that do not have trophic state classifications
fits_summary <- coeffs_wide |>
  left_join(glance_wide, by = c("lake_id", "trophic_group")) |>
  left_join(pred |> dplyr::select(lake_id, trophic_state, ecozone), 
            by = "lake_id") |>
  arrange(trophic_group, lake_id) |>
  mutate(plankton = ifelse(trophic_group %in% c("Mixotroph","Autotroph"),
                           "Phytoplankton", "Zooplankton"),
         trophic_state = if_else(trophic_state == "Hypereutrophic", 
                                 "Eutrophic", trophic_state),
         trophic_state = factor(trophic_state, levels = c("Oligotrophic","Mesotrophic","Eutrophic"))) |>
  filter(!is.na(trophic_state)) 

#----------------------------------------------------------------------------#
# manuscript figs

#split into phyto and zoop dfs
fits_phyto <- fits_summary |> filter(plankton == "Phytoplankton")
fits_zoop  <- fits_summary |> filter(plankton == "Zooplankton")

#two ART ANOVAs for phytos vs. zoops (no interaction bc missing mixotrophs in eutrophic lakes)
mod_art_phyto_fg <- art(slope ~ trophic_group , data = fits_phyto)
mod_art_phyto_ts <- art(slope ~ trophic_state , data = fits_phyto)

anova(mod_art_phyto_fg) #autotroph and mixotroph slopes are significantly different
anova(mod_art_phyto_ts) #slopes are sig different across trophic states

#post-hoc - main effects
art.con(mod_art_phyto_ts, "trophic_state")
#oligo and meso lakes have sig lower slope than eutrophic

mod_art_zoop_fg <- art(slope ~ trophic_group, data = fits_zoop)
mod_art_zoop_ts <- art(slope ~ trophic_state, data = fits_zoop)

anova(mod_art_zoop_fg) #herbivore and non-herbivore slopes are significantly different
anova(mod_art_zoop_ts) #not significant

# Combine phyto and zoop
fits_all <- bind_rows(
  fits_phyto |> mutate(plankton = "Phytoplankton"),
  fits_zoop  |> mutate(plankton = "Zooplankton"))

# Summarize slopes for Figure 3
slope_summary_ts <- fits_all |>
  group_by(plankton, trophic_state) |>
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    se = sd(slope, na.rm = TRUE) / sqrt(n()),
    n = n(), .groups = "drop") |>
  mutate(Letter = case_when(
    plankton == "Phytoplankton" & trophic_state %in% c(
      "Oligotrophic", "Mesotrophic") ~ "a",
    plankton == "Phytoplankton" &  trophic_state %in% c(
      "Eutrophic") ~ "b",
    TRUE ~ as.character(NA))) |>
  mutate(trophic_state = factor(trophic_state, levels = c(
    "Oligotrophic","Mesotrophic","Eutrophic")))

p_ts <- ggplot(slope_summary_ts, 
               aes(x = trophic_state, y = mean_slope, color = trophic_state)) +
  geom_point(position = position_dodge(width = 0.6), 
             size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_slope - se, ymax = mean_slope + se),
                position = position_dodge(width = 0.6), width = 0.2, 
                color = "black") + ylim(c(-0.7,1.0)) +
  geom_text(aes(y = mean_slope + se + 0.03, label = Letter),
            position = position_dodge(width = 0.6),
            vjust = 0, color = "black", size = 3) +
  facet_wrap(~plankton) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  theme_minimal() + xlab("") + ylab("Mean slope") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(),
        text = element_text(size=6),
        strip.text = element_text(face = "bold")) 

slope_summary_fg <- fits_all |>
  group_by(plankton, trophic_group) |>
  summarise(mean_slope = mean(slope, na.rm = TRUE),
            se = sd(slope, na.rm = TRUE) / sqrt(n()),
            n = n(), .groups = "drop") |>
  mutate(Letter = case_when(trophic_group == "Autotroph" ~ "a",
                            trophic_group == "Mixotroph" ~ "b",
                            trophic_group == "Herbivore" ~ "B",
                            trophic_group == "Non-herbivore" ~ "A",
                            TRUE ~ as.character(NA))) |>
  mutate(trophic_group = factor(trophic_group, levels = c(
    "Autotroph","Mixotroph","Herbivore","Non-herbivore")))

p_fg <- ggplot(slope_summary_fg, 
               aes(x = trophic_group, y = mean_slope, color = trophic_group)) +
  geom_point(position = position_dodge(width = 0.6), 
             size = 3, show.legend = FALSE) + ylim(c(-0.7, 1.4)) +
  geom_errorbar(aes(ymin = mean_slope - se, ymax = mean_slope + se),
                position = position_dodge(width = 0.6), 
                width = 0.2, color = "black") +
  geom_text(aes(y = mean_slope + se + 0.1, label = Letter),
            position = position_dodge(width = 0.6),
            vjust = 0, color = "black", size = 3) +
  theme_minimal() + xlab("") + ylab("Mean slope") +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  scale_color_manual(values = c(
    "Autotroph" = "#3E6E66",
    "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C",
    "Non-herbivore" = "#F68A4D")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=6),
        axis.ticks = element_line()) 

ggarrange(p_ts, p_fg, ncol = 1, nrow = 2,       
          heights = c(1.4, 1),  
          widths = c(2,1), common.legend = FALSE)
#ggsave("figs/plankton_sig_slopes.jpg", width = 5, height = 3)

# slope across functional groups (Figure 2)
spectra_summary <- spectra_lake_bins |>
  left_join(pred |> dplyr::select(trophic_state, lake_id), by = "lake_id") |>
  group_by(trophic_group, trophic_state, log_size) |> 
  summarise(mean_log_nbss = mean(log_nbss, na.rm = TRUE),
            sd_log_nbss = sd(log_nbss, na.rm = TRUE),
            n = n()) |>
  mutate(se_log_nbss = sd_log_nbss / sqrt(n),
         trophic_state = if_else(trophic_state == "Hypereutrophic", 
                                 "Eutrophic", trophic_state),
         trophic_state = factor(trophic_state, levels = c(
           "Oligotrophic","Mesotrophic","Eutrophic"))) |>
  filter(!is.nan(mean_log_nbss),
         !is.na(trophic_group),
         !is.na(trophic_state))

ggplot(spectra_summary, aes(x = log_size, y = mean_log_nbss, 
                            color = trophic_group)) +
  geom_point(show.legend = T) +
  geom_smooth(method = "lm", se = FALSE,  show.legend = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, 
              color = "grey50", linetype = "dashed") +
  scale_x_continuous("Log10 mean size (µm)") +
  scale_y_continuous("Log10 NBSS") +
  scale_color_manual(values = c(
    "Autotroph" = "#3E6E66",
    "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C",
    "Non-herbivore" = "#F68A4D")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggsave("figs/slopes_by_group.jpg", width = 5, height = 4)

#slopes faceted by trophic state (Figure S3)
ggplot(spectra_summary, aes(x = log_size, y = mean_log_nbss, color = trophic_group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +      
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), 
              se = FALSE, color = "black") +
  facet_wrap(~ trophic_state, nrow=1) +                                         
  scale_x_continuous("Log10 mean size (µm)") +
  scale_y_continuous("Log10 NBSS") +
  scale_color_manual(values = c(
    "Autotroph" = "#3E6E66",
    "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C",
    "Non-herbivore" = "#F68A4D")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggsave("figs/ssa_functional_groups_by_ts.jpg", width = 6, height = 5)

#------------------
#manuscript tables (Table S2)
group_slopes <- spectra_summary |>
  filter(!is.na(trophic_group), !is.na(log_size), !is.na(mean_log_nbss)) |>
  group_by(trophic_group) |>
  nest() |>
  mutate(n = map_int(data, nrow),
         model = map(data, ~ lm(mean_log_nbss ~ log_size, data = .x)),
         tidy = map(model, broom::tidy),
         glance = map(model, broom::glance)) |>
  unnest(tidy) |>
  filter(term == "log_size") |>
  transmute(trophic_group,
            estimate = estimate,
            std.error = std.error,
            t_value = statistic,
            p_value = p.value,
            n, r.squared = map_dbl(glance, "r.squared")) |>
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         p_value = signif(p_value, 3),
         t_value = round(t_value, 3),
         r.squared = round(r.squared, 3)) |>
  ungroup() |>
  mutate(p_adj = signif(p.adjust(p_value, method = "BH"),3)) |>
  dplyr::select(-p_value)
#write.csv(group_slopes, "output/fg_slopes.csv", row.names = FALSE)

#fg and ts slopes (Table S3)
group_by_state_slopes <- spectra_summary |>
  filter(!is.na(trophic_state), !is.na(trophic_group)) |>
  group_by(trophic_state, trophic_group) |>
  nest() |>
  mutate(n = map_int(data, nrow),
         model = map(data, ~ lm(mean_log_nbss ~ log_size, data = .x)),
         tidy = map(model, broom::tidy),
         glance = map(model, broom::glance)) |>
  unnest(tidy) |>
  filter(term == "log_size") |>
  transmute(trophic_state,
            trophic_group,
            estimate = estimate,
            std.error = std.error,
            t_value = statistic,
            p_value = p.value,
            n, r.squared = map_dbl(glance, "r.squared")) |>
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         p_value = signif(p_value, 3),
         t_value = round(t_value, 3),
         r.squared = round(r.squared, 3)) |>
  ungroup() |>
  mutate(p_adj = signif(p.adjust(p_value, method = "BH"),3)) |>
  dplyr::select(-p_value)
#write.csv(group_by_state_slopes, "output/fg_ts_slopes.csv", row.names = FALSE)

# overall quadratic 
quad_all <- lm(mean_log_nbss ~ log_size + I(log_size^2), data = spectra_summary)
summary(quad_all)        
broom::tidy(quad_all)     
broom::glance(quad_all)   

#Figure S4
lake_slope_summary <- fits_collapsed |>
  group_by(trophic_state, trophic_group) |>
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),        
    sd_slope   = sd(slope, na.rm = TRUE),          
    n_lakes    = n(),                              
    se_slope   = sd_slope / sqrt(n_lakes),         
    .groups = "drop")

ggplot(lake_slope_summary, aes(x = trophic_group, 
                               y = mean_slope, color = trophic_state)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope),
                position = position_dodge(width = 0.6), width = 0.2) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  labs(y = "Mean NBSS slope ± SE", x = "")
#ggsave("figs/slope_ts_fg.jpg", width = 6, height = 5)

#raw plankton biomass by group
biomass_abs <- all_plankton |>
  left_join(fun_groups, by = "target_taxon") |>
  left_join(pred |> dplyr::select(lake_id, trophic_state), by = "lake_id") |>
  group_by(lake_id, trophic_group, trophic_state) |>
  summarise(total_biomass = sum(biomass, na.rm = TRUE),
            mean_size = mean(length_um, na.rm = TRUE), #gald_um
            .groups = "drop") |>
  mutate(trophic_group = factor(trophic_group, levels = c(
    "Autotroph", "Mixotroph", "Herbivore", "Non-herbivore")))

# compute group means per trophic_state + trophic_group
dat <- biomass_abs |> 
  filter(!is.na(trophic_group),
         !is.na(trophic_state)) |>
  mutate(log_size = log10(mean_size),
         plankton = ifelse(trophic_group %in% c("Herbivore","Non-herbivore"),
                           "zooplankton","phytoplankton"))

means <- dat |>
  group_by(trophic_state, trophic_group) |> 
  summarize(mu = median(mean_size, na.rm = TRUE), .groups = "drop") |>
  mutate(plankton = ifelse(trophic_group %in% c("Herbivore","Non-herbivore"),
                           "zooplankton","phytoplankton"))

#histogram of plankton size across trophic states (Figure 1)
ggplot(dat |> filter(), aes(x = mean_size, fill = trophic_group)) +
  geom_density(alpha = 0.8) +
  geom_vline(data = means, aes(xintercept = mu, color = trophic_group, group = plankton), 
             linetype = "dashed", size = 0.6, show.legend = FALSE) +
  facet_wrap(~trophic_state + plankton, nrow = 3, scales = "free",
             labeller = labeller(trophic_state = label_value, 
                                 plankton = function(x) "")) +
  theme_minimal(base_size = 10) +
  scale_fill_manual(values = c(
    "Autotroph" = "#3E6E66",
    "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C",
    "Non-herbivore" = "#F68A4D")) +
  scale_color_manual(values = c(
    "Autotroph" = "#3E6E66",
    "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C",
    "Non-herbivore" = "#F68A4D")) +
  labs(x = "Mean size (µm)", y = "Density", fill = "") +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggsave("figs/ts_size_density_plots.jpg", width = 6, height = 5)

#values for results text
median(dat$mean_size[dat$trophic_group=="Autotroph"])
median(dat$mean_size[dat$trophic_group=="Mixotroph"])
median(dat$mean_size[dat$trophic_group=="Herbivore"])
median(dat$mean_size[dat$trophic_group=="Non-herbivore"])

(median(dat$mean_size[dat$trophic_group=="Autotroph"]) - 
    median(dat$mean_size[dat$trophic_group=="Mixotroph"])) /
  median(dat$mean_size[dat$trophic_group=="Mixotroph"]) *100

median(dat$mean_size[dat$trophic_group=="Autotroph" & 
                       dat$trophic_state %in% "Oligotrophic"])
median(dat$mean_size[dat$trophic_group=="Autotroph" & 
                       dat$trophic_state %in% "Mesotrophic"])
median(dat$mean_size[dat$trophic_group=="Autotroph" & 
                       dat$trophic_state %in% "Eutrophic"])

median(dat$mean_size[dat$trophic_group=="Mixotroph" & 
                       dat$trophic_state %in% "Oligotrophic"])
median(dat$mean_size[dat$trophic_group=="Mixotroph" & 
                       dat$trophic_state %in% "Mesotrophic"])
median(dat$mean_size[dat$trophic_group=="Mixotroph" & 
                       dat$trophic_state %in% "Eutrophic"])

median(dat$mean_size[dat$trophic_group=="Herbivore" & 
                       dat$trophic_state %in% "Oligotrophic"])
median(dat$mean_size[dat$trophic_group=="Herbivore" & 
                       dat$trophic_state %in% "Mesotrophic"])
median(dat$mean_size[dat$trophic_group=="Herbivore" & 
                       dat$trophic_state %in% "Eutrophic"])

median(dat$mean_size[dat$trophic_group=="Non-herbivore" & 
                       dat$trophic_state %in% "Oligotrophic"])
median(dat$mean_size[dat$trophic_group=="Non-herbivore" & 
                       dat$trophic_state %in% "Mesotrophic"])
median(dat$mean_size[dat$trophic_group=="Non-herbivore" & 
                       dat$trophic_state %in% "Eutrophic"])

# Compute medians (slope and chla) across trophic state and functional groups
median_points <- fits_collapsed |>
  left_join(pred |> dplyr::select(lake_id, chla_day), by = "lake_id") |>
  group_by(trophic_group, trophic_state) |>
  summarise(median_slope = median(slope, na.rm = TRUE),
            median_logchla = median(log10(chla_day), na.rm = TRUE),
            .groups = "drop")

#chla vs slope (Figure 5)
fits_collapsed |>
  left_join(pred |> dplyr::select(lake_id, chla_day), by = "lake_id") |>
  mutate(log_chla = log10(chla_day),
         slope_dir = ifelse(slope >= 0, "positive", "negative")) |>
  ggplot(aes(x = slope, y = log_chla, color = trophic_state, fill = trophic_state)) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.6) +  
  geom_point(data = median_points,
             aes(x = median_slope, y = median_logchla, fill = trophic_state),
             shape = 21, size = 2, stroke = 0.8, color = "black",inherit.aes = FALSE) +
  facet_wrap(~trophic_group, ncol = 4) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  scale_fill_manual(values = c("Oligotrophic" = "#21908CFF", 
                               "Mesotrophic" = "#FDE725FF",
                               "Eutrophic" = "#440154FF")) +
  labs(x = "Slope", y = "Log10 Chla (µg/L)", color = "", fill = "") +
  theme_minimal()  +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggsave("figs/chla_vs_slope.jpg", width = 5, height = 3)

#add log chla to dat
dat <- dat |> 
  left_join(pred |> dplyr::select(lake_id, chla_day), by = "lake_id") |> 
  mutate(log_chla = log10(chla_day))

#fit a GAM for each functional group
models <- dat |>
  group_by(trophic_group) |>
  group_map(~ gam(log_size ~ s(log_chla, k = 6), data = .x))

#name each model
names(models) <- unique(dat$trophic_group)

#rearrange panel order 
dat <- dat |>
  mutate(trophic_group = factor(trophic_group, levels = c(
    "Herbivore","Non-herbivore","Autotroph","Mixotroph"))) |>
  mutate(facet_row = case_when(
    trophic_group %in% c("Herbivore", "Non-herbivore") ~ "Row1",
    trophic_group %in% c("Autotroph", "Mixotroph") ~ "Row2"),
    facet_col = case_when(
      trophic_group %in% c("Herbivore", "Autotroph") ~ "Left",
      trophic_group %in% c("Non-herbivore", "Mixotroph") ~ "Right"))

#panel labels
panel_labs <- tibble::tibble(
  trophic_group = c("Herbivore", "Non-herbivore",
                    "Autotroph", "Mixotroph"),
  facet_row = c("Row1", "Row1", "Row2", "Row2"),
  facet_col = c("Left", "Right", "Left", "Right"),
  label = c("Herbivore", "Non-herbivore",
            "Autotroph", "Mixotroph"))

#Figure 6
ggplot(dat ,aes(x = log_chla, y = log_size, color = trophic_state)) +
  geom_point(alpha = 0.4, size = 2, show.legend = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
              se = TRUE, size = 0.9) +
  geom_text_npc(data = panel_labs,
                aes(npcx = 0.5, npcy = 1.4, label = label), inherit.aes = FALSE, size = 3.5) +
  facet_grid(facet_row ~ facet_col, scales = "free_y") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  labs(x = "Log10 chla", y = "Log10 mean size (µm)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank())
#ggsave("figs/gam_size_vs_chla.jpg", width = 5, height = 4)

#pull out stats for results (Table S4) NEEDS TO BE UPDATED IN RESULTS!!!!
extract_gam_stats <- function(model) {
  s <- summary(model)
  tibble(
    edf = round(s$s.table[1, "edf"],2),
    F = round(s$s.table[1, "F"],2),
    p_value = signif(s$s.table[1, "p-value"],3), 
    #r2 = round(s$r.sq,2),
    dev_expl = round(s$dev.expl * 100,2))}

results_table <- map_dfr(models, extract_gam_stats, .id = "trophic_group") |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  dplyr::select(-p_value) |>
  rename("p-value" = "p_adj")
#write.csv(results_table, "output/functional_group_GAMs.csv", row.names = F)

# Fit interactive GAMs per functional group
gams_interactive <- dat |>
  group_by(trophic_group) |>
  group_map(~ gam(
    log_size ~ trophic_state + s(log_chla, by = trophic_state, k = 6),
    data = .x,
    method = "REML"))

names(gams_interactive) <- unique(dat$trophic_group)

#extract smooth f-stats (approximate given gamma smooths)
extract_smooths <- function(model) {
  s_tab <- as.data.frame(summary(model)$s.table)
  s_tab <- s_tab |>
    tibble::rownames_to_column("term") |>
    dplyr::select(term, edf, F, p_value = `p-value`) |>
    dplyr::mutate(term = sub("s\\(log_chla\\):trophic_state", "", term)) |>
    dplyr::rename("trophic_state" = "term")
  return(s_tab)}

# Table S5 (NEEDS TO BE UPDATED IN RESULTS AND SI)
smooths_table <- imap_dfr(gams_interactive, ~ {
  extract_smooths(.x) |>
    dplyr::mutate(trophic_group = .y)}) |>
  dplyr::select(trophic_group, trophic_state, edf, F, p_value) |>
  dplyr::mutate(p_adj = round(p.adjust(p_value, method = "BH"),3),
                edf = round(edf,2),
                F = round(F, 2),
                p_value = round(p_value, 3)) |>
  dplyr::select(-p_value)
#write.csv(smooths_table, "output/functional_group_ts_GAMs.csv", row.names = F)

# spearman rank correlations to test for monotonic size relationships between functional group and trophic state
# Table S6 NEDDS TO BE UPDATED IN SI
cor_table <- dat |>
  dplyr::group_by(trophic_group, trophic_state) |>
  dplyr::summarise(cor = cor(log_size, log_chla, method = "spearman"),
                   p_value = cor.test(log_size, log_chla, method = "spearman")$p.value,
                   .groups = "drop") |>
  dplyr::mutate(p_value = round(p_value, 3),
                cor = round(cor,2))
#write.csv(cor_table, "output/size_correlation_table.csv", row.names = F)

#calculate percent positive slope
pos_summary <- fits_ordered |>
  group_by(trophic_state, trophic_group) |>
  summarize(n = n(), n_pos = sum(slope > 0, na.rm = TRUE),
            pct_pos = 100 * n_pos / n) |> 
  ungroup()

# visualize positive slopes too (Figure 4)
ggplot(fits_ordered, aes(x = trophic_group, y = slope, color = trophic_group)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_summary(fun = median, geom = "point", size = 4, shape = 18, color = "black") +
  facet_wrap(~trophic_state, ncol = 1,
             labeller = labeller(trophic_state = label_map)) +
  scale_color_manual(values = c(
    "Autotroph" = "#3E6E66", "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C", "Non-herbivore" = "#F68A4D")) +
  labs(x = "", y = "NBSS slope") +
  theme_minimal() +
  geom_text(data = pos_summary,
            aes(x = trophic_group, y = 2,
                label = paste0(round(pct_pos, 1), "%")),
            inherit.aes = FALSE,size = 3, vjust = 0) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(),
        text = element_text(size=9),
        legend.position = "none",
        legend.direction = "horizontal") 
#ggsave("figs/all_lake_slopes.jpg", width = 5, height = 5)

#values for results text
sum(pos_summary$n_pos[pos_summary$trophic_group=="Mixotroph"])/
  sum(pos_summary$n[pos_summary$trophic_group=="Mixotroph"]) * 100

sum(pos_summary$n_pos[pos_summary$trophic_group=="Autotroph"])/
  sum(pos_summary$n[pos_summary$trophic_group=="Autotroph"]) * 100

sum(pos_summary$n_pos[pos_summary$trophic_group=="Herbivore"])/
  sum(pos_summary$n[pos_summary$trophic_group=="Herbivore"]) * 100

sum(pos_summary$n_pos[pos_summary$trophic_group=="Non-herbivore"])/
  sum(pos_summary$n[pos_summary$trophic_group=="Non-herbivore"]) * 100

sum(pos_summary$n_pos[pos_summary$trophic_state=="Oligotrophic"])/
  sum(pos_summary$n[pos_summary$trophic_state=="Oligotrophic"]) * 100

sum(pos_summary$n_pos[pos_summary$trophic_state=="Mesotrophic"])/
  sum(pos_summary$n[pos_summary$trophic_state=="Mesotrophic"]) * 100

sum(pos_summary$n_pos[pos_summary$trophic_state=="Eutrophic"])/
  sum(pos_summary$n[pos_summary$trophic_state=="Eutrophic"]) * 100

# Get Canada map Figure S2
canada <- ne_countries(scale = "medium", country = "Canada", returnclass = "sf")

# SSA data with coordinates
ssa_df <- fits_summary |>
  left_join(pred |> dplyr::select(lake_id,latitude,longitude), by = "lake_id") |>
  filter(!is.na(latitude) & !is.na(longitude))

ggplot() +
  geom_sf(data = canada, fill = "gray90", color = "black") +
  geom_point(data = ssa_df, aes(x = longitude, y = latitude, 
                                color = trophic_state), size = 2) +
  labs(color = "") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal")
#ggsave("figs/ts_map.jpg", width = 5, height = 5)
