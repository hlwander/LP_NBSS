# calculating size spectra slopes using functional feeding groups

pacman::p_load(tidyverse, dplyr, ggplot2, ggmap,
               rnaturalearth, rnaturalearthdata, 
               ARTool, ggpubr, mgcv, ggpp, emmeans, 
               multcomp, lmtest, sandwich)

# calculate trophic state using chla 
pred <- read.csv("data/env.csv") |>
  mutate(tsi_sd   = 60 - 14.41 * log(secchi), 
         tsi_chl  = 9.81 * log(chla) + 30.6,
         tsi_tp   = 14.42 * log(tp) + 4.15,
         tsi_mean = rowMeans(data.frame(tsi_sd, tsi_chl, tsi_tp), na.rm = TRUE),
         trophic_state = case_when(tsi_mean < 40 ~ "Oligotrophic",
                                   tsi_mean < 50 ~ "Mesotrophic",
                                   tsi_mean < 70 ~ "Eutrophic",
                                   is.nan(tsi_mean) ~ NA,
                                   TRUE ~ "Hypereutrophic")) |>
  filter(!is.na(trophic_state)) |> #removing lakes with NA values
  mutate(trophic_state = if_else(trophic_state == "Hypereutrophic", 
                                 "Eutrophic", trophic_state),
         trophic_state = factor(trophic_state, levels = c(
           "Oligotrophic", "Mesotrophic", "Eutrophic"))) 

# read in supplemental phyto length data
phyto_supp_lengths <- read.csv("data/supp_lengths/collated_mean_phyto_colony_lengths.csv") |>
  rename(supp_length_um = avg_length_um)

#read in plankton df
plankton_with_strategy <- read.csv("data/plankton_all.csv")

#define bin breaks
breaks <- seq(floor(min(plankton_with_strategy$log_length, na.rm = TRUE) / 0.1) * 0.1,
              ceiling(max(plankton_with_strategy$log_length, na.rm = TRUE) / 0.1) * 0.1, 
              by = 0.1)

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
  group_by(key, trophic_group, bin_mid, size_lower, size_upper, bin_width) |>
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
  group_by(key, trophic_group) |>  
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
  dplyr::select(key, trophic_group, n_bins, tidy, glance) |>
  unnest(tidy) |>
  pivot_wider(id_cols = c(key, trophic_group, n_bins),
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
  dplyr::select(key, trophic_group, glance) |>
  unnest_wider(glance)

#this drops 11 lakes that do not have trophic state classifications
fits_summary <- coeffs_wide |>
  left_join(glance_wide, by = c("key", "trophic_group")) |>
  left_join(pred |> dplyr::select(key, trophic_state), 
            by = "key") |>
  arrange(trophic_group, key) |>
  mutate(plankton = ifelse(trophic_group %in% c("Mixotroph","Autotroph"),
                           "Phytoplankton", "Zooplankton"),
         trophic_state = if_else(trophic_state == "Hypereutrophic", 
                                 "Eutrophic", trophic_state),
         trophic_state = factor(trophic_state, levels = c("Oligotrophic","Mesotrophic","Eutrophic"))) |>
  filter(!is.na(trophic_state)) 

#------------------------------------------------------------------------------#
#also a df for the whole community
community_bins <- plankton_with_strategy |>
  mutate(size_bin = cut(log_length, breaks = breaks, include.lowest = TRUE),
         bin_id = as.integer(size_bin),
         log_lower = breaks[bin_id],
         log_upper = breaks[bin_id + 1],
         size_lower = 10^log_lower,
         size_upper = 10^log_upper,
         bin_width = size_upper - size_lower,
         bin_mid = (log_lower + log_upper) / 2) |>
  filter(!is.na(biomass)) |>
  group_by(key, bin_mid, size_lower, size_upper, bin_width) |> #no trophic group here!
  summarise(total_biomass = sum(biomass, na.rm = TRUE), .groups = "drop") |>
  mutate(nbss = total_biomass / bin_width,
         log_nbss = ifelse(nbss > 0, log10(nbss), NA_real_),
         log_size = bin_mid)

#calculate NBSS as slope of log NBSS vs log size bin
community_fits <- spectra_lake_bins |>
  filter(!is.na(log_nbss), !is.na(log_size))  |>
  group_by(key) |>  
  summarise(n_bins = n(), model = list(
    if (n() >= 5) {
      tryCatch(lm(log_nbss ~ log_size, data = cur_data()),
               error = function(e) NULL)
    } else NULL), .groups = "drop") |>
  mutate(model_class = map(model, ~ {if (is.null(.x)) {"NULL"
  } else {paste(class(.x), collapse = "|")}
  }) |> unlist(),
  ok = map_lgl(model, ~ !is.null(.x) && inherits(.x, "lm")))

# Keep only successful fits
good_comm <- community_fits |> filter(ok) |> #keeping at least 5 size class bins
  mutate(tidy = map(model, broom::tidy), # 599/659 lakes 
         glance = map(model, broom::glance))

# Extract tidy coefficients and glance stats, then pivot coefficients to columns
coeffs_wide_comm <- good_comm |>
  dplyr::select(key, n_bins, tidy, glance) |>
  unnest(tidy) |>
  pivot_wider(id_cols = c(key, n_bins),
              names_from = term,
              values_from = c(estimate, std.error, p.value),
              names_glue = "{.value}__{term}")

#rename cols
names(coeffs_wide_comm) <- names(coeffs_wide_comm) %>%
  gsub("^estimate__\\(Intercept\\)$", "intercept", .) %>%
  gsub("^std.error__\\(Intercept\\)$", "intercept_se", .) %>%
  gsub("^p.value__\\(Intercept\\)$", "intercept_p", .) %>%
  gsub("^estimate__log_size$", "slope", .) %>%
  gsub("^std.error__log_size$", "slope_se", .) %>%
  gsub("^p.value__log_size$", "slope_p", .)

# Unpack glance stats into columns and join with coeffs
glance_wide_comm <- good_comm |>
  dplyr::select(key, glance) |>
  unnest_wider(glance)

#this drops 11 lakes that do not have trophic state classifications
fits_summary_comm <- coeffs_wide_comm |>
  left_join(glance_wide_comm, by = c("key")) |>
  left_join(pred |> dplyr::select(key, trophic_state), 
            by = "key") |>
  arrange(key) |>
  mutate(trophic_state = factor(trophic_state, levels = c(
    "Oligotrophic","Mesotrophic","Eutrophic"))) |>
  filter(!is.na(trophic_state)) 

#----------------------------------------------------------------------------#
# manuscript figs

#first try linear model with transformation and check if this violates assumptions
mod <- lm(slope ~ trophic_state * trophic_group, data = fits_summary)
summary(mod)
#indicates significant interaction 

shapiro.test(residuals(mod)) 
#p < 0.05 so residuals are not normal, but this test is sensitive so it's okay given ~600 lakes
car::leveneTest(slope ~ trophic_state * trophic_group, data = fits_summary)
#p < 0.05; so suggests unequal variance...

#diagnostic plots
plot(mod)

#ANOVA for significance
car::Anova(mod, type = 3)

#because levene's test indicates heterogeneity, checking hc3 coefficients here
hc3_results <- coeftest(mod, vcov = vcovHC(mod, type = "HC3"))

# Convert to dataframe and format (Table S4)
hc3_df <- as.data.frame(hc3_results[,]) |>
  tibble::rownames_to_column("term") |>
  rename(estimate  = Estimate,
         se        = `Std. Error`,
         t_value   = `t value`,
         p_value   = `Pr(>|t|)`) |>
  mutate(term = case_when(
      term == "(Intercept)"                                         ~ "Intercept (Oligotrophic × Autotroph)",
      term == "trophic_stateMesotrophic"                           ~ "Mesotrophic",
      term == "trophic_stateEutrophic"                             ~ "Eutrophic",
      term == "trophic_groupMixotroph"                             ~ "Mixotroph",
      term == "trophic_groupHerbivore"                             ~ "Herbivore",
      term == "trophic_groupNon-herbivore"                         ~ "Non-herbivore",
      term == "trophic_stateMesotrophic:trophic_groupMixotroph"    ~ "Mesotrophic × Mixotroph",
      term == "trophic_stateEutrophic:trophic_groupMixotroph"      ~ "Eutrophic × Mixotroph",
      term == "trophic_stateMesotrophic:trophic_groupHerbivore"    ~ "Mesotrophic × Herbivore",
      term == "trophic_stateEutrophic:trophic_groupHerbivore"      ~ "Eutrophic × Herbivore",
      term == "trophic_stateMesotrophic:trophic_groupNon-herbivore"~ "Mesotrophic × Non-herbivore",
      term == "trophic_stateEutrophic:trophic_groupNon-herbivore"  ~ "Eutrophic × Non-herbivore",
      TRUE ~ term),
    estimate = round(estimate, 3),
    se       = round(se, 3),
    t_value  = round(t_value, 3),
    p_value  = case_when(p_value < 0.001 ~ "<0.001", TRUE
                         ~ as.character(round(p_value, 3))))
#write.csv(hc3_df, "output/hc3_coefficients.csv", row.names = FALSE)

#create table
emm <- emmeans(mod, pairwise ~ trophic_state | trophic_group)

emm_df <- as.data.frame(emm$emmeans)

#table 1
contr_df <- as.data.frame(emm$contrasts) |> dplyr::select(-c(df, t.ratio)) |>
  mutate(estimate = round(estimate, 3), SE = round(SE, 3),
    p.value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))
#write.csv(contr_df, "output/contrasts_table.csv", row.names = FALSE)

#------------------------------------------------------------------------------#
#same as above but for the whole plankton community
mod_full <- lm(slope ~ trophic_state, data = fits_summary_comm)
summary(mod_full)

shapiro.test(residuals(mod_full)) 
#p < 0.05 so residuals are not normal, but this test is sensitive so it's okay given ~600 lakes
car::leveneTest(slope ~ trophic_state, data = fits_summary)
#p < 0.12; so suggests equal variance!!

#diagnostic plots
plot(mod_full)

#ANOVA for significance
car::Anova(mod_full, type = 3)
#no sig differences across trophic states

#------------------------------------------------------------------------------#
# slope across functional groups (Figure 2)
spectra_summary <- spectra_lake_bins |>
  left_join(pred |> dplyr::select(trophic_state, key), by = "key") |>
  group_by(trophic_group, trophic_state, log_size) |> 
  summarise(mean_nbss = mean(nbss, na.rm = TRUE),
            sd_nbss   = sd(nbss, na.rm = TRUE),
            n = n(), .groups = "drop") |>
  mutate(se_nbss = sd_nbss / sqrt(n),
         log_mean_nbss = ifelse(mean_nbss > 0, log10(mean_nbss), NA_real_), #transform after averaging
         trophic_state = if_else(trophic_state == "Hypereutrophic", 
                                 "Eutrophic", trophic_state),
         se_log_nbss = ifelse(mean_nbss > 0,
                              se_nbss / (mean_nbss * log(10)),
                              NA_real_)) |>
  mutate(trophic_state = factor(trophic_state, levels = c(
           "Oligotrophic","Mesotrophic","Eutrophic"))) |>
  filter(!is.nan(log_mean_nbss),
         !is.na(trophic_group),
         !is.na(trophic_state))

community_summary <- community_bins |>
  left_join(pred |> dplyr::select(trophic_state, key), by = "key") |>
  filter(!is.na(trophic_state))

ggplot(spectra_summary, aes(x = log_size, y = log_mean_nbss, 
                            color = trophic_group)) +
  geom_point(show.legend = T) +
  geom_smooth(method = "lm", se = FALSE,  show.legend = FALSE) +
  geom_smooth(data = community_summary,
              aes(x = log_size, y = log_nbss),
              method = "lm",
              color = "black",
              linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(expression(Log[10]~~size~(bin~midpoint))) +
  scale_y_continuous(expression(Log[10]~normalized~biomass)) +
  scale_color_manual(values = c(
    "Autotroph" = "#3E6E66",
    "Mixotroph"  = "#739A88",
    "Herbivore"    = "#DE482C",
    "Non-herbivore" = "#F68A4D")) +
  theme_minimal() + 
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
#ggsave("figs/slopes_by_group.jpg", width = 5, height = 4)

#slopes faceted by trophic state (Figure S3)
ggplot(spectra_summary, aes(x = log_size, y = log_mean_nbss, color = trophic_group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +      
  geom_smooth(data = community_summary, aes(x = log_size, y = log_nbss), 
              method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  facet_wrap(~ trophic_state, nrow=1) +                                         
  scale_x_continuous(expression(Log[10]~~size~(bin~midpoint))) +
  scale_y_continuous(expression(Log[10]~normalized~biomass)) +
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
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
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

#Figure 3
lake_slope_summary <- fits_summary |>
  group_by(trophic_state, trophic_group) |>
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),        
    sd_slope   = sd(slope, na.rm = TRUE),          
    n_lakes    = n(),                              
    se_slope   = sd_slope / sqrt(n_lakes),         
    .groups = "drop") |>
  mutate(letters = ifelse(trophic_group %in% "Autotroph" & 
                            trophic_state == "Eutrophic", "b",
                          ifelse(trophic_group %in% "Autotroph" & 
                                   !trophic_state == "Eutrophic", "a", "")))

ggplot(lake_slope_summary, aes(x = trophic_group, 
                               y = mean_slope, color = trophic_state)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope),
                position = position_dodge(width = 0.6), width = 0.2) +
  geom_text(aes(label = letters, group = trophic_state), 
            position = position_dodge(width = 0.6),color = "black",
            vjust = -2.2, size = 4, show.legend = FALSE) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  labs(y = "Mean NBSS slope", x = "")
#ggsave("figs/slope_ts_fg.jpg", width = 6, height = 5)

#Figure S4 all plankton
lake_slope_summary_comm <- fits_summary_comm |>
  group_by(trophic_state) |>
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),        
    sd_slope   = sd(slope, na.rm = TRUE),          
    n_lakes    = n(),                              
    se_slope   = sd_slope / sqrt(n_lakes),         
    .groups = "drop") 

ggplot(lake_slope_summary_comm, aes(x = trophic_state, 
                               y = mean_slope, color = trophic_state)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope),
                position = position_dodge(width = 0.6), width = 0.2) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  labs(y = "Mean NBSS slope", x = "")
#ggsave("figs/slope_ts_all_plankton.jpg", width = 6, height = 5)

#raw plankton biomass by group
biomass_abs <- plankton_with_strategy |>
  left_join(pred |> dplyr::select(key, trophic_state), by = "key") |>
  group_by(key, trophic_group, trophic_state) |>
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
  labs(x = "Mean size (µm)", y = "Frequency", fill = "") +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
#ggsave("figs/ts_size_freq_plots.jpg", width = 6, height = 5)

#values for results text
median(dat$mean_size[dat$trophic_group=="Autotroph"])
median(dat$mean_size[dat$trophic_group=="Mixotroph"])
median(dat$mean_size[dat$trophic_group=="Herbivore"])
median(dat$mean_size[dat$trophic_group=="Non-herbivore"])

#calculate IQRs
autotroph_medians <- dat |>
  filter(trophic_group == "Autotroph") |>
  group_by(key) |>
  summarize(med_size = median(mean_size, na.rm = TRUE), .groups = "drop")

IQR(autotroph_medians$med_size, na.rm = TRUE)

mixotroph_medians <- dat |>
  filter(trophic_group == "Mixotroph") |>
  group_by(key) |>
  summarize(med_size = median(mean_size, na.rm = TRUE), .groups = "drop")

IQR(mixotroph_medians$med_size, na.rm = TRUE)

herbivore_medians <- dat |>
  filter(trophic_group == "Herbivore") |>
  group_by(key) |>
  summarize(med_size = median(mean_size, na.rm = TRUE), .groups = "drop")

IQR(herbivore_medians$med_size, na.rm = TRUE)

nonherb_medians <- dat |>
  filter(trophic_group == "Non-herbivore") |>
  group_by(key) |>
  summarize(med_size = median(mean_size, na.rm = TRUE), .groups = "drop")

IQR(nonherb_medians$med_size, na.rm = TRUE)


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
median_points <- fits_summary |>
  left_join(pred |> dplyr::select(key, chla, secchi, tp), by = "key") |>
  group_by(trophic_group, trophic_state) |>
  summarise(median_slope = median(slope, na.rm = TRUE),
            median_logchla = median(log10(chla), na.rm = TRUE),
            median_logsecchi = median(log10(secchi), na.rm = TRUE),
            median_logtp = median(log10(tp), na.rm = TRUE),
            .groups = "drop")

#slope vs chla (Figure 5)
fits_summary |>
  left_join(pred |> dplyr::select(key, chla), by = "key") |>
  mutate(log_chla = log10(chla),
         slope_dir = ifelse(slope >= 0, "positive", "negative")) |>
  ggplot(aes(x = log_chla, y = slope, color = trophic_state, fill = trophic_state)) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.6) +  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, 
              color = "black", show.legend = FALSE) +
  geom_point(data = median_points,
             aes(y = median_slope, x = median_logchla, fill = trophic_state),
             shape = 21, size = 2, stroke = 0.8, color = "red",inherit.aes = FALSE) +
  facet_wrap(~trophic_group, ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  scale_fill_manual(values = c("Oligotrophic" = "#21908CFF", 
                               "Mesotrophic" = "#FDE725FF",
                               "Eutrophic" = "#440154FF")) +
  labs(y = "NBSS slope", x = expression(Log[10]~chlorophyll~italic(a)), color = "", fill = "") +
  theme_minimal()  +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.ticks = element_line(color = "black", linewidth = 0.2))
#ggsave("figs/slope_vs_chla.jpg", width = 5, height = 4)

#figure S6 (tp)
fits_summary |>
  left_join(pred |> dplyr::select(key, tp), by = "key") |>
  mutate(log_chla = log10(tp),
         slope_dir = ifelse(slope >= 0, "positive", "negative")) |>
  ggplot(aes(x = log_chla, y = slope, color = trophic_state, fill = trophic_state)) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.6) +  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, 
              color = "black", show.legend = FALSE) +
  geom_point(data = median_points,
             aes(y = median_slope, x = median_logtp, fill = trophic_state),
             shape = 21, size = 2, stroke = 0.8, color = "red",inherit.aes = FALSE) +
  facet_wrap(~trophic_group, ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  scale_fill_manual(values = c("Oligotrophic" = "#21908CFF", 
                               "Mesotrophic" = "#FDE725FF",
                               "Eutrophic" = "#440154FF")) +
  labs(y = "NBSS slope", x = expression(Log[10]~TP), color = "", fill = "") +
  theme_minimal()  +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.ticks = element_line(color = "black", linewidth = 0.2))
#ggsave("figs/slope_vs_tp.jpg", width = 5, height = 4)

#figure S7 (Secchi)
fits_summary |>
  left_join(pred |> dplyr::select(key, secchi), by = "key") |>
  mutate(log_chla = log10(secchi),
         slope_dir = ifelse(slope >= 0, "positive", "negative")) |>
  ggplot(aes(x = log_chla, y = slope, color = trophic_state, fill = trophic_state)) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.6) +  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, 
              color = "black", show.legend = FALSE) +
  geom_point(data = median_points,
             aes(y = median_slope, x = median_logsecchi, fill = trophic_state),
             shape = 21, size = 2, stroke = 0.8, color = "red",inherit.aes = FALSE) +
  facet_wrap(~trophic_group, ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  scale_fill_manual(values = c("Oligotrophic" = "#21908CFF", 
                               "Mesotrophic" = "#FDE725FF",
                               "Eutrophic" = "#440154FF")) +
  labs(y = "NBSS slope", x = expression(Log[10]~Secchi), color = "", fill = "") +
  theme_minimal()  +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.ticks = element_line(color = "black", linewidth = 0.2))
#ggsave("figs/slope_vs_secchi.jpg", width = 5, height = 4)

#calculate percent positive slope
pos_summary <- fits_summary |>
  group_by(trophic_state, trophic_group) |>
  summarize(n = n(), n_pos = sum(slope > 0, na.rm = TRUE),
            pct_pos = 100 * n_pos / n) |> 
  ungroup()

#count hopw many lakes in each trophic state
lake_counts <- fits_summary  |>
  distinct(key, trophic_state) |>  
  count(trophic_state) |>              
  mutate(label = paste0(trophic_state, " (n = ", n, ")"))

label_map <- setNames(lake_counts$label, lake_counts$trophic_state)

# visualize positive slopes too (Figure 4)
ggplot(fits_summary, aes(x = trophic_group, y = slope, color = trophic_group)) +
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
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        text = element_text(size=9),
        legend.position = "none",
        legend.direction = "horizontal") 
#ggsave("figs/all_lake_slopes.jpg", width = 5, height = 5)

#positive slopes for all plankton (Figure S5)
pos_summary_comm <- fits_summary_comm |>
  group_by(trophic_state) |>
  summarize(n = n(), n_pos = sum(slope > 0, na.rm = TRUE),
            pct_pos = 100 * n_pos / n) |> 
  ungroup()

ggplot(fits_summary_comm, aes(x = trophic_state, y = slope, color = trophic_state)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_summary(fun = median, geom = "point", size = 4, shape = 18, color = "black") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  labs(x = "", y = "NBSS slope") +
  theme_minimal() +
  geom_text(data = pos_summary_comm,
            aes(x = trophic_state, y = 2,
                label = paste0(round(pct_pos, 1), "%")),
            inherit.aes = FALSE,size = 3, vjust = 0) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        text = element_text(size=9),
        legend.position = "none",
        legend.direction = "horizontal") 
#ggsave("figs/all_lake_slopes_all_plankton.jpg", width = 5, height = 5)

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

# lake coordinates
lat_long_df <- read.csv("data/map_data.csv") |>
  filter(!is.na(latitude) & !is.na(longitude))

ggplot() +
  geom_sf(data = canada, fill = "gray90", color = "black") +
  geom_point(data = lat_long_df, aes(x = longitude, y = latitude, 
                                color = trophic_state), size = 2) +
  labs(color = "") +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal")
#ggsave("figs/ts_map.jpg", width = 5, height = 5)
