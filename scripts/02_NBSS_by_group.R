# calculating size spectra slopes using functional feeding groups

pacman::p_load(tidyverse, dplyr, ggplot2, ggmap,
               rnaturalearth, rnaturalearthdata, 
               ARTool, ggpubr, mgcv, ggpp)

# calculate trophic state using chla (file modified from CP where NAs were replaced using miss forest or ecozonal means)
pred <- read.csv("data/chla.csv") |>
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
        axis.line = element_line(color = "black",linewidth = 0.2),
        axis.ticks = element_line(color = "black",linewidth = 0.2),
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
        axis.line = element_line(color = "black",linewidth = 0.2),
        axis.ticks = element_line(color = "black",linewidth = 0.2)) 

ggarrange(p_ts, p_fg, ncol = 1, nrow = 2,       
          heights = c(1.4, 1),  
          widths = c(2,1), common.legend = FALSE)
#ggsave("figs/plankton_sig_slopes.jpg", width = 4, height = 4)

# slope across functional groups (Figure 2)
spectra_summary <- spectra_lake_bins |>
  left_join(pred |> dplyr::select(trophic_state, key), by = "key") |>
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
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(expression(Log[10]~mean~size)) +
  scale_y_continuous(expression(Log[10]~NBSS)) +
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
ggplot(spectra_summary, aes(x = log_size, y = mean_log_nbss, color = trophic_group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +      
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
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

# overall quadratic 
quad_all <- lm(mean_log_nbss ~ log_size + I(log_size^2), data = spectra_summary)
summary(quad_all)        
broom::tidy(quad_all)     
broom::glance(quad_all)   

#Figure S4
lake_slope_summary <- fits_summary |>
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
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("Oligotrophic" = "#21908CFF", 
                                "Mesotrophic" = "#FDE725FF",
                                "Eutrophic" = "#440154FF")) +
  labs(y = "Mean NBSS slope ± SE", x = "")
#ggsave("figs/slope_ts_fg.jpg", width = 6, height = 5)

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

weighted_peaks <- dat |>
  group_by(trophic_state, trophic_group) |>
  group_modify(~{
    d <- density(
      x = .$mean_size,
      weights = .$total_biomass,
      na.rm = TRUE)
    tibble(mu = d$x[which.max(d$y)])}) |>
  ungroup() |>
  mutate(plankton = ifelse(trophic_group %in% c("Herbivore","Non-herbivore"),
                      "zooplankton","phytoplankton"))

#histogram of plankton size across trophic states (Figure 1)
ggplot(dat |> filter(), aes(x = mean_size, fill = trophic_group)) +
  geom_density(aes(weight = total_biomass), alpha = 0.7) +
  geom_vline(data = weighted_peaks, aes(xintercept = mu, color = trophic_group, 
                                        group = plankton),
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
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
#ggsave("figs/ts_size_density_plots_weighted.jpg", width = 6, height = 5)

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
  left_join(pred |> dplyr::select(key, chla_day), by = "key") |>
  group_by(trophic_group, trophic_state) |>
  summarise(median_slope = median(slope, na.rm = TRUE),
            median_logchla = median(log10(chla_day), na.rm = TRUE),
            .groups = "drop")

#slope vs chla (Figure 5)
fits_summary |>
  left_join(pred |> dplyr::select(key, chla_day), by = "key") |>
  mutate(log_chla = log10(chla_day),
         slope_dir = ifelse(slope >= 0, "positive", "negative")) |>
  ggplot(aes(x = log_chla, y = slope, color = trophic_state, fill = trophic_state)) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.6) +  
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
  labs(y = "Slope", x = expression(Log[10]~chlorophyll~italic(a)), color = "", fill = "") +
  theme_minimal()  +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.ticks = element_line(color = "black", linewidth = 0.2))
#ggsave("figs/slope_vs_chla.jpg", width = 5, height = 3)

#add log chla to dat
dat <- dat |> 
  left_join(pred |> dplyr::select(key, chla_day), by = "key") |> 
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
  labs(x = expression(Log[10]~chlorophyll~italic(a)), y = expression(Log[10]~mean~size)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4))
#ggsave("figs/gam_size_vs_chla.jpg", width = 5, height = 4)

#pull out stats for results (Table S4) 
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
  mutate(trophic_state = factor(trophic_state)) |>
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

# Table S5 
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
# Table S6 
cor_table <- dat |>
  dplyr::group_by(trophic_group, trophic_state) |>
  dplyr::summarise(n = sum(complete.cases(log_size, log_chla)),
                   cor = cor(log_size, log_chla, method = "spearman"),
                   p_value = cor.test(log_size, log_chla, method = "spearman")$p.value,
                   .groups = "drop") |>
  dplyr::mutate(p_value = round(p_value, 3),
                cor = round(cor,2))
#write.csv(cor_table, "output/size_correlation_table.csv", row.names = F)
#warning is fine because sample size is large enough

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
