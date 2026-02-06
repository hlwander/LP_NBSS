# LP_NBSS

This code reproduces figures from the analysis of plankton biomass-size relationships using data collected from \>550 Canadian lakes as part of the NSERC Lake Pulse Network project.

# Repository structure

This repository contains four folders: 1) `data` contains all the csv files required to generate figures in the analysis of biomass-size relationships across functional groups and lake trophic state, as well as open access plankton datasets for which a subset of phytoplankton taxon lengths were utilized (see the `supp_lengths` folder), 2) `figs` contains the final manuscript and supplemental figures for submission, 3) `output` contains all the supplemental tables (Table S2-S6) that summarize normalized biomass size spectra slope statistics, general additive model statistics, and correlations for plankton size among functional groups and lake trophic state, and 4) `scripts` contains all the code necessary to reproduce figures and analyses described above.

# Instructions to reproduce figures and analyses

1.  Run `01_supp_phyto_colony_lengths.R` to pull mean coloniy size data for 17 autotrophic phytoplankton taxa from 4 open access datasets. This script generates a csv file (`data/supp_lengths/collated_mean_phyto_lengths.csv`) that contains mean sizes for the 17 taxa to be used in the biomass-size analysis.
2.  Run `02_NBSS_by_group.R` to calculate normalized biomass for each lake and functional group, visualize normalized biomass size spectra relationships, and generate the final manuscript figures and tables.
