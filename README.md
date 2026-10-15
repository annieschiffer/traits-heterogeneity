# Testing the stress gradient hypothesis in invaded communities

## Overview

In this analysis, I test the stress gradient hypothesis along a cold stress gradient in eastern Idaho, 
where plant communities are invaded by cheatgrass (*Bromus tectorum*). My aims were to 1) quantify species interactions
between cheatgrass and neighboring shrubs and annuals in warm, low elevation sites vs. cold, high elevation sites; and 
2) determine how cheatgrass functional traits vary across competitive neighborhoods and environmental conditions.

## Methods

I established 6 sites at the Sheep Experiment Station, 3 at low elevation and 3 at high elevation. At each site,
I set up a 2x2 factorial experiment, in which I planted individual cheatgrass seeds in 2 patch types 
(under shrubs or in open spaces) and in 2 treatments (annual neighbors removed or annual neighbors present). This 
allowed me to test effects of shrubs and annuals on individual cheatgrass. I measured the demography (establishment
and aboveground biomass) and functional traits (SLA, SRL, height, establishment phenology, flowering phenology) of all
planted cheatgrass individuals.

## Analysis

I built separate Bayesian hierarchical models for each of the 2 demographic response variables and 5 functional trait
response variables. All models include elevation, patch type, treatment, and their interactions as predictor variables.
I also ran a PCA to test for differences in functional dispersion between favorable low elevation sites and stressful 
high elevation sites.

## Repository

All scripts and analyses can be run with `RUN_EVERYTHING.R`. This script calls all cleaning scripts, analyses, and
scripts for making figures that appear in the manuscript. All scripts are in the `code` folder, raw data are in the 
`raw_data` folder, cleaned data in the `clean_data` folder, and final results in the `outputs/final` folder.