# Treatment Pathway Exploration

This repository contains the R code used to explore patient treatment pathways using time series clustering.

It contains the following folders:

-   [code](./code)

Includes functions and code used in the documentation notebooks.

-   [data](./data)

Here is where you need to store any relevant files for the analysis. Original files for project are not available.

The original files which are absent were called:

`pat_level_info.csv` for patient level information (patient characteristics)

`GLP-1 RA first-use as anti-obesity medication (no diabetes)_Event.csv` for the event log

`Prescription Duration Rules Table_data_20240624.csv` for prescription rules for missing prescription duration

-   [documentation](./documentation)

A series of R notebooks which guide you through all of the analysis of the project from tranforming the event log into patient trajectories to stratifying patient characteristics by cluster.

-   [models](./models)

This folder stores any models you decide to save down for later use.

-   [results](./results)

This folder stores any results produced from the analyses including image files.

# R packages used

Many different R packages will need to be installed to run the R notebooks and code, including (but not exhaustive):

-   `TraMineR` - for creating sequence objects and most of the visualisations

-   `WeightedCluster` - for average silhouette width measure

-   `Entropy` - for measuring Shannon's Entropy

-   `R6` - for using the R6 class for creating models
