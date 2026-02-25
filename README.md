# Biodiversity databases as underutilized resources for pathogen discovery: a quantitative synthesis of bat and rodent tissue collections in natural history museums
This repository contains code for data cleaning and analysis associated with <a href="https://doi.org/10.1101/2025.11.28.691153" target="_blank">this study</a>.

# Overview and files
This study investigates the utility of biodiversity databases for locating museum samples for pathogen detection work. Here, we filtered all bat and rodent records in GBIF to locate tissue samples for pathobiology work. The unfiltered rodent and bat GBIF data are stored <a href="https://zenodo.org/records/18663923" target="_blank">here</a>, along with the downloaded filtered tissue dataset. The (large) unfiltered files can be read into the "cleaning.R" script to generate the filtered dataset. To associate latitude, longitude, and region with countries for visualization and statistical models, the "LatLong.csv" and "georegion.csv" files here can be read into the scripts.

The filtered tissue dataset can be read into the "analysis.R" script to generate the figures and analyses included in the study. Note that for the temporal distribution figure, "dat_years.csv" (an intermediate file) is also required. This file is also generated in "cleaning.R" and contains the date of collection of *all* samples, not just those with viable tissue.

# Contributors
Maya Juman (mmj38@cam.ac.uk)
Colleen Cronin (colleencronin2.0@gmail.com)
Alex Richardson (richardsonalex@gmail.com)
