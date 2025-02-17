# Qikiqtaruk Phenology Time-Series code repository

This code repository contains the analytical scripts and tabular data required to generate the figures, tables and statistics reported in *Assmann et al. 2020 - Drone data reveal heterogeneity in tundra greenness and phenology not captured by satellites*. The execution of the code also requires the spatial data files found in the accompanying data repository on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3904258.svg)](https://doi.org/10.5281/zenodo.3904258), as well as the retrival of Landsat data as detailed in the manuscript and the scripts contained in this repository.

# Content
1. [Authors](#Authors)
2. [Citation](#Citation)
3. [Acknowledgements](#Ackowledgements)
4. [Getting started](#Getting-started)
5. [License](#License)
---
# Authors
Jakob J. Assmann, Isla H. Myers-Smith, Jeffrey T. Kerby, Andrew M. Cunliffe, Gergana N. Daskalova

Contact: j.assmann@bios.au.dk

# Citation
Jakob J. Assmann, Isla H. Myers-Smith, Jeffrey T. Kerby, Andrew M. Cunliffe and Gergana N. Daskalova. 2020. Drone data reveal heterogeneity in tundra greenness and phenology not captured by satellites. Environmental Research Letters. 15. 125002. https://doi.org/10.1088/1748-9326/abbf7d

[\[back to top\]](#Qikiqtaruk-Phenology-Time-Series-Code-Repository)

# Acknowledgements 
(from the manuscript)

We would like to thank the Team Shrub field crews of the 2016 and 2017 field seasons for their hard work and effort invested in collecting the data presented in this research, this includes Will Palmer, Santeri Lehtonen, Callum Tyler, Sandra Angers-Blondin and Haydn Thomas. Furthermore, we would like to thank Tom Wade and Simon Gibson-Poole from the University of Edinburgh Airborne GeoSciences Facility, as well as Chris McLellan and Andrew Gray from the NERC Field Spectroscopy Facility for their support in our drone endeavours. We also want to express our gratitude to Ally Phillimore, Ed Midchard, Toke Høye and two anonymous reviewers for providing feedback on earlier versions of this manuscript. Lastly, JJA would like to thank IMS, Ally Phillimore and Richard Ennos for academic mentorship throughout his PhD.

We thank the Herschel Island—Qikiqtaruk Territorial Park Team and Yukon Government for providing logistical support for our field research on Qikiqtaruk including: Richard Gordon, Cameron Eckert and the park rangers Edward McLeod, Sam McLeod, Ricky Joe, Paden Lennie and Shane Goosen. We thank the research group of Hugues Lantuit at the Alfred Wegener Institute and the Aurora Research Institute for logistical support. Research permits include Yukon Researcher and Explorer permits (16-48S&E and 17-42S&E) and Yukon Parks Research permits (RE-Inu-02-16 and 17-RE-HI-02). All airborne activities were licensed under the Transport Canada special flight operations certificates ATS 16-17-00008441 RDIMS 11956834 (2016) and ATS 16-17-00072213 RDIMS 12929481 (2017).

Funding for this research was provided by NERC through the ShrubTundra standard grant (NE/M016323/1), a NERC E3 Doctoral Training Partnership PhD studentship for Jakob Assmann (NE/L002558/1), a research grant from the National Geographic Society (CP-061R-17), a Parrot Climate Innovation Grant, the Aarhus University Research Foundation, and the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement (754513) for Jeffrey Kerby, a NERC support case for use of the NERC Field Spectroscopy Facility (738.1115), equipment loans from the University of Edinburgh Airborne GeoSciences Facility and the NERC Geophysical Equipment Facility (GEF 1063 and 1069).

Finally, we would like to thank the Inuvialuit people for the opportunity to conduct research in the Inuvialuit Settlement Region.

[\[back to top\]](#Qikiqtaruk-Phenology-Time-Series-Code-Repository)

# Getting started

## Disclaimer

This code repository was designed to be as self-contained as possible. However, due to the limiations of git and GitHub in storing and version controlling large files we had to outsource the spatial data to a Zenodo repository. To complete the analysis these data need to be downloaded and all local references to the spatial data folder need to be updated. The relevant path variables are declared at the top of each script. 

Following suggestions from the reviewers, we added Landsat 8 Collection 1 Surface Reflectance images to the anlysis. These were too large for the repositories, but a list of all scenes included in the analysis can be found [here](data/auxillary/ls8_cloud_data.csv). This list should allow for retrieval of the scenes from the [USGS Eearth Explorer](https://earthexplorer.usgs.gov/). Furthermore, the VI values extracted from these Landsat 8 scenes are stored in tabular data on this code repository as indicated in the relevant scripts. 

## Software
The original analysis and raster manipulations were carried out using R v. 3.6.0 (R Core Team 2019) and the raster package v. 3.0-12 (Hijmans 2016). Visualisations were produced using the packages rasterVis v. 0.45 (Perpiñán and Hijmans 2018), ggplot v. 3.2.1 (Wickham 2016), cowplot v. 0.9.4 (Wilke 2019), and QGIS v. 3.10 (QGIS Development Team 2020).

## Repository folder structure
```
.
├── data                     Data produced or required by the analysis not in the Zenodo repository.
│   ├── auxillary               Flight log summaries and solar noon tables.    
│   ├── fig_*                   Intermediate and final data outputs associated with a figure. 
│   ├── modis                   MODIS time-series for the pixels covering the eight plots.
│   └── site_boundaries         Study site boundaries.
├── figures                  Figures generated by the analysis scripts and QGIS layout files.
├── log                      Log files for scripts with parallel processing.
└── scripts                  Scripts to prepare the data and conduct the analysis.
```

## Reproducing the analyses, figures and statistics
1. Run the three data preparation scripts (`scripts/prep_*`).
2. Run the relevant scripts (`scripts/fig_*`) to generate individual plots, figure panels or associated statistics.
3. Assemble the final figures from the individual plots / figure panels using the QGIS Layout Manager and the QGIS project files (`*.qgz`) in the `figures/fig_*/` folders.

Note: There are two ways of generating Figure 4: 1) aggregation then interpolation (used in the manuscript) or 2) interpolation and then aggregation (not used in the manuscript). A detailed explanation can be found in the scripts. The interpretation of the results is not affected by the method, but both are included in this repository for transparency. 

[\[back to top\]](#Qikiqtaruk-Phenology-Time-Series-Code-Repository)

# License 
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />The content of this repository is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

[\[back to top\]](#Qikiqtaruk-Phenology-Time-Series-Code-Repository)
