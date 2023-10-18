# Flood Grid Rescaler

[![DOI](https://zenodo.org/badge/547351392.svg)](https://zenodo.org/badge/latestdoi/547351392)

A QGIS Processing script for rescaling flood hazard rasters (water depth, water surface elevation (WSE), and DEM) that includes the following tools:
- Aggregation via Averaging: Two methods for coarsening flood hazard grids described in [Bryant et. al., (2022)]
- Downscaling: Three methods for interpolating a high resolution WSE grid using a high resolution DEM and simple hydraulic assumptions.

The algorithms are tested against QGIS 3.28.11 and require [WhiteboxTools for QGIS plugin](https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html).

## 1 Installation Instructions

### download the scripts
download the algorithms of interest from the [processing folder](floodrescaler/processing) to your local machine

### add to your QGIS profile
In the QGIS [Processing Toolbox](https://docs.qgis.org/3.22/en/docs/user_manual/processing/toolbox.html#the-toolbox), select the python icon drop down ![Scripts](/assets/mIconPythonFile.png) , and `Add Script to Toolbox...`. This should load new algorithms to the `Scripts/FloodRescaling` group on the Processing Toolbox as shown below:

![screen capture](/assets/processingToolbox_screengrab.png)

### install and configure WhiteboxTools
Ensure the [WhiteboxTools for QGIS plugin](https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html) is installed and configured. 
Note this is only required for the Downscaling/CostGrow algo. 

## 2 Use
Instructions are provided on the algorithm dialog

## 3 Example Data

Example data is provided in the [examples](/examples) folder

## 4 Attribution

If you use these tools for your work, please cite the following:

### Aggregation via Averaging

```
@article{bryant_bias_2023,
	title = {Bias in {Flood} {Hazard} {Grid} {Aggregation}},
	url = {https://onlinelibrary.wiley.com/doi/abs/10.1029/2023WR035100},
	doi = {10.1029/2023WR035100}, 
	journal = {Water Resources Research},
	author = {Bryant, Seth and Kreibich, Heidi and Merz, Bruno},
	year = {2023},
}
```


### Downscaling

```
@article{bryant_resolution_2023,
	title = {Resolution {Enhancement} of {Flood} {Inundation} {Grids} [in review]},
	doi = {10.5194/hess-2023-156},
	journal = {Hydrology and Earth System Sciences},
	author = {Bryant, Seth and Guy, Schumann and Heiko, Apel and Heidi, Kreibich and Bruno, Merz},
	year = {2023},
}


@inproceedings{lindsay_whitebox_2014,
	title = {The whitebox geospatial analysis tools project and open-access {GIS}},
	url = {https://jblindsay.github.io/ghrg/pubs/LindsayGISRUK2014.pdf},
	urldate = {2022-04-11},
	booktitle = {Proceedings of the {GIS} {Research} {UK} 22nd {Annual} {Conference}, {The} {University} of {Glasgow}},
	author = {Lindsay, JB},
	year = {2014},
	pages = {16--18},
}

```

## 5 Development

create a virtual environment from the supported QGIS version and the `./requirements.txt` file. 

add a ./definitions.py file similar to the below:

```
import pathlib, os
src_dir = os.path.dirname(os.path.abspath(__file__))

wrk_dir = os.path.expanduser('~')

#add your WhiteBoxTools executable file
wbt_exe = r'l:\06_SOFT\whitebox\v2.2.0\whitebox_tools.exe'
```
