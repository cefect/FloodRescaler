# Flood Grid Rescaler

[![DOI](https://zenodo.org/badge/547351392.svg)](https://zenodo.org/badge/latestdoi/547351392)

QGIS Processing scripts for rescaling or changing the resolution of flood hazard rasters.

Tested against QGIS 3.34.5

## UPDATES
- **2024 05 17**: added *Water Grid from Complement* tool to easily convert between WSE and WSH grids
- **2024 01 19**: added *pixel* based isolated filter for the *Resolution Enhancer* tool


## 01 Description
Includes the following tools:
- **Resolution Enhancer (WSE)**: Three methods for interpolating a high resolution WSE grid using a high resolution DEM and simple hydraulic assumptions described in [Bryant et. al., (2024)](https://doi.org/10.5194/hess-28-575-2024). For more advanced applications and features, see the [FloodDownscaler2](https://github.com/cefect/FloodDownscaler2) project. Requires [WhiteboxTools for QGIS plugin](https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html).
- **Aggregation via Averaging**: Two methods for coarsening flood hazard grids described in [Bryant et. al., (2023)](https://doi.org/10.1029/2023WR035100)
- **Water Grid from Complement**: Two tools to generate a Water Surface Elevation (WSE) or Water Surface Height (WSH) grid from its complement and a DEM.

 

## 02 Installation Instructions

### download the scripts
Download the *floodrescaler.zip* archive from the [latest release](https://github.com/cefect/FloodRescaler/releases) then unzip to your local machine.

Alternatively, download the processing scripts of interest from the [processing folder](floodrescaler/processing) source code directly.  

### add to your QGIS profile
In the QGIS [Processing Toolbox](https://docs.qgis.org/3.22/en/docs/user_manual/processing/toolbox.html#the-toolbox), select the python icon drop down ![Scripts](/assets/mIconPythonFile.png) , and `Add Script to Toolbox...`. This should load new algorithms to the `Scripts/FloodRescaling` group on the Processing Toolbox as shown below:

![screen capture](/assets/processingToolbox_screengrab.png)

Alternatively, add the unzipped folder as a *Scripts folder* to your QGIS profile: Settings >  Processing > Scripts > Scripts folders

### install and configure WhiteboxTools
Ensure the [WhiteboxTools for QGIS plugin](https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html) is installed and configured. 
Note this is only required for the Downscaling/CostGrow algo.  Tested against v2.2.0.

## 03 Use
Instructions are provided on the algorithm dialog

## 04 Example Data

Example data is provided in the [examples](/examples) folder

## 05 Attribution

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


### Resolution Enhancement

```
@Article{hess-28-575-2024,
	AUTHOR = {Bryant, S. and Schumann, G. and Apel, H. and Kreibich, H. and Merz, B.},
	TITLE = {Technical Note: Resolution enhancement of flood inundation grids},
	JOURNAL = {Hydrology and Earth System Sciences},
	VOLUME = {28},
	YEAR = {2024},
	NUMBER = {3},
	PAGES = {575--588},
	URL = {https://hess.copernicus.org/articles/28/575/2024/},
	DOI = {10.5194/hess-28-575-2024}
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

## 06 Development

- create a virtual environment from the supported QGIS version and the `./requirements.txt` file. 
- add a ./definitions.py file similar to the below
- add the QGIS python plugins directory (`C:\Users\cef\AppData\Local\QGIS\QGIS3\profiles\dev\python\plugins`) to your PYTHONPATH. Needed by the downscaler to mimic the dependency on WBT.

### definitions.py

```
import pathlib, os
src_dir = os.path.dirname(os.path.abspath(__file__))

wrk_dir = os.path.expanduser('~')

#add your WhiteBoxTools executable file
wbt_exe = r'l:\06_SOFT\whitebox\v2.2.0\whitebox_tools.exe'
```
