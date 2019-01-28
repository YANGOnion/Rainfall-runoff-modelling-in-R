# Rainfall-runoff-modelling-in-R

This repository provides functions for rainfall-runoff modelling in R, including watershed delineation, climatic data preparation, the setup of a lumped hydrological model, and hydrograph visualization. In addition, an example for the complete process of rainfall-runoff modelling in a medium-size watershed in China using all functions is documented as follows.

## Watershed delineation

__Prerequisites__
* The R package required is `raster`. We also need the `arcpy` module in python which will be installed if you have the ArcGIS software in your computer.

* The `wts_extract` directory in this repository contains codes and data for watershed delineation. Please set the workspace to be this location.

__1. Creating a spatial point of the watershed outlet__

We first create a shapefile of a spatial point for the outlet of the watershed. Here we use the hydrological gauge coded 66193 located at (113.45, 27.66667) as the `create_point.R` shows.
```r
library(raster)
df=data.frame(id=66193,lon=113.45,lat=27.66667)
p=SpatialPoints(df, proj4string=CRS('+proj=longlat +datum=WGS84'))
p=SpatialPointsDataFrame(p, data=df)
shapefile(p, paste('input/point/',i,'.shp',sep=''), overwrite=TRUE)
```

__2. Extracting the watershed polygon__

ArcGIS is an excellent tool for watershed delineation. Note that the following codes are written in python in order to use the `arcpy` module in ArcGIS. As written in the `wts_extract.py`, the function of extracting the polygon of a watershed from flow direction map and flow accumulation map is:
```python
import arcpy

def wts_extrat(dir_raster,acc_raster,point,snap_output,Watersh_as_d1,RasterT_Watersh1):
	"""Extract the polygon of a watershed from flow direction map and flow accumulation map.
    Args:
        dir_raster: the input file of a direction map
        acc_raster: the input file of an accumulation map
        point: the input file of a point object
		snap_output: the output file of the snap pour point
		Watersh_as_d1: the output file of the watershed raster
		RasterT_Watersh1: the output file of the watershed polygon
    Returns:
        None
    """
	# Process: 捕捉倾泻点
	arcpy.gp.SnapPourPoint_sa(point, acc_raster, snap_output, ".015", "id")
	# Process: 分水岭
	arcpy.gp.Watershed_sa(dir_raster, snap_output, Watersh_as_d1, "VALUE")
	# Process: 栅格转面
	arcpy.RasterToPolygon_conversion(Watersh_as_d1, RasterT_Watersh1, "SIMPLIFY", "VALUE")
	return None
```

We need a [flow direction map](http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/flow-direction.htm) and a [flow accumulation map](http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/flow-accumulation.htm). These maps can be obtained from [HydroSHEDS](https://hydrosheds.cr.usgs.gov/dataavail.php). Here we use the data in Asia with 15 arcmin resolution named `as_dir_15s` and `as_acc_15s` for flow direction map and flow accumulation map respectively. We specify the parameters of the `wts_extrat` function, i.e. the input & output file paths, as follows:

```python
## input gis objects
dir_raster = "input\\as_dir_15s\\as_dir_15s"
acc_raster = "input\\as_acc_15s\\as_acc_15s"
i=66193 # id (or file name) of the watershed
point = "input\\point\\"+str(i)+".shp"

## output gis objects
snap_output = "output\\snap\\"+str(i)
Watersh_as_d1 = "output\\wts_raster\\"+str(i)
RasterT_Watersh1 = "output\\wts_polygon\\"+str(i)+".shp"
```

And then we call the function:
```python
wts_extrat(dir_raster,acc_raster,point,snap_output,Watersh_as_d1,RasterT_Watersh1)
```

The shapefile of the watershed polygon is saved in the `output/wts_polygon` directory.

## Rainfall-runoff modelling

__Prerequisites__
* 




