# Rainfall-runoff-modelling-in-R

This repository provides functions for rainfall-runoff modelling in R, including watershed delineation, climatic data preparation, the setup of a lumped hydrological model, and hydrograph visualization. In addition, an example for the complete process of rainfall-runoff modelling in a medium-size watershed in China using all functions is documented as follows.

## 1. Watershed delineation

__Prerequisites__
* The R package required is `raster`. We also need the `arcpy` module in python which will be installed if you have the ArcGIS software in your computer.

* The `wts_extract` directory in this repository contains codes and data for watershed delineation. Please set the workspace to be this location.

__1.1 Creating a spatial point of the watershed outlet__

We first create a shapefile of a spatial point for the outlet of the watershed. Here we use the hydrological gauge coded 66193 located at (113.45, 27.66667) as the `create_point.R` shows.
```r
library(raster)
df=data.frame(id=66193,lon=113.45,lat=27.66667)
p=SpatialPoints(df, proj4string=CRS('+proj=longlat +datum=WGS84'))
p=SpatialPointsDataFrame(p, data=df)
shapefile(p, paste('input/point/',i,'.shp',sep=''), overwrite=TRUE)
```

__1.2 Extracting the watershed polygon__

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

## 2. Rainfall-runoff modelling

__Prerequisites__
* The R packages required are `raster`,`velox`,`data.table`, and `airGR`. The `rasterVis` package for visualizing raster forcing data is optional. For hydrograph plotting, `ggplot2`,`gridExtra`, and `grid`.
* The `hydro_modelling` directory in this repository contains codes and data for watershed delineation. Please set the workspace to be this location.

__2.1 Extracting areal forcing data for the watershed__

The codes can be found in `modelling.R` in the `hydro_modelling` directory in this repository.

We read the watershed polygon shapefile created by the watershed delineation process:
```r
library(raster)
poly_file='input/wts_polygon/66193.shp'
wts_poly=shapefile(poly_file)
```

We extract areal mean daily temperature of the watershed from 1960-01-01 to 1960-01-31 with 0.25 degree resolution from [GLDAS_CLSM025_D](https://disc.gsfc.nasa.gov/datasets/GLDAS_CLSM025_D_V2.0/summary?keywords=GLDAS) as an example. The temperature raster data have been processed to be a `.tif` file including 31 layers. Although the `raster` package provides functions for raster extraction, we recommand the [`velox`](http://philipphunziker.com/velox/) package as a high-performance package for raster extraction and manipulation written in C++ but implemented in R.
```r
library(velox)
forcing_file='input/temperature/1960.01.tif'
vx_forcing=velox(forcing_file)
areal_value=vx_forcing$extract(sp=wts_poly,fun=mean,small)[1,] # the return value is a matrix
head(areal_value)
# [1] 6.500975 5.566156 7.507562 6.176838 5.147455 3.996281
```

We can visulize the watershed boundary and surrounding rasters of the forcing data:
```r
library(rasterVis)
ras_forcing=brick('input/temperature/1960.tif')
levelplot(crop(ras_forcing,wts_poly,snap='out')[[1]],margin=list(draw=F),colorkey=list(space='right'),
          par.settings=YlOrRdTheme(),main='Mean temperature 1960-01-01')+layer(sp.polygons(wts_poly))
```
![alt text](https://github.com/YANGOnion/Rainfall-runoff-modelling-in-R/blob/master/hydro_modelling/output/areal_tem.png)

We can also output a data frame containing date and forcing time series.
```r
library(data.table)
dt=data.table(date=seq(as.Date('1960-01-01'),as.Date('1960-01-31'),1),tasmean=areal_value)
head(dt)
#          date  tasmean
# 1: 1960-01-01 6.500975
# 2: 1960-01-02 5.566156
# 3: 1960-01-03 7.507562
# 4: 1960-01-04 6.176838
# 5: 1960-01-05 5.147455
# 6: 1960-01-06 3.996281
```

__2.2 Model calibration and simulation__

The codes can be found in `modelling.R` in the `hydro_modelling` directory in this repository.

Provided that we already have the data frame of daily forcing data for the watershed, i.e., date, rainfall, observed runoff, temperature and potential evaporation. The data for watershed 66193 is saved in `hydro_modelling/input/66193.csv` in this repository. We read the data via `data.table` package. The data are from 1960-01-01 to 2000-12-31.

```r
library(data.table)
dt=fread('input/66193.csv')
dt[,date:=as.Date(date)]
```
We choose the [GR4J](https://webgr.irstea.fr/en/modeles/journalier-gr4j-2/) lumped hydrological model since it can be implemented from the well-developed package `airGR`. The [CemaNeige](https://webgr.irstea.fr/en/modeles/modele-de-neige/) snow module can be coupled in the GR4J. This repository combine all functions for building a model and wrap them into three functions in the `hydro_modelling/gr_model.R`:
* `grcab` function for calibration;
* `grsim` function for simulation;
* `SumOutputGr4j` function for extracting useful water flux variables related to snow and soil water.

We call these functions one by one with inout forcing and observed runoff data `dt` mentioned above, being a complete process of modelling. The data from 1960-01-01 to 1979-12-31 are for calibration and the remaining are for validation.

```r
source('gr_model.R')

## calibration
cab=grcab(dt=dt,start='1960-01-01',end='1979-12-31',snow=F)
cab$crit # calibration NSE
# [1] 0.8807795

## simulation
sim=grsim(dt=dt,start='1980-01-01',end='2000-12-31',param=cab$param,snow=F)
sim$crit # validation NSE
# [1] 0.9017213

## extracting simulated series
out=SumOutputGr4j(sim$output)
head(out)
#        Prod       Qsim
# 1: 176.2113 0.08535246
# 2: 177.7764 0.09517187
# 3: 177.6512 0.09136123
# 4: 179.1433 0.09034428
# 5: 180.1686 0.10051123
# 6: 180.6574 0.09378710

```
> Setting `snow=T` will activate the snow module.
> `Prod` is the prodction storage as well as the soil water storage and `Qsim` is the simulated runoff.

From the reaults we can konw that the NSE values for calibration and validation are 0.88 and 0.90 respectively, indicating good model performance.

__2.3 Hydrograph visualization__

We use `ggplot2` package to visualize the hydrograph with the function `hydrograph` in the `hydro_modelling/hydrograph.R`.

```r
source('hydrograph.R')

## the data frame of simulated results
dt_sim=cbind(dt[date>=as.Date('1980-01-01')],out)
head(dt_sim)
#          date    runoff  rainfall  tasmean       pet     Prod       Qsim
# 1: 1980-01-01 0.4855511 1.3600000 7.398309 0.5385461 176.2113 0.08535246
# 2: 1980-01-02 0.5683155 2.6400000 7.709050 0.8354931 177.7764 0.09517187
# 3: 1980-01-03 0.4689982 0.1766667 4.118336 0.3519521 177.6512 0.09136123
# 4: 1980-01-04 0.5517626 2.0066667 1.639363 0.2802955 179.1433 0.09034428
# 5: 1980-01-05 0.5876272 1.7783333 3.547596 0.5796179 180.1686 0.10051123
# 6: 1980-01-06 0.5186569 1.1733333 2.312646 0.5848638 180.6574 0.09378710

hydrograph(dt=dt_sim,rain_var='rainfall',flow_var=c('runoff','Qsim'),start='2000-01-01',end='2000-12-31')
```
![Alt text](https://github.com/YANGOnion/Rainfall-runoff-modelling-in-R/blob/master/hydro_modelling/output/hydrograph.png)
