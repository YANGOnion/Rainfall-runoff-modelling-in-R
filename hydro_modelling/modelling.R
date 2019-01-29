
############################### extracting areal forcing data for the watershed ####################

## to read the watershed polygon
library(raster)
poly_file='input/wts_polygon/66193.shp'
wts_poly=shapefile(poly_file)

area(wts_poly)/10^6
plot(wts_poly)

## to extract areal mean values of forcing data within the watershed polygon
library(velox)
forcing_file='input/temperature/1960.01.tif'
vx_forcing=velox(forcing_file)
areal_value=vx_forcing$extract(sp=wts_poly,fun=mean,small=T)[1,] # the return value is a matrix
head(areal_value)

## an example for plotting raster forcing data and watershed polygon
library(rasterVis)
ras_forcing=brick('input/temperature/1960.tif')
levelplot(crop(ras_forcing,wts_poly,snap='out')[[1]],margin=list(draw=F),colorkey=list(space='right'),
          par.settings=YlOrRdTheme(),main='Mean temperature 1960-01-01')+layer(sp.polygons(wts_poly))


## to output the data.table with date and forcing variables 
library(data.table)
dt=data.table(date=seq(as.Date('1960-01-01'),as.Date('1960-01-31'),1),tasmean=areal_value)
head(dt)

################################### model calibration and simulation ################################

source('gr_model.R')

library(data.table)
dt=fread('input/66193.csv')
dt[,date:=as.Date(date)]

## calibration
cab=grcab(dt=dt,start='1960-01-01',end='1979-12-31',snow=F)
cab$crit # calibration NSE

## simulation
sim=grsim(dt=dt,start='1980-01-01',end='2000-12-31',param=cab$param,snow=F)
sim$crit # validation NSE

## extracting simulated series
out=SumOutputGr4j(sim$output)
head(out)

##################################### hydrograph visualization ######################################

source('hydrograph.R')

## the data frame of simulated results
dt_sim=cbind(dt[date>=as.Date('1980-01-01')],out)
head(dt_sim)

hydrograph(dt=dt_sim,rain_var='rainfall',flow_var=c('runoff','Qsim'),start='2000-01-01',end='2000-12-31')



