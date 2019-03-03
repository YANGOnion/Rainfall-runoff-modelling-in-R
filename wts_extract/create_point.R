

library(raster)
df=data.frame(lon=113.45,lat=27.66667)
p=SpatialPoints(df, proj4string=CRS('+proj=longlat +datum=WGS84'))
p=SpatialPointsDataFrame(p, data=df)
shapefile(p, 'input/point/66193.shp', overwrite=TRUE)

