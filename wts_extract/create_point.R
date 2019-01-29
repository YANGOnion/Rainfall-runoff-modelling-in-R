

library(raster)
df=data.frame(id=66193,lon=113.45,lat=27.66667)
p=SpatialPoints(df, proj4string=CRS('+proj=longlat +datum=WGS84'))
p=SpatialPointsDataFrame(p, data=df)
shapefile(p, paste('input/point/',i,'.shp',sep=''), overwrite=TRUE)

