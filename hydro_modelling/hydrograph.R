
#' to plot the hydrograph of a watershed
#' @param dt a data.table of date and hydrological time series, e.g., rainfall and runoff
#' @param rain_var the name of the rainfall column in the dt
#' @param flow_var the names of the observed and simulated runoff columns in the dt
#' @param start the starting date to be ploted
#' @param end the ending date to be ploted
#' @return to plot the hydrograph and return NULL
hydrograph=function(dt,rain_var,flow_var,start,end){
  require(data.table)
  require(gridExtra)
  require(ggplot2)
  require(grid)
  dt_daily=copy(dt)
  dt_daily=dt_daily[,date:=as.Date(date)][date>=as.Date(start)][date<=as.Date(end)]
  dt_daily=melt(dt_daily,measure.vars=flow_var)
  dt_daily[,rainfall:=dt_daily[,rain_var,with=F][[1]]]
  ## rain plot
  g.top <- ggplot(dt_daily, aes(x = date, y = rainfall)) +
    geom_bar(stat = 'identity', fill = "steelblue1") +
    scale_y_continuous(trans = "reverse") +
    theme_bw(base_size = 15) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          plot.margin = unit(c(1,5,-10,6),units="points"),
          axis.title.y = element_text(vjust =0.3))+
    labs(y = "Rainfall (mm)")
  ## flow plot
  g.bottom <- ggplot(dt_daily, aes(x = date, y = value, color=variable)) +
    geom_line() +
    theme_bw(base_size = 15) +
    theme(plot.margin = unit(c(0,5,1,1),units="points"),legend.position = c(0.9, 0.9),
          legend.title = element_blank()) +
    labs(x = "Date", y = "Runoff (mm)")
  ## layouting plots
  g.top <- ggplotGrob(g.top)
  g.bottom <- ggplotGrob(g.bottom)
  maxWidths <- unit.pmax(g.top$widths, g.bottom$widths)
  g.top$widths=maxWidths
  g.bottom$widths=maxWidths
  grid.arrange(g.top,g.bottom, heights = c(0.25, 0.75))
  return(NULL)
}