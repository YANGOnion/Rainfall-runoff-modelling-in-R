
#' Calibration of GR4J model for one watershed
#' @param dt a data.table of input: 'date' for date, 'rainfall' for precipitation, 'pet' for potential evaporation, 'runoff' for observed runoff
#' @param start the date of the starting series formatting as 'yyyy-mm-dd'
#' @param end the date of the ending series formatting as 'yyyy-mm-dd'
#' @param plot whether to plot results
#' @return a list of: output series, parameters, and NSE
gr4jcab=function(dt,start,end,plot=F){
  require(airGR)
  require(data.table)
  dt_daily=copy(dt)
  dt_daily[,date:=as.Date(date)]
  dt_daily[,date:=as.POSIXct(date)]
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = dt_daily$date,
                                   Precip = dt_daily$rainfall, PotEvap = dt_daily$pet)
  Ind_Run <- seq(which(format(dt_daily$date, format = "%Y-%m-%d")==start),
                 which(format(dt_daily$date, format = "%Y-%m-%d")==end))
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                 InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                                 RunOptions = RunOptions, Qobs = dt_daily$runoff[Ind_Run])
  CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)
  OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                     InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                     FUN_MOD = RunModel_GR4J, FUN_CRIT = ErrorCrit_NSE)
  OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = OutputsCalib$ParamFinalR)
  if(plot) plot(OutputsModel, Qobs = dt_daily$runoff[Ind_Run])
  return(list(output=OutputsModel,param=OutputsCalib$ParamFinalR,crit=OutputsCalib$CritFinal))
}

#' Simulation of GR4J model for one MOPEX watershed
#' @param dt a data.table of input: 'date' for date, 'rainfall' for precipitation, 'pet' for potential evaporation
#' @param start the date of the starting series formatting as 'yyyy-mm-dd'
#' @param end the date of the ending series formatting as 'yyyy-mm-dd'
#' @param Param the parameters of GR4J model
#' @param plot whether to plot results
#' @return a list of: output series and NSE
gr4jsim=function(dt,start,end,param,plot=F){
  require(data.table)
  require(airGR)
  dt_daily=copy(dt)
  dt_daily[,date:=as.Date(date)]
  dt_daily[,date:=as.POSIXct(date)]
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = dt_daily$date,
                                   Precip = dt_daily$rainfall, PotEvap = dt_daily$pet)
  Ind_Run <- seq(which(format(dt_daily$date, format = "%Y-%m-%d")==start),
                 which(format(dt_daily$date, format = "%Y-%m-%d")==end))
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                 InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param)
  if(plot==T) plot(OutputsModel, Qobs = dt_daily$runoff[Ind_Run])
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
                                 RunOptions = RunOptions, Qobs = dt_daily$runoff[Ind_Run])
  OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)
  return(list(output=OutputsModel,crit=OutputsCrit$CritValue))
}



#' Calibration of CemaNeige-GR4J model for one watershed
#' @param dt a data.table of input: 'date' for date, 'rainfall' for precipitation, 'pet' for potential evaporation, 'runoff' for observed runoff,
#           'tasmean' for temperature. 'tasmax' for maximum and 'tasmin' for minimum daily temperature can be used 
#           if 'tasmean' misses.
#' @param start the date of the starting series formatting as 'yyyy-mm-dd'
#' @param end the date of the ending series formatting as 'yyyy-mm-dd'
#' @param elev the mean elevation of the watershed
#' @param elev_band  the vector of min, q01 to q99 and max of catchment elevation distribution 
#' @param plot whether to plot results
#' @return a list of: output series, parameters, and NSE
CemaNeigeGR4jcab=function(dt,start,end,elev,elev_band=NULL,plot=F){
  require(airGR)
  require(data.table)
  dt_daily=copy(dt)
  dt_daily[,date:=as.Date(date)]
  dt_daily[,date:=as.POSIXct(date)]
  if(!'tasmean'%in%names(dt_daily)) dt_daily[,tasmean:=(tasmax+tasmin)/2]
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, DatesR = dt_daily$date,
                                   Precip = dt_daily$rainfall, PotEvap = dt_daily$pet, TempMean = dt_daily$tasmean, 
                                   ZInputs = elev, HypsoData = elev_band, NLayers = 5)
  Ind_Run <- seq(which(format(dt_daily$date, format = "%Y-%m-%d")==start),
                 which(format(dt_daily$date, format = "%Y-%m-%d")==end))
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J,
                                 InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                                 RunOptions = RunOptions, Qobs = dt_daily$runoff[Ind_Run])
  CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_CemaNeigeGR4J, FUN_CALIB = Calibration_Michel)
  OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                     InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                     FUN_MOD = RunModel_CemaNeigeGR4J, FUN_CRIT = ErrorCrit_NSE)
  OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = OutputsCalib$ParamFinalR)
  if(plot) plot(OutputsModel, Qobs = dt_daily$runoff[Ind_Run])
  return(list(output=OutputsModel,param=OutputsCalib$ParamFinalR,crit=OutputsCalib$CritFinal))
}


#' Simulation of CemaNeige-GR4J model for one MOPEX watershed
#' @param dt a data.table of input: 'date' for date, 'rainfall' for precipitation, 'pet' for potential evaporation, 'runoff' for observed runoff,
#           'tasmean' for temperature. 'tasmax' for maximum and 'tasmin' for minimum daily temperature can be used 
#           if 'tasmean' misses.
#' @param start the date of the starting series formatting as 'yyyy-mm-dd'
#' @param end the date of the ending series formatting as 'yyyy-mm-dd'
#' @param Param the parameters of GR4J model
#' @param elev the mean elevation of the watershed
#' @param elev_band  the vector of min, q01 to q99 and max of catchment elevation distribution 
#' @param plot whether to plot results
#' @param plot whether plot results
#' @return list of: output series, NSE
CemaNeigeGR4jsim=function(dt,start,end,param,elev,elev_band=NULL,plot=F){
  require(data.table)
  require(airGR)
  dt_daily=copy(dt)
  dt_daily[,date:=as.Date(date)]
  dt_daily[,date:=as.POSIXct(date)]
  if(!'tasmean'%in%names(dt_daily)) dt_daily[,tasmean:=(tasmax+tasmin)/2]
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, DatesR = dt_daily$date,
                                   Precip = dt_daily$rainfall, PotEvap = dt_daily$pet, TempMean = dt_daily$tasmean, 
                                   ZInputs = elev, HypsoData = elev_band, NLayers = 5)
  Ind_Run <- seq(which(format(dt_daily$date, format = "%Y-%m-%d")==start),
                 which(format(dt_daily$date, format = "%Y-%m-%d")==end))
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J,
                                 InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param)
  if(plot==T) plot(OutputsModel, Qobs = dt_daily$runoff[Ind_Run])
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
                                 RunOptions = RunOptions, Qobs = dt_daily$runoff[Ind_Run])
  OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)
  return(list(output=OutputsModel,crit=OutputsCrit$CritValue))
}


#' Combine calibration of GR4j & CemaNeige-GR4j
grcab=function(dt,start,end,snow=F,elev=0,elev_band=NULL,plot=F){
  if(snow) return(CemaNeigeGR4jcab(dt=dt,start=start,end=end,elev=elev,elev_band=elev_band,plot=plot))
  else return(gr4jcab(dt=dt,start=start,end=end,plot=plot))
}

#' Combine simulation of GR4j & CemaNeige-GR4j
grsim=function(dt,start,end,param,snow=F,elev=0,elev_band=NULL,plot=F){
  if(snow) return(CemaNeigeGR4jsim(dt=dt,start=start,end=end,param=param,elev=elev,elev_band=elev_band,plot=plot))
  else return(gr4jsim(dt=dt,start=start,end=end,param=param,plot=plot))
}


#' To sum up snow-related variables of all layers
#' @param output the 'OutputsModel' object obtained by the calibration
#' @param elev_band the vector of min, q01 to q99 and max of catchment elevation distribution
#' @return a data.table of snowmelt, liquid precipitation, snowpack, production storage, and simulated runoff
SumOutputGr4j=function(output,elev_band=NULL){
  if(is.null(output$CemaNeigeLayers)) return(data.table(Prod=output$Prod,Qsim=output$Qsim))
  if(is.null(elev_band)){
    Melt=output$CemaNeigeLayers$Layer01$Melt
    Pliq=output$CemaNeigeLayers$Layer01$Pliq
    SnowPack=output$CemaNeigeLayers$Layer01$SnowPack
  }else{
    Melt=(output$CemaNeigeLayers$Layer01$Melt+output$CemaNeigeLayers$Layer02$Melt+
            output$CemaNeigeLayers$Layer03$Melt+output$CemaNeigeLayers$Layer04$Melt+
            output$CemaNeigeLayers$Layer05$Melt)/5
    Pliq=(output$CemaNeigeLayers$Layer01$Pliq+output$CemaNeigeLayers$Layer02$Pliq+
            output$CemaNeigeLayers$Layer03$Pliq+output$CemaNeigeLayers$Layer04$Pliq+
            output$CemaNeigeLayers$Layer05$Pliq)/5
    SnowPack=(output$CemaNeigeLayers$Layer01$SnowPack+output$CemaNeigeLayers$Layer02$SnowPack+
            output$CemaNeigeLayers$Layer03$SnowPack+output$CemaNeigeLayers$Layer04$SnowPack+
            output$CemaNeigeLayers$Layer05$SnowPack)/5
  }
  return(data.table(Melt=Melt,Pliq=Pliq,SnowPack=SnowPack,Prod=output$Prod,Qsim=output$Qsim))
}

