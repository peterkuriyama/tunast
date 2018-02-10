#' Function run VAST
#' Specify model characteristics outside of the function. Currently has the ability
#' to save the output and write the results to the directory

#' @param dat_index Specify the dataset to use numerically. Options are 1-annual counts; 
#' 2-annual numbers by size; 3 - annual numbers by size and season
#' @param ObsModel ObsModel configuration
#' @export


run_vast <- function(dat_index, ObsModel){

  #Specify data set
  dat <- the_data[[dat_index]]  
  
  DateFile = paste0(getwd(),"/VAST_output_", paste0("Dat_", dat_index,"_Obs_", 
    paste(ObsModel, collapse = "_"), "/"))
  dir.create(DateFile)

  #Save configurations
  Record = ThorsonUtilities::bundlelist( c("Version","Method","grid_size_km","n_x","FieldConfig","RhoConfig","OverdispersionConfig","ObsModel","Kmeans_Config") )
  save( Record, file=file.path(DateFile,"Record.RData"))
  capture.output( Record, file=paste0(DateFile,"Record.txt"))

  #Specify the geostatistical data  
  Data_Geostat <- data.frame(spp = as.character(dat$spp), 
                            Catch_num = dat$cpue,
                            Year = dat$year, 
                            Vessel = "missing", 
                            Lat = dat$lat,
                            AreaSwept_km2 = 1, 
                            Lon = dat$lon)
  
  #Specify extrapolation list
  Extrapolation_List <- Grid_Fn(Data_Geostat = Data_Geostat,
    zone = 12, grid_dim_km = 10)
  save(Extrapolation_List, file = paste0(DateFile, "Extrapolation_List.RData"))

  #Extrapolate the grids  
  Spatial_List = Spatial_Information_Fn( grid_size_km = grid_size_km, n_x = n_x, 
    Method = Method, Lon = Data_Geostat[,'Lon'], Lat = Data_Geostat[,'Lat'], 
    Extrapolation_List = Extrapolation_List, 
    randomseed = Kmeans_Config[["randomseed"]], nstart = Kmeans_Config[["nstart"]],  
   iter.max = Kmeans_Config[["iter.max"]], DirPath = DateFile, Save_Results = TRUE)
  save(Spatial_List, file = paste0(DateFile, "Spatial_List.RData"))
  
  # Add knots to Data_Geostat
  Data_Geostat = cbind(Data_Geostat, "knot_i"=Spatial_List$knot_i)
  save(Data_Geostat, file = paste0(DateFile, "Data_Geostat.RData"))
  SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List = Extrapolation_List,
   Spatial_List = Spatial_List, Data_Geostat = Data_Geostat, PlotDir = DateFile)

  #---------------------------------
  #Run the model
  TmbData = VAST::Data_Fn("Version" = Version, "FieldConfig" = FieldConfig, 
    "OverdispersionConfig" = OverdispersionConfig,  
    "RhoConfig" = RhoConfig, "ObsModel" = ObsModel, 
    "c_i" = as.numeric(Data_Geostat$spp) - 1, 
    "b_i" = Data_Geostat$Catch_num, 
    "a_i" = Data_Geostat$AreaSwept_km2, 
    "v_i" = as.numeric(Data_Geostat$Vessel) - 1, 
    "s_i" = Data_Geostat$knot_i - 1, 
    "t_i" = Data_Geostat$Year, 
    "a_xl" = Spatial_List$a_xl, 
    "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, 
    "Method"=Spatial_List$Method, "Options"=Options )
  
  TmbList = VAST::Build_TMB_Fn("TmbData" = TmbData, 
    "RunDir" = DateFile, 
    "Version" = Version, 
    "RhoConfig" = RhoConfig, 
    "loc_x" = Spatial_List$loc_x, "Method" = Method)
  save(TmbList, file = paste0(DateFile, "TmbList.RData"))
  Obj = TmbList[["Obj"]]
  
  Opt = TMBhelper::Optimize(obj = Obj, 
    newtonsteps = 1,
    lower = TmbList[["Lower"]], upper = TmbList[["Upper"]], getsd = TRUE, 
    savedir = DateFile, bias.correct = FALSE)

  Report = Obj$report()
  Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
  save(Save, file=paste0(DateFile,"Save.RData"))
  save.image(paste0(DateFile,"Result.RData"))
}
