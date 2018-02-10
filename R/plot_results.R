#' Function to plot VAST results
#' Specify model characteristics outside of the function. Currently has the ability
#' to save the output and write the results to the directory

#' @param dat_index Specify the dataset to use numerically. Options are 1-annual counts; 
#' 2-annual numbers by size; 3 - annual numbers by size and season
#' @param ObsModel ObsModel configuration
#' @export


plot_results <- function(dat_index, ObsModel){


  #------------Load and parse the model results
  DateFile = paste0(getwd(),"/VAST_output_", paste0("Dat_", dat_index,"_Obs_", 
    paste(ObsModel, collapse = "_"), "/"))
  load(paste0(DateFile, "Save.Rdata"))
  Opt <- Save$Opt
  Report <- Save$Report
  ParHat <- Save$ParHat
  TmbData <- Save$TmbData
  rm(Save) #remove the data

  #Load Data types
  list.files(DateFile)
  load(paste0(DateFile, "Extrapolation_List.Rdata"))
  load(paste0(DateFile, "Spatial_List.Rdata"))
  load(paste0(DateFile, "Data_Geostat.Rdata"))
  #Encounter probability
  # Enc_prob = SpatialDeltaGLMM::Check_encounter_prob( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)

  #------------Process the QQ plot
  Q = SpatialDeltaGLMM::QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))
  print("QQ plots done")
  #------------Plot the extrapolation Map
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
  MapDetails_List[["Cex"]] = 2
  print("Maps done")
  #------------Decide which years to plot
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  
  #------------Residuals
  residuals <- SpatialDeltaGLMM:::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)
  print("residuals done")
  #------------Option to plot anisotropy
  # SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )
  
  #------------Plot spatial results
  # # pre <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  # # # pos <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=2, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  density_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  presence_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  positive_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=2, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  print("Maps done")
  # # e_presence_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=6, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  # # e_positive_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=7, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
  
  #------------Index of abundance
  Index = SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, strata_names="all_areas", use_biascorr=TRUE )
  print("Index done")

  #------------Range Shifts
  SpatialDeltaGLMM::Plot_range_shifts(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=DateFile, Year_Set=Year_Set)
  
  #------------Factors
  Plot_factors(Report = Report, ParHat = Obj$env$parList(), Data = TmbData, SD = Opt$SD, mapdetails_list = MapDetails_List,
               Year_Set = Year_Set, category_names = c("L1","L2","L3","L4","L5"), plotdir = DateFile)
  print("Factors done")

  #------------Covariance List
  Cov_List = Summarize_Covariance(Report = Report, ParHat = Obj$env$parList(), Data = TmbData, SD = Opt$SD, plot_cor = TRUE,
                                  category_names = levels(Data_Geostat[, "spp"]), plotdir = DateFile, plotTF = FieldConfig,
                                  mgp = c(2, 0.5, 0), tck = -0.02, oma = c(0, 5, 2, 2))
  print("covariance done")

}