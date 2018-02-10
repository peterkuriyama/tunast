#-----------------------------------------------------------------------------------------------------
#Set working directory based on computer I'm using
setwd("C:\\Users\\Peter\\Dropbox\\postdoc\\tunast")
if(Sys.info()[1] == 'Windows') setwd('\\Users\\peter.kuriyama\\Desktop\\tunast')
if(Sys.info()[1] != 'Windows') setwd('/Users/peterkuriyama/Dropbox/postdoc/tunast')

#-----------------------------------------------------------------------------------------------------
#Install packages
# install.packages("TMB")

# devtools::install_github("james-thorson/VAST",force = TRUE)
# devtools::install_github("james-thorson/utilities",force = TRUE)
# devtools::install_github("nwfsc-assessz/geostatistical_delta-GLMM",force = TRUE)

#Load Packages
library(TMB)
library(plyr)
library(VAST)
library(tidyverse)
library(SpatialDeltaGLMM)
library(reshape2)
library(devtools)
source("Grid_Fn.R")
source("Spatial_Information_Fn.R")
source("Calc_Polygon_Areas_and_Polygons_Fn.R")
source("Calc_Anisotropic_Mesh.R")

devtools::install_github('peterkuriyama/tunast')
library(tunast)
#-----------------------------------------------------------------------------------------------------
#Load the Data
the_data <- process_bigeye_data(nspp = 3)

#Check the proportion zeroes in the data
# the_data[[2]] %>% group_by(spp, year, lat, lon) %>% summarize(nzeroes = length(which(cpue == 0))) %>% 
#   ungroup %>% distinct(nzeroes)
# hist(the_data[[2]]$cpue, breaks = 50)

# head(the_data[[2]])

#-----------------------------------------------------------------------------------------------------
#VAST Model

#Link for handling zero or 100 percent encounter probability
# https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/What-to-do-with-a-species-with-0%25-or-100%25-encounters-in-any-year

#---------------------------------
#Specify VAST Model Configurations
Version <- "VAST_v4_0_0"
Method <- c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km <- 10
n_x <- c(10, 50, 100, 200, 1000)[2] #number of stations
Kmeans_Config <-  list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )

#Specify number of species
nspp <- 3
#------------FieldConfig
#Omega refers to spatial variation
  #Omega1 - variation in encounter probability
  #Omega2 - variation in positive catch rates; 0 is off
#Epsilon is spatiotemporal variation
#AR1 is an AR1 process
#>0 is number of elements in factor-analysis covariance
FieldConfig <- c("Omega1"= nspp, "Epsilon1"=nspp, "Omega2"=nspp, "Epsilon2"=nspp) #

#------------Rhoconfig
#Specify structure of time intervals (fixed effect or random effect)
# 1: 
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) #Parameters among years

#------------Overdispersion Config
#Optional, vector governing correlated overdispersion among categories v_i (vessels)
OverdispersionConfig <- c("Vessel"=0, "VesselYear"=0)  #Different vessel catchabilities, 0 means they the same

#------------ObsModel
#First element specifies distribution of positive catch rates
  #0 - normal
  #1 - lognormal
  #2 - Gamma
  #5 - negative binomial
  #7 - poisson
#Second element specifies functional form of encounter probabilities
  #0 - Conventional delta-model, logit-link encounter and log-link positive catch rates
  #1 - Alternative delta-model, log-link for numbers-density, log-link biomass per number
  #2 - Link for tweedie observation
  #3 - Conventional delta-model, fixes encounter probability to 1
ObsModel <- c(1, 0) 

#Specify options
Options <- c("SD_site_density"=0, "SD_site_logdensity"=0, "Calculate_Range"=1,
             "Calculate_evenness"=0, "Calculate_effective_area"=1, "Calculate_Cov_SE"=0,
             'Calculate_Synchrony'=0, 'Calculate_Coherence'=0)

strata.limits <- data.frame('STRATA'="EPO")
Region <- "Other"

#Look at proportion of zeroes
# the_data$bet_numbers_annual %>% group_by(year, lat, lon) 

#-----------------------------------------------------------------------------------------------------
#Run example with annual counts and one species
nspp <- 1
FieldConfig <- c("Omega1"= nspp, "Epsilon1"=nspp, "Omega2"=nspp, "Epsilon2"=nspp) #
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) #Parameters among years

#Run the model
run_vast(dat_index = 1, ObsModel = c(1, 0))

#analyze the results
plot_results(dat_index = 1, ObsModel = c(1, 0))


#Run the length comp example
nspp <- 3
FieldConfig <- c("Omega1"= nspp, "Epsilon1"=nspp, "Omega2"=nspp, "Epsilon2"=nspp) #
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) #Parameters among years
run_vast(dat_index = 2, ObsModel = c(1, 3))


#---------------------------------



# For function: Specify data and ObsModel



# Enc_prob = SpatialDeltaGLMM::Check_encounter_prob( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)
Q = SpatialDeltaGLMM::QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
MapDetails_List[["Cex"]] = 2
# Decide which years to plot
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

residuals <- SpatialDeltaGLMM:::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

# SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )

# # pre <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
# # # pos <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=2, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
density_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
presence_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
positive_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=2, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
# # e_presence_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=6, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
# # e_positive_xt <- SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=7, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
# #
Index = SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, strata_names="all_areas", use_biascorr=TRUE )
# #
SpatialDeltaGLMM::Plot_range_shifts(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=DateFile, Year_Set=Year_Set)
# #
Plot_factors(Report = Report, ParHat = Obj$env$parList(), Data = TmbData, SD = Opt$SD, mapdetails_list = MapDetails_List,
             Year_Set = Year_Set, category_names = c("L1","L2","L3","L4","L5"), plotdir = DateFile)

Cov_List = Summarize_Covariance(Report = Report, ParHat = Obj$env$parList(), Data = TmbData, SD = Opt$SD, plot_cor = TRUE,
                                category_names = levels(Data_Geostat[, "spp"]), plotdir = DateFile, plotTF = FieldConfig,
                                mgp = c(2, 0.5, 0), tck = -0.02, oma = c(0, 5, 2, 2))






