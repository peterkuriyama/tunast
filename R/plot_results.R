#' Function to plot VAST results
#' Specify model characteristics outside of the function. Currently has the ability
#' to save the output and write the results to the directory

#' @param dat_index Specify the dataset to use numerically. Options are 1-annual counts; 
#' 2-annual numbers by size; 3 - annual numbers by size and season
#' @param ObsModel ObsModel configuration
#' @export


plot_results <- function(dat_index, ObsModel){

  browser
  DateFile = paste0(getwd(),"/VAST_output_", paste0("Dat_", dat_index,"_Obs_", 
    paste(ObsModel, collapse = "_"), "/"))
}