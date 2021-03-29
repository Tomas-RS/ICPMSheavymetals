#' First step in ICP-MS R analysis
#'
#' Installs necessary packages and loads ICP-MS result.
#' Enter name of saved ICP-MS results file e.g. setup_all1("Sample_run_HeavyMetalstest_modified.csv")
#'
#' @param x name of ICP-MS results file
#'
#' @return Loaded packages and a ICP-MS results file made into the "intensity" dataframe
#' @export

setup_all1 <- function(x){
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!require("basicTrendline")) install.packages("basicTrendline")
library(basicTrendline)

Intensities <<- read.csv(x)
View(Intensities)

paste("The setup_all1 function has installed/loaded tidyverse and basicTrendline
      and loaded your ICPMS data analysis file as a dataframe named 'Intensities'.
      Prior to using this you should have set your working directory
      and prepared the data analysis file by adding Sequence_position and
      Sequence_name column titles, underscores to metals (Cd 111 becomes Cd_111),
      removing columns for acquisition time, dataset file and method file,
      and deleting all worksheets except the one with both Intensities and RSDs.
      Save this worksheet as a .csv file (Sample_run_XYZ) in same folder as working directory/project")
}
