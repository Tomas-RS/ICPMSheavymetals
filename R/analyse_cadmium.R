#' Third and last step in ICP-MS R analysis for cadmium_111 without adjusting for drift.
#'
#' Transforms cadmium (Cd_111) cps outputs to parts per billion, accounting for In_115 internal standard,
#' plotting a graph for baseline drift using 1 ppb standards but not adjusting values according to drift,
#' using blank values to deduct background signal, plotting a cadmium standard curve,
#' accounting for 1/50 dilution during sample prep, and summarising cadmium results.
#' Results output includes a cadmium_results dataframe
#' and a .csv file in a new 'Results' folder in the R working directory
#'
#' Enter Cadmium_analysis, then names of standard for drift and standards for standard curve
#' e.g. analyse_cadmium3(Cadmium_analyis, "1 ppb standard", "0.0004 ppb metal standard", "0.002 ppb metal standard",
#' "0.02 ppb metal standard", "0.1 ppb metal standard", "0.5 ppb metal standard", "1 ppb metal standard",
#' "5 ppb metal standard", "10 ppb metal standard")
#'
#' Using this exact notation ("___ ppb metal standard") for standard curve samples is critical
#'
#' @param a Enter Cadmium_analysis
#' @param b Name of standards used for monitoring drift e.g. "1 ppb standard"
#' @param c Name of first standard curve sample e.g. "0.0004 ppb metal standard"
#' @param d Name of second standard curve sample e.g. "0.002 ppb metal standard"
#' @param e Name of third standard curve sample e.g. "0.02 ppb metal standard"
#' @param f Name of fourth standard curve sample e.g. "0.1 ppb metal standard"
#' @param g Name of fifth standard curve sample e.g. "0.5 ppb metal standard"
#' @param h Name of sixth standard curve sample e.g. "1 ppb metal standard"
#' @param i Name of seventh standard curve sample e.g. "5 ppb metal standard"
#' @param j Name of eighth standard curve sample e.g. "10 ppb metal standard"
#'
#' @return Dataframe and .csv file with all cadmium results
#' @export

analyse_cadmium3 <- function(a,b,c,d,e,f,g,h,i,j){
  a$Z01_Infactor <- as.numeric(a$In_115) / 200000
  a$Z02_Inadjust <- as.numeric(a$Cd_111) / a$Z01_Infactor

  Cadmium_standards <- filter(a, Sample_name == b)
  trendline(Cadmium_standards$Sequence_position, Cadmium_standards$Z02_Inadjust,model = "line2P",
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3,
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift of standards",
            xlab = "Sequence position", ylab = "Signal (Cd)", CI.fill = FALSE, CI.lty = "blank")


  Cadmium_blanks <- filter(a, Sample_name == "Blank")
  plot(Cadmium_blanks$Sequence_position, Cadmium_blanks$Z02_Inadjust,
       main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Cd)")

  a$Z04_Cd_blankadjust <- a$Z02_Inadjust - mean(Cadmium_blanks[["Z02_Inadjust"]])

  Cadmium_standardcurve <- filter(a, Sample_name == c
                                  | Sample_name == d | Sample_name == e
                                  | Sample_name == f | Sample_name == g
                                  | Sample_name == h | Sample_name == i
                                  | Sample_name == j)
  Cadmium_stdconccharacters <- gsub(" ppb metal standard", "", Cadmium_standardcurve$Sample_name)
  Cadmium_standardcurve$Concentration <- as.numeric(Cadmium_stdconccharacters)
  trendline(Cadmium_standardcurve$Concentration, Cadmium_standardcurve$Z04_Cd_blankadjust,model = "line2P",
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3,
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main="Cd standard curve",
            xlab = "Cd concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
  Cadmium_stdplot <- lm(Cadmium_standardcurve$Z04_Cd_blankadjust ~ Cadmium_standardcurve$Concentration)
  a$Z05_Cd_ppb <- a$Z04_Cd_blankadjust / coef(Cadmium_stdplot)[2]


  a$Z07_Sample_Cd_ppb <- a$Z05_Cd_ppb * 50

  a$Z08Cd_RSD_ppb <- a$Z07_Sample_Cd_ppb * (as.numeric(a$RSD_percent) / 100)

  Cadmium_detectionlimit_ppb <- (3 * sd(Cadmium_blanks[["Z02_Inadjust"]])) / coef(Cadmium_stdplot)[2]

  Cadmium_results <- select(a, Sample_name, Z07_Sample_Cd_ppb, RSD_percent, Z08Cd_RSD_ppb, Z05_Cd_ppb)

  Cadmium_results <<- data.frame(Sample_name = Cadmium_results$Sample_name,
                                Cd_ppb_original_sample = Cadmium_results$Z07_Sample_Cd_ppb,
                                RSD_percent = Cadmium_results$RSD_percent,
                                RSD_ppb = Cadmium_results$Z08Cd_RSD_ppb,
                                Cd_ppb_detected_in_measured_sample = Cadmium_results$Z05_Cd_ppb,
                                Detection_limit_ppb = Cadmium_detectionlimit_ppb)
  Cadmium_analysis <<- a
  View(Cadmium_analysis)
  View(Cadmium_results)

  if(!dir.exists("Results")) dir.create("Results")
  write.csv(Cadmium_results, file = file.path("Results", "Cadmium_results.csv"), row.names = FALSE)
}
