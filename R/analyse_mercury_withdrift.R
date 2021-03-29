#' Third and last step in ICP-MS R analysis for mercury_202 when adjusting for drift.
#'
#' Transforms mercury (Hg_202) cps outputs to parts per billion, accounting for Ir_193 internal standard,
#' plotting a graph for baseline drift using 1 ppb standards and adjusting values according to drift,
#' using blank values to deduct background signal, plotting a mercury standard curve,
#' accounting for 1/50 dilution during sample prep, and summarising mercury results.
#' Results output includes a mercury_results dataframe
#' and a .csv file in a new 'Results' folder in the R working directory
#'
#' Enter Mercury_analysis, then names of standard for drift and standards for standard curve
#' e.g. analyse_mercury_withdrift3(Mercury_analyis, "1 ppb standard", "0.0004 ppb metal standard", "0.002 ppb metal standard",
#' "0.02 ppb metal standard", "0.1 ppb metal standard", "0.5 ppb metal standard", "1 ppb metal standard",
#' "5 ppb metal standard", "10 ppb metal standard")
#'
#' Using this exact notation ("___ ppb metal standard") for standard curve samples is critical
#'
#' @param a Enter Mercury_analysis
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
#' @return Dataframe and .csv file with all mercury results
#' @export

analyse_mercury_withdrift3 <- function(a,b,c,d,e,f,g,h,i,j){
  a$Z01_Irfactor <- as.numeric(a$Ir_193) / 200000
  a$Z02_Iradjust <- as.numeric(a$Hg_202) / a$Z01_Irfactor

  Mercury_standards <- filter(a, Sample_name == b)
  trendline(Mercury_standards$Sequence_position, Mercury_standards$Z02_Iradjust,model = "line2P",
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3,
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift of standards",
            xlab = "Sequence position", ylab = "Signal (Hg)", CI.fill = FALSE, CI.lty = "blank")


  Mercury_driftplot <- lm(Mercury_standards$Z02_Iradjust ~ Mercury_standards$Sequence_position)
  a$Z03_Hg_driftadjust <- a$Z02_Iradjust - coef(Mercury_driftplot)[2] *
    (a$Sequence_position - Mercury_standards[1,1])

  Mercury_blanks <- filter(a, Sample_name == "Blank")
  plot(Mercury_blanks$Sequence_position, Mercury_blanks$Z03_Hg_driftadjust,
       main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Hg)")

  a$Z04_Hg_blankadjust <- a$Z03_Hg_driftadjust - mean(Mercury_blanks[["Z03_Hg_driftadjust"]])

  Mercury_standardcurve <- filter(a, Sample_name == c
                                  | Sample_name == d | Sample_name == e
                                  | Sample_name == f | Sample_name == g
                                  | Sample_name == h | Sample_name == i
                                  | Sample_name == j)
  Mercury_stdconccharacters <- gsub(" ppb metal standard", "", Mercury_standardcurve$Sample_name)
  Mercury_standardcurve$Concentration <- as.numeric(Mercury_stdconccharacters)
  trendline(Mercury_standardcurve$Concentration, Mercury_standardcurve$Z04_Hg_blankadjust,model = "line2P",
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3,
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main="Hg standard curve",
            xlab = "Hg concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
  Mercury_stdplot <- lm(Mercurystandardcurve$Z04_Hg_blankadjust ~ Mercury_standardcurve$Concentration)
  a$Z05_Hg_ppb <- a$Z04_Hg_blankadjust / coef(Mercury_stdplot)[2]


  a$Z07_Sample_Hg_ppb <- a$Z05_Hg_ppb * 50

  a$Z08Hg_RSD_ppb <- a$Z07_Sample_Hg_ppb * (as.numeric(a$RSD_percent) / 100)

  Mercury_detectionlimit_ppb <- (3 * sd(Mercury_blanks[["Z02_Iradjust"]])) / coef(Mercury_stdplot)[2]

  Mercury_results <- select(a, Sample_name, Z07_Sample_Hg_ppb, RSD_percent, Z08Hg_RSD_ppb, Z05_Hg_ppb)

  Mercury_results <<- data.frame(Sample_name = Mercury_results$Sample_name,
                                 Hg_ppb_original_sample = Mercury_results$Z07_Sample_Hg_ppb,
                                 RSD_percent = Mercury_results$RSD_percent,
                                 RSD_ppb = Mercury_results$Z08Hg_RSD_ppb,
                                 Hg_ppb_detected_in_measured_sample = Mercury_results$Z05_Hg_ppb,
                                 Detection_limit_ppb = Mercury_detectionlimit_ppb)
  Mercury_analysis <<- a
  View(Mercury_analysis)
  View(Mercury_results)

  if(!dir.exists("Results")) dir.create("Results")
  write.csv(Mercury_results, file = file.path("Results", "Mercury_results.csv"), row.names = FALSE)
}
