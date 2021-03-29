#' Third and last step in ICP-MS R analysis for lead_208 when adjusting for drift.
#'
#' Transforms lead (Pb_208) cps outputs to parts per billion, accounting for Ir_193 internal standard,
#' plotting a graph for baseline drift using 1 ppb standards and adjusting values according to drift,
#' using blank values to deduct background signal, plotting a lead standard curve,
#' accounting for 1/50 dilution during sample prep, and summarising lead results.
#' Results output includes a lead_results dataframe
#' and a .csv file in a new 'Results' folder in the R working directory
#'
#' Enter Lead_analysis, then names of standard for drift and standards for standard curve
#' e.g. analyse_lead_withdrift3(Lead_analyis, "1 ppb standard", "0.0004 ppb metal standard", "0.002 ppb metal standard",
#' "0.02 ppb metal standard", "0.1 ppb metal standard", "0.5 ppb metal standard", "1 ppb metal standard",
#' "5 ppb metal standard", "10 ppb metal standard")
#'
#' Using this exact notation ("___ ppb metal standard") for standard curve samples is critical
#'
#' @param a Enter Lead_analysis
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
#' @return Dataframe and .csv file with all lead results
#' @export

analyse_lead_withdrift3 <- function(a,b,c,d,e,f,g,h,i,j){
  a$Z01_Irfactor <- as.numeric(a$Ir_193) / 200000
  a$Z02_Iradjust <- as.numeric(a$Pb_208) / a$Z01_Irfactor

  Lead_standards <- filter(a, Sample_name == b)
  trendline(Lead_standards$Sequence_position, Lead_standards$Z02_Iradjust,model = "line2P",
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3,
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main="Baseline Drift of standards",
            xlab = "Sequence position", ylab = "Signal (Pb)", CI.fill = FALSE, CI.lty = "blank")

  Lead_driftplot <- lm(Lead_standards$Z02_Iradjust ~ Lead_standards$Sequence_position)
  a$Z03_Pb_driftadjust <- a$Z02_Iradjust - coef(Lead_driftplot)[2] *
    (a$Sequence_position - Lead_standards[1,1])

  Lead_blanks <- filter(a, Sample_name == "Blank")
  plot(Lead_blanks$Sequence_position, Lead_blanks$Z03_Pb_driftadjust,
       main="Trend of blank values", xlab = "Sequence position", ylab = "Signal (Pb)")

  a$Z04_Pb_blankadjust <- a$Z03_Pb_driftadjust - mean(Lead_blanks[["Z03_Pb_driftadjust"]])

  Lead_standardcurve <- filter(a, Sample_name == c
                                  | Sample_name == d | Sample_name == e
                                  | Sample_name == f | Sample_name == g
                                  | Sample_name == h | Sample_name == i
                                  | Sample_name == j)
  Lead_stdconccharacters <- gsub(" ppb metal standard", "", Lead_standardcurve$Sample_name)
  Lead_standardcurve$Concentration <- as.numeric(Lead_stdconccharacters)
  trendline(Lead_standardcurve$Concentration, Lead_standardcurve$Z04_Pb_blankadjust,model = "line2P",
            Pvalue.corrected = TRUE, linecolor = "blue", lty = 1, lwd = 3,
            show.equation = TRUE, show.Rpvalue = TRUE, Rname = 1, Pname = 0, xname = "x", yname = "y",
            yhat = FALSE, summary = TRUE, text.col = "black", main="Pb standard curve",
            xlab = "Pb concentration (ppb)", ylab = "Signal", CI.fill = FALSE, CI.lty = "blank")
  Lead_stdplot <- lm(Leadstandardcurve$Z04_Pb_blankadjust ~ Lead_standardcurve$Concentration)
  a$Z05_Pb_ppb <- a$Z04_Pb_blankadjust / coef(Lead_stdplot)[2]


  a$Z07_Sample_Pb_ppb <- a$Z05_Pb_ppb * 50

  a$Z08Pb_RSD_ppb <- a$Z07_Sample_Pb_ppb * (as.numeric(a$RSD_percent) / 100)

  Lead_detectionlimit_ppb <- (3 * sd(Lead_blanks[["Z02_Iradjust"]])) / coef(Lead_stdplot)[2]

  Lead_results <- select(a, Sample_name, Z07_Sample_Pb_ppb, RSD_percent, Z08Pb_RSD_ppb, Z05_Pb_ppb)

  Lead_results <<- data.frame(Sample_name = Lead_results$Sample_name,
                                 Pb_ppb_original_sample = Lead_results$Z07_Sample_Pb_ppb,
                                 RSD_percent = Lead_results$RSD_percent,
                                 RSD_ppb = Lead_results$Z08Pb_RSD_ppb,
                                 Pb_ppb_detected_in_measured_sample = Lead_results$Z05_Pb_ppb,
                                 Detection_limit_ppb = Lead_detectionlimit_ppb)
  Lead_analysis <<- a
  View(Lead_analysis)
  View(Lead_results)

  if(!dir.exists("Results")) dir.create("Results")
  write.csv(Lead_results, file = file.path("Results", "Lead_results.csv"), row.names = FALSE)
}
