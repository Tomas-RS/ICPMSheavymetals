#' Second step in ICP-MS R analysis for lead_208 using iridium_193 as internal standard.
#'
#' Arranges lead (Pb_208) data to prepare for analysis by analyse_lead3 or analyse_lead_withdrift3
#' Enter Cd_111 e.g. setup_lead2("Pb_208", Intensities)
#'
#' @param a Enter "Pb_208"
#' @param b Enter Intensities
#'
#' @return Dataframe with lead data including RSD and internal standard
#' @export

setup_lead2 <- function(a, b){

  seq_pos <- b$Sequence_position
  seq_pos <- seq_pos[!is.na(seq_pos)]

  Lead_analysis <<- data.frame()
  for (n in seq_pos) {
    idx <- which(b$Sequence_position == n)
    df_rsd <- b[idx+1, a]
    df_rsd <- lapply(df_rsd, function(x) {
      gsub("%", "", x)
    })
    df_rsd <- as.data.frame(df_rsd, check.names = FALSE)
    colnames(df_rsd) <- "RSD_percent"
    df <- cbind(data.frame(Sequence_position = n,
                           Sample_name = b[idx,"Sample_name"]),
                "Pb_208 (cps)" = b$Pb_208[idx],
                "RSD_percent" = df_rsd,
                "Ir_193 (cps)" = b$Ir_193[idx])

    Lead_analysis <<- rbind(Lead_analysis, df)
  }
  View(Lead_analysis)
}

