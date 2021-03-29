#' Second step in ICP-MS R analysis for mercury_202 using iridium_193 as internal standard.
#'
#' Arranges mercury (Hg_202) data to prepare for analysis by analyse_mercury3 or analyse_mercury_withdrift3
#' Enter Hg_202 e.g. setup_mercury2("Hg_202", Intensities)
#'
#' @param a Enter "Hg_202"
#' @param b Enter Intensities
#'
#' @return Dataframe with mercury data including RSD and internal standard
#' @export

setup_mercury2 <- function(a, b){

  seq_pos <- b$Sequence_position
  seq_pos <- seq_pos[!is.na(seq_pos)]

  Mercury_analysis <<- data.frame()
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
                "Hg_202 (cps)" = b$Hg_202[idx],
                "RSD_percent" = df_rsd,
                "Ir_193 (cps)" = b$Ir_193[idx])

    Mercury_analysis <<- rbind(Mercury_analysis, df)
  }
  View(Mercury_analysis)
}

