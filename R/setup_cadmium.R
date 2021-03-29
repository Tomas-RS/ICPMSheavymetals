#' Second step in ICP-MS R analysis for cadmium_111 using indium_115 as internal standard.
#'
#' Arranges cadmium (Cd_111) data to prepare for analysis by analyse_cadmium3 or analyse_cadmium_withdrift3
#' Enter Cd_111 e.g. setup_cadmium2("Cd_111", Intensities)
#'
#' @param a Enter "Cd_111"
#' @param b Enter Intensities
#'
#' @return Dataframe with cadmium data including RSD and internal standard
#' @export

setup_cadmium2 <- function(a, b){

seq_pos <- b$Sequence_position
seq_pos <- seq_pos[!is.na(seq_pos)]

Cadmium_analysis <<- data.frame()
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
              "Cd_111 (cps)" = b$Cd_111[idx],
              "RSD_percent" = df_rsd,
              "In_115 (cps)" = b$In_115[idx])

  Cadmium_analysis <<- rbind(Cadmium_analysis, df)
}
View(Cadmium_analysis)
}

