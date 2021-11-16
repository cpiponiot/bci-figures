#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Functions used for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


agb_bci <- function(dbh,
                    wd,
                    method = "chave14",
                    use_height_allom = FALSE) {
  # with generic height allometry from Martinez Cano et al 2019
  if (use_height_allom) {
    # Chave et al 2005 - moist forests, with height (in Mg)
    h <- 58.0 * dbh ^ 0.73 / (21.8 + dbh ^ 0.73)
    if (method == "chave05") {
      agb <- 0.0509 * wd * dbh ^ 2 * h / 1000
    }
    if (method == "chave14") {
      # Chave et al. 2014, equation 4 with the BIOMASS package 
      agb <- BIOMASS::computeAGB(D = dbh, WD = wd, H = h)
    }
  } else {
    # without any height information
    # Chave et al 2005 - moist forests, without height (in Mg)
    if (method == "chave05") {
      agb <- wd * exp(-1.499 + 2.148 * log(dbh) + 
                        0.207 * log(dbh) ^ 2 - 0.0281 * log(dbh) ^ 3) /1000
    }
    if (method == "chave14") {
      # Chave et al. 2014, equation 7 with the BIOMASS package (transform into kg)
      agb <- BIOMASS::computeAGB(D = dbh, WD = wd, coord = c(-79.8461, 9.1543))
    }
  }
  
  return(agb)
}
