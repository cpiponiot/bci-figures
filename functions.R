#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Functions used for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# estimate agb values at BCI with allometric equations from Chave et al. 2014
# (method = "chave14") or Chave et al., 2005 (method = "chave05"); the height
# equation from Martinez Cano et al 2019 is used when use_height_allom = TRUE
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


# function to interpolate missing dbh value ####

# param dbh: all dbhs measurements of one tree or stem
# param year: the years corresponding to those measurements
# param DFstatus: the status of the tree or stem for all those measurements

# returns a dbh vector with all missing values interpolated (only when there is
# at least one dbh measurement before and after the missing dbh)

interpolate_missing <- function(dbh, year, DFstatus) {
  if (any(DFstatus == "missing" & !is.na(DFstatus))) { 
    # order dbh and year vector in chronological order
    year_ord <- year[order(year)]
    dbh_ord <- dbh[order(year)]
    
    # years of all missing values 
    for (y0 in year[DFstatus == "missing" & !is.na(DFstatus)]) {
      # last (non NA) dbh measurement prior to the missing measurement
      dbh1 <- data.table::last(dbh_ord[!is.na(dbh_ord) & year_ord < y0])
      # first (non NA) dbh measurement after the missing measurement
      dbh2 <- data.table::first(dbh_ord[!is.na(dbh_ord) & year_ord > y0])
      
      # check that dbh1 and dbh2 have values (otherwise: no data to interpolate from)
      if (length(dbh1) > 0 & length(dbh2) > 0) {
        # last year with (non NA) dbh measurement prior to the missing measurement
        year1 <- max(year_ord[!is.na(dbh_ord) & year_ord < y0])
        # first year with (non NA) dbh measurement after the missing measurement
        year2 <- min(year_ord[!is.na(dbh_ord) & year_ord > y0])
        
        # slope of the regression between the measurements before and after the
        # missing value
        slope <- (dbh2 - dbh1) / (year2 - year1)
        # interpolation:
        dbh_ord[year_ord == y0] <- slope * (y0 - year1) + dbh1
      } 
    }
    
    return(dbh_ord[order(order(year))])
    
  } else {
    return(dbh)
  }
}
# Kohyama correction function ####
# instanteneous biomass (or other) fluxes as recommended by Kohyama 2019 (eq 1-2 in Table 1)

kohyama_correction <- function(dt, vars) {
  dt_sub = dt[variable %in% vars, .(value = sum(value * weight) / sum(weight),
                                    weight = sum(weight)),
              .(dT, year, variable, group, site)]
  dt_new = dcast(dt_sub, year + dT + site + weight + group ~ variable)
  dt_new$B0 = dt_new[, vars[1], with = FALSE]
  dt_new$gain = dt_new[, vars[2], with = FALSE]
  dt_new$loss = dt_new[, vars[3], with = FALSE]
  dt_new[, BS0 := B0 - loss * dT]
  dt_new[, BT := B0 + (gain - loss) * dT]
  dt_new[, vars[2] := (log(BT / BS0) * (BT - B0)) / (dT * log(BT / B0))]
  dt_new[, vars[3] := (log(B0 / BS0) * (BT - B0)) / (dT * log(BT / B0))]
  
  dt_new = melt(dt_new,
                id.vars = c("site", "weight", "group"),
                measure.vars = vars)
  dt_corr = rbind(dt[!variable %in% vars, colnames(dt_new), with = FALSE], dt_new)
  dt_corr = subset(dt_corr,!is.na(value))
  return(dt_corr)
}

